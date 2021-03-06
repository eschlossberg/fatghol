#! /usr/bin/env python
#
"""Classes and functions to deal with fatgraphs.
"""
#
#   Copyright (C) 2008-2012 Riccardo Murri <riccardo.murri@gmail.com>
#   All rights reserved.
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
__docformat__ = 'reStructuredText'

# import cython

## stdlib imports

from collections import defaultdict, Iterator
import logging
import os.path

## application-local imports

import fatghol
from fatghol.cache import (
    Caching,
    ocache0,
    ocache_contract,
    ocache_eq,
    ocache_isomorphisms,
)
from fatghol.combinatorics import Permutation
from fatghol.cyclicseq import CyclicList, CyclicTuple
from fatghol.iterators import (
    BufferingIterator,
)
from fatghol.loadsave import load, save
from fatghol.runtime import runtime
import fatghol.timing as timing
from fatghol.utils import (
    concat,
    ltranslate,
    maybe,
    sign,
)


## main

class BoundaryCycle(frozenset):
    """A boundary cycle of a Fatgraph.

    Boundary cycles are a cyclic sequence of 'corners': a corner
    consists of a vertex `v` and (an unordered pair of) two
    consecutive indices (in the cyclic order at `v`, so, either
    `j == i+1` or `i` and `j` are the starting and ending indices).

    Two boundary cycles are equal if they comprise the same
    corners.

    A `BoundaryCycle` instance is constructed from a sequence of
    triples `(v, i, j)`, where `i` and `j` are consecutive (in the
    cyclic order sense) indices at a vertex `v`.  (Although no check
    is performed in the constructor code.)
    """

    __slots__ = ()

    # no code to be added to the `frozenset` base class
    pass


class Vertex(CyclicList):
    """A (representative of) a vertex of a ribbon graph.

    A vertex is represented by the cyclically ordered list of its
    (decorated) edges.  The edge colorings may be accessed through a
    (read-only) sequence interface.

    At init time, the number of loops attached to this vertex is
    computed and stored in the `.num_loops` attribute::

      >>> Vertex([0,1,2]).num_loops
      0
      >>> Vertex([1,1,0,0]).num_loops
      2
    """
    # FIXME: reimplement `Vertex` as a `tuple`/`CyclicTuple` subclass.

    __slots__ = ('num_loops')

    def __init__(self, seq):
        CyclicList.__init__(self, seq)
        self.num_loops = len(self) - len(set(self))

    def __str__(self):
        return repr(self)


class Edge(object):
    """An edge of a fatgraph.

    An edge is represented by its two endpoints; each endpoint has the form
    `(v_idx, a_idx)`, where:
      - `v_idx` is the index of the endpoint vertex (within the fatgraph), and
      - `a_idx` is the index at which this edge appears within the vertex `v_idx`.

    It is guaranteed that endpoints are stored in increasing vertex
    index order::

      >>> x = Edge((3, 0), (1, 2))
      >>> x.endpoints
      ((1, 2), (3, 0))
      >>> x.endpoints[0][0] < x.endpoints[1][0]
      True

    `Edge` objects bear no reference to a particular `Fatgraph`
    instance, to allow sharing the same edge instance among
    `Fatgraph`s that are created by contraction or other geometrical
    operations.
    """

    __slots__ = ('endpoints')

    def __init__(self, va1, va2):
        if va1[0] < va2[0]:
            self.endpoints = (va1, va2)
        else:
            self.endpoints = (va2, va1)

    def is_loop(self):
        """Return `True` if this `Edge` instance represents looping edge.

        Examples::

          >>> Edge((0, 1), (0, 2)).is_loop()
          True
          >>> Edge((0, 1), (1, 3)).is_loop()
          False
        """
        return self.endpoints[0][0] == self.endpoints[1][0]

    def meets(self, v):
        """Return `True` if vertex `v` is one of the endpoints.

        Example::

          >>> x = Edge((3, 1), (1, 2))
          >>> x.meets(3)
          True
          >>> x.meets(1)
          True
          >>> x.meets(0)
          False          
        """
        return (v == self.endpoints[0][0]) or (v == self.endpoints[1][0])

    def other_end(self, v, a):
        """Return the endpoint opposed to `(v, a)`.

        Example::

          >>> l = Edge((0, 1), (1, 3))
          >>> l.other_end(0,1)
          (1, 3)
          >>> l.other_end(1,3)
          (0, 1)
        """
        if self.endpoints[0] == (v, a):
            return self.endpoints[1]
        else:
            return self.endpoints[0]


class Isomorphism(object):
    """An isomorphism of `Fatgraphs`.
    """

    __slots__ = (
        'pe',
        'pv',
        'rot',
        'source',
        'target',
    )

    def __init__(self, source, target, pv, rot, pe):
        # sanity checks
        assert len(pv) == source.num_vertices
        assert len(pv) == target.num_vertices
        assert len(pe) == source.num_edges
        assert len(pe) == target.num_edges
        assert set(pv.keys()) == set(range(source.num_vertices))
        assert set(pv.values()) == set(range(target.num_vertices)), \
            ("Got pv.values=%r but target.num_vertices=%r"
             " (and target=%r): number of elements do not match!"
             % (pv.values(), target.num_vertices, target))
        assert set(pe.keys()) == set(range(source.num_edges))
        assert set(pe.values()) == set(range(target.num_edges))

        self.source = source
        self.target = target
        self.pe = pe
        self.rot = rot
        self.pv = pv

    def __str__(self):
        return "(%s, %s, %s)" % (self.pv, self.rot, self.pe)

    def compare_orientations(self):
        """Return +1 or -1 depending on whether the orientations of
        the target Fatgraph pulls back to the orientation of the
        source Fatgraph via this `Isomorphism`.
        """
        image_edge_numbering = Permutation((self.source.edge_numbering[x],
                                            self.target.edge_numbering[self.pe[x]])
                                           for x in xrange(self.source.num_edges))
        return image_edge_numbering.sign()

    def is_identity(self):
        """
        Return `True` if this is the identity morphism.

        Only makes sense for automorphisms, so an `AssertionException`
        is raised if source and target differ.
        """
        assert self.source == self.target
        for v in xrange(len(self.pv)):
            if v != self.pv[v]:
                return False
            l = len(self.source.vertices[v])
            if self.rot[v] % l != 0:
                return False
        for e in xrange(len(self.pe)):
            if e != self.pe[e]:
                return False
        return True

    def is_orientation_reversing(self):
        """Return `True` if this `Isomorphism` reverses orientation on
        the source and target `Fatgraph` instances."""
        return (-1 == self.compare_orientations())

    def transform_boundary_cycle(self, bcy):
        """Return a new `BoundaryCycle` instance, obtained by
        transforming each corner according to a graph isomorphism.
        """
        triples = []
        for (v, i, j) in bcy:
            l = len(self.source.vertices[v])
            # create transformed triple 
            v_ = self.pv[v]
            i_ = (i + self.rot[v]) % l  # XXX: is it `-` or `+`?
            j_ = (j + self.rot[v]) % l
            # ensure the contract is honored, that `j` is the
            # index _following_ `i` in the cyclic order
            if i_ == 0 and j_ == l:
                i_, j_ = j_, i_
            triples.append((v_, i_, j_))
        return BoundaryCycle(triples)


class EqualIfIsomorphic(Caching):
    """Instances of this class will compare equal if there is an
    isomorphism mapping one to the other.
    """

    __slots__ = ('invariants')

    def __init__(self, invariants):
        self.invariants = invariants
        # set this instance's ``persistent id''
        Caching.__init__(self)

    @ocache_eq
    def __eq__(self, other):
        """Return `True` if `self` and `other` are isomorphic."""

        assert isinstance(other, type(self)), \
            "EqualIfIsomorphic.__eq__:" \
            " called with incompatible type arguments: `%s` and` %s`" % (type(self), type(other))

        # shortcuts
        if self is other:
            return True
        if self.invariants != other.invariants:
            return False

        # go the long way: try to find an explicit isomorphims
        # between `self` and `other`
        try:
            # if there is any morphism, then return `True`
            self.isomorphisms(other).next()
            return True
        except StopIteration:
            # list of morphisms is empty, objects are not equal.
            return False

    # both `__eq__` and `__ne__` are needed for testing equality of objects;
    # see `<http://www.voidspace.org.uk/python/articles/comparison.shtml>`
    def __ne__(self, other):
        """The opposite of `__eq__` (which see)."""
        return not self.__eq__(other)


class Fatgraph(EqualIfIsomorphic):
    """A fully-decorated ribbon graph.

    Several attributes of this object are computed at init time and
    store characteristics of the fatgraph:

      `.boundary_cycles`
        List of the boundary components of this `Fatgraph` object.
        Each boundary component is represented by the list of
        (colored) edges.

      `.genus`
        The genus of the Riemann surface this fatgraph lies on.
      
      `.num_boundary_cycles`
        Number of boundary cycles of this `Fatgraph`.

      `.num_edges`
        The number of edges of this `Fatgraph` object.

      `.num_vertices`
        Number of vertices of this `Fatgraph` object.

      `.vertices`
        List of vertices of this `Fatgraph`; each one is an instance
        of the `Vertex` class; see `Fatgraph.__init__` for examples.

        
    Examples::
      >>> Fatgraph([Vertex([2,1,0]), Vertex([2,1,0])]).num_boundary_cycles
      1
      >>> Fatgraph([Vertex([2,1,0]), Vertex([2,0,1])]).num_boundary_cycles
      3

    Two `Fatgraph`s compare equal if they are isomorphic::

      >>> Fatgraph([Vertex([1,0,0,1])]) == Fatgraph([Vertex([1,1,0,0])])
      True

      >>> Fatgraph([Vertex([2,0,0]), Vertex([2,1,1])]) \
            == Fatgraph([Vertex([2,2,0]), Vertex([1,1,0])])
      True

      >>> Fatgraph([Vertex([2,0,1]), Vertex([2,0,1])]) \
            == Fatgraph([Vertex([2,1,0]), Vertex([2,0,1])])
      False

      >>> Fatgraph([Vertex([2,0,1]), Vertex([2,0,1])]) \
            == Fatgraph([Vertex([2,0,0]), Vertex([2,1,1])])
      False

      >>> Fatgraph([Vertex([2,0,0]), Vertex([2,1,1])]) \
            == Fatgraph([Vertex([1,1,0,0])])
      False
    """

    __slots__ = (
        '__weakref__',
        'boundary_cycles',
        'edges',
        'edge_numbering',
        'genus',
        'invariants',
        'num_boundary_cycles',
        'num_edges',
        'num_vertices',
        'vertices',
    )

    def __init__(self, g_or_vs, **kwargs):
        """Construct a `Fatgraph` instance, taking list of vertices.

        Argument `g_or_vs` can be either:

          - a sequence of `Vertex` class instances::  

              >>> g1 = Fatgraph([Vertex([2,0,1]), Vertex([2,1,0])])

            Note that the list of vertices is assigned, *not copied*
            into the instance variable.

        or:

          - another `Fatgraph` instance, in which case this acts as a
            copy-constructor (and any remaining arguments are
            ignored), returning a new `Fatgraph` instance::

              >>> g2 = Fatgraph(g1)
              >>> g2 is g1
              False
              >>> g2 == g1
              True

            The returned `Fatgraph` instance *shares* all attributes
            with the instance given as argument::

              >>> g2.vertices is g1.vertices
              True
              >>> g2.edge_numbering is g1.edge_numbering
              True
              >>> g2.edges is g1.edges
              True
        """
        # dispatch based on type of arguments passed
        if isinstance(g_or_vs, Fatgraph):
            # copy-constructor used by class `NumberedFatgraph`
            self.boundary_cycles = g_or_vs.boundary_cycles
            self.edges = g_or_vs.edges
            self.edge_numbering = g_or_vs.edge_numbering
            self.genus = g_or_vs.genus
            self.num_boundary_cycles = g_or_vs.num_boundary_cycles
            self.num_edges = g_or_vs.num_edges
            self.num_vertices = g_or_vs.num_vertices
            self.vertices = g_or_vs.vertices

        else:  # initialize *new* instance

            #: list of vertices
            self.vertices = g_or_vs

            #: Number of edge colors
            self.num_edges = kwargs.get('num_edges',
                                        sum(len(v) for v in self.vertices) / 2)

            #: Number of vertices  XXX: why is this settable with kwarg???
            self.num_vertices = kwargs.get('num_vertices', len(g_or_vs))

            if 'edges' in kwargs:
                self.edges = kwargs.get('edges')
            else:
                #: Adjacency list of this graph.  For each edge, store a pair
                #: `(v1, v2)` where `v1` and `v2` are (indices of)
                #: endpoint vertices of an edge, and a corresponding pair
                #: `(i1, i2)` where `i1` and `i2` are indices of the given
                #: edge in vertices `v1` and `v2`.
                #: Each pair `(v, i)` represents a flag by the endpoint
                #: vertex and the index of the edge in the vertex.  (The
                #: vertex index alone is not enough for representing the
                #: edge arrow for loops.)
                endpoints = [[] for dummy in xrange(self.num_edges)]
                for current_vertex_index in xrange(self.num_vertices):
                    for (edge_index_in_vertex, edge) \
                            in enumerate(self.vertices[current_vertex_index]):
                        assert edge in range(self.num_edges), \
                            "Fatgraph.__init__:" \
                            " edge number %d not in range 0..%d" \
                            % (edge, self.num_edges)
                        endpoints[edge].append((current_vertex_index, edge_index_in_vertex))
                # now wrap endpoints into `Edge` objects
                self.edges = [Edge(*e) for e in endpoints]

            ## Orientation is given by an ordering of the edges,
            ## which directly translates into an orientation of the
            ## associated cell.  
            if 'orientation' in kwargs:
                self.edge_numbering = kwargs.get('orientation')
            else:
                self.edge_numbering = [x for x in xrange(self.num_edges)]

            self.boundary_cycles = self.compute_boundary_cycles()
            self.num_boundary_cycles = len(self.boundary_cycles)

            # by Euler, V-E+n=2-2*g
            self.genus = (self.num_edges - self.num_vertices
                          - self.num_boundary_cycles + 2) / 2

        # before computing invariants, check that internal data
        # structures are in a consistent state
        assert self.__ok()

        # used for isomorphism testing
        EqualIfIsomorphic.__init__(self, (
            self.num_vertices,
            self.num_edges,
            self.num_boundary_cycles,
        ))

    def __ok(self):
        """Perform coherency checks on internal state variables of
        `Fatgraph` instance and return `True` if they all pass.
        """
        assert self.num_edges > 0, \
            "Fatgraph `%s` has 0 edges." % (self)
        # check edge endpoints
        for (edgeno, edge) in enumerate(self.edges):
            assert len(edge.endpoints) == 2
            assert isinstance(edge.endpoints[0], tuple)
            assert isinstance(edge.endpoints[1], tuple)
            assert isinstance(edge.endpoints[0][0], int)
            assert isinstance(edge.endpoints[0][1], int)
            assert isinstance(edge.endpoints[1][0], int)
            assert isinstance(edge.endpoints[0][1], int)
            assert (0 <= edge.endpoints[0][0] < self.num_vertices)
            assert (0 <= edge.endpoints[1][0] < self.num_vertices)
            assert (0 <= edge.endpoints[0][1] < len(self.vertices[edge.endpoints[0][0]]))
            assert (0 <= edge.endpoints[1][1] < len(self.vertices[edge.endpoints[1][0]]))
            assert (edgeno in self.vertices[edge.endpoints[0][0]]), \
                "Invalid endpoint %s for edge %d of graph `%s`" \
                % (edge.endpoints[0], edgeno, self)
            assert (edgeno in self.vertices[edge.endpoints[1][0]]), \
                "Invalid endpoint %s for edge %d of graph `%s`" \
                % (edge.endpoints[1], edgeno, self)
        # check that each edge occurs exactly two times in vertices
        cnt = [0 for x in xrange(self.num_edges)]
        for v in self.vertices:
            for edgeno in v:
                cnt[edgeno] += 1
        for edgeno, cnt in enumerate(cnt):
            if edgeno < self.num_edges:
                assert cnt == 2, \
                    "Regular edge %d appears in %d vertices" \
                    % (edgeno, cnt)
            else:
                assert cnt == 1, \
                    "External edge %d appears in %d vertices" \
                    % (edgeno, cnt)

        assert self.edge_numbering is not None

        return True

    def __repr__(self):
        if not hasattr(self, 'vertices'):
            return "Fatgraph(<Initializing...>)"
        else:
            return "Fatgraph(%s)" % repr(self.vertices)

    def __str__(self):
        return repr(self)

    def automorphisms(self):
        """Enumerate automorphisms of this `Fatgraph` object.

        See `.isomorphisms()` for details of how a `Fatgraph`
        isomorphism is represented.
        """
        return self.isomorphisms(self)

    # @cython.locals(vertices=list, corners=list,
    #               v=cython.int, i=cython.int, corner=tuple)
    @staticmethod
    def _first_unused_corner(vertices, corners):
        # fast-forward to the first unused corner
        for v in xrange(len(vertices)):
            for i in xrange(len(vertices[v])):
                corner = corners[v][i]
                if corner is not None:
                    return (v, i, corner)
        return (v, i, None)

    def compute_boundary_cycles(self):
        """Return a list of boundary cycles of this `Fatgraph` object.

        Boundary cycles are represented as a cyclic list of 'corners':
        a corner is a triple `(v, i, j)` consisting of a vertex and
        two consecutive indices (in the cyclic order, so, either `j ==
        i+1` or `i` and `j` are the starting and ending indices)::
        
          >>> Fatgraph([Vertex([2,1,0]),Vertex([2,0,1])]).compute_boundary_cycles()
          [BoundaryCycle([(1, 2, 0), (0, 0, 1)]),
           BoundaryCycle([(0, 1, 2), (1, 1, 2)]),
           BoundaryCycle([(0, 2, 0), (1, 0, 1)])]

        This verbose representation allows one to distinguish the
        boundary cycles made from the same set of edges::

          >>> Fatgraph([Vertex([0,1,2,0,1,2])]).compute_boundary_cycles()
          [BoundaryCycle([(0, 2, 3), (0, 4, 5), (0, 0, 1)]),
           BoundaryCycle([(0, 1, 2), (0, 3, 4), (0, 5, 0)])]
        """

        # Build the collection of "corners" of `graph`,
        # structured just like the set of vertices.
        # By construction, `corners[v][i]` has the the
        # form `(v,i,j)` where `j` is the index following
        # `i` in the cyclic order.
        corners = [[(v, i, (i + 1) % len(self.vertices[v]))
                    for i in xrange(len(self.vertices[v]))]
                   for v in xrange(self.num_vertices)]

        result = []
        while True:
            # fast-forward to the first unused corner
            (v, i, corner) = Fatgraph._first_unused_corner(self.vertices, corners)
            if corner is None:
                break

            # build a list of corners comprising the same boundary
            # cycle: start with one corner, follow the edge starting
            # at the second delimiter of the corner to its other
            # endpoint, and repeat until we come back to the starting
            # point.  
            start = (v, i)
            triples = []
            while (v, i) != start or len(triples) == 0:
                assert corners[v][i] is not None
                corner = corners[v][i]
                corners[v][i] = None
                triples.append(corner)
                assert v == corner[0]
                assert i == corner[1]
                j = corner[2]
                edgeno = self.vertices[v][j]
                (v, i) = self.edges[edgeno].other_end(v, j)
            result.append(BoundaryCycle(triples))

        return result

    def bridge(self, edge1, side1, edge2, side2):
        """Return a new `Fatgraph`, formed by inserting trivalent
        vertices in the middle of edges `edge1` and `edge2` and
        connecting them with a new edge.

          >>> g = Fatgraph([Vertex([0,1,2]), Vertex([0,2,1])])
          >>> g1 = g.bridge(0, 0, 1, 1)
          >>> g1 is g
          False
          >>> g1 == g
          False
          
        Arguments `side1` and `side2` control which side the new edge
        is attached to (valid values are 0 or 1), i.e., which of the
        two inequivalent cyclic orders the new trivalent vertices will
        be given::
        
          >>> g = Fatgraph([Vertex([0,1,2]), Vertex([0,2,1])])
          >>> g1 = g.bridge(0, 0, 1, 0)
          >>> g2 = g.bridge(0, 1, 1, 0)
          >>> g1 == g2
          False

        In more detail: let 0,1,2 be the indices of the edges attached
        to the new vertex in the middle of `edge1`, where 0,1 denote
        the two halves of `edge1`.  If `side1` is `0`, then the new
        trivalent vertex will have the cyclic order [0,1,2]; if
        `side1` is `1`, then 0,1 are swapped and the new trivalent
        vertex gets the cyclic order [1,0,2]::

          >>> g1 == Fatgraph([Vertex([0,1,2]), Vertex([4,2,5]), Vertex([0,4,3]), Vertex([1,5,3])])
          True
          >>> g2 == Fatgraph([Vertex([0,1,2]), Vertex([4,2,5]), Vertex([4,0,3]), Vertex([1,5,3])])
          True

        It is worth noting that this procedure involves 5 edges in
        total, 3 of which need new indices.

          >>> g.num_edges
          3
          >>> g1.num_edges
          6
          
        This function is obviously symmetric: the pairs `edge1, side1`
        and `edge2, side2` can be swapped and the result stays the
        same::

          >>> g3 = g.bridge(1, 0, 0, 0)
          >>> g1 == g3
          True

          >>> g4 = g.bridge(1, 0, 0, 1)
          >>> g2 == g4
          True

        Examples::
        
        1) Bridging different sides of the same edge may yield
        different results::
        
          >>> g.bridge(0, 0, 1, 0)  == Fatgraph([Vertex([0,1,2]), Vertex([4,2,5]), Vertex([0,4,3]), Vertex([1,5,3])])
          True
          
          >>> g.bridge(0, 1, 1, 0) == Fatgraph([Vertex([0,1,2]), Vertex([4,2,5]), Vertex([4,0,3]), Vertex([1,5,3])])
          True

        2) One can connect an edge to itself on different sides::
        
          >>> g.bridge(0, 0, 0, 1) == Fatgraph([Vertex([0,1,2]), Vertex([5,2,1]), Vertex([0,4,3]), Vertex([5,4,3])])
          True

        3) And also with both ends on the same side::
        
          >>> g.bridge(0, 0, 0, 0) == Fatgraph([Vertex([0,1,2]), Vertex([5,2,1]), Vertex([0,4,3]), Vertex([4,5,3])])
          True
          
        """
        assert side1 in [0, 1], \
            "Fatgraph.bridge: Invalid value for `side1`: '%s' - should be 0 or 1" % side1
        assert side2 in [0, 1], \
            "Fatgraph.bridge: Invalid value for `side2`: '%s' - should be 0 or 1" % side2

        opposite_side1 = 0 if side1 == 1 else 1
        opposite_side2 = 0 if side2 == 1 else 1

        ## assign edge indices
        connecting_edge = self.num_edges
        ## break `edge1` in two halves: if `v1a` and `v1b` are the
        ## endpoints of `edge1`, then the "one_half" edge extends from
        ## the `v1a` endpoint of `edge1` to the new vertex
        ## `midpoint1`; the "other_half" edge extends from the
        ## `midpoint1` new vertex to `v1b`.
        one_half1 = edge1
        other_half1 = self.num_edges + 1
        ## break `edge2` in two halves; if `edge2` is the same edge as
        ## `edge1`, then we are breaking the second half of `edge1` in
        ## two parts.  Otherwise, proceed as above.  In any case, the
        ## "other half" of `edge2` needs a new edge index.
        if edge2 == edge1:
            one_half2 = other_half1
        else:
            one_half2 = edge2
        other_half2 = self.num_edges + 2

        ## assign new vertex indices
        midpoint1_index = self.num_vertices
        midpoint2_index = self.num_vertices + 1

        if side1 == 1:
            midpoint1 = Vertex([other_half1, one_half1, connecting_edge])
        else:  # side2 == 0
            midpoint1 = Vertex([one_half1, other_half1, connecting_edge])

        if side2 == 1:
            midpoint2 = Vertex([other_half2, one_half2, connecting_edge])
        else:  # side2 == 0
            midpoint2 = Vertex([one_half2, other_half2, connecting_edge])

        ## two new vertices are added: the mid-points of the connected edges.
        new_vertices = self.vertices + [midpoint1, midpoint2]

        ## three new edges are added (constructed below)
        new_edges = self.edges + [
            None,  # new edge: connecting_edge
            None,  # new edge: other_half1
            None,  # new edge: other_half2
        ]

        ## the connecting edge has endpoints in the mid-points of
        ## `edge1` and `edge2`, and is *always* in third position.
        new_edges[connecting_edge] = Edge((midpoint1_index, 2), (midpoint2_index, 2))

        ((v1a, pos1a), (v1b, pos1b)) = self.edges[edge1].endpoints
        new_edges[one_half1] = Edge((v1a, pos1a), (midpoint1_index, side1))
        if edge1 != edge2:
            # replace `edge1` with new `other_half1` in the second endpoint
            new_vertices[v1b] = Vertex(new_vertices[v1b][:pos1b]
                                       + [other_half1]
                                       + new_vertices[v1b][pos1b + 1:])
            new_edges[other_half1] = Edge((midpoint1_index, opposite_side1), (v1b, pos1b))
        else:
            # same edge, "other half" ends at the second endpoint
            new_edges[other_half1] = Edge((midpoint1_index, opposite_side1), (midpoint2_index, side2))

        # replace `edge2` with new `other_half2` in the second
        # endpoint; again we need to distinguish the special case when
        # `edge1` and `edge2` are the same edge.
        ((v2a, pos2a), (v2b, pos2b)) = self.edges[edge2].endpoints
        if edge1 != edge2:
            new_edges[one_half2] = Edge((v2a, pos2a), (midpoint2_index, side2))
        else:
            # `edge1 == edge2`, so `one_half2 == other_half1`
            pass  # new_edges[one_half2] = new_edges[other_half1]
        # "other half" of second edge *always* ends at the previous
        # edge endpoint, so replace `edge2` in `v2b`.
        new_vertices[v2b] = Vertex(new_vertices[v2b][:pos2b]
                                   + [other_half2]
                                   + new_vertices[v2b][pos2b + 1:])
        new_edges[other_half2] = Edge((midpoint2_index, opposite_side2), (v2b, pos2b))

        ## inherit orientation, and add the three new edges in the order they were created
        new_edge_numbering = self.edge_numbering + \
                             [connecting_edge, other_half1, other_half2]

        # build new graph 
        return Fatgraph(new_vertices,
                        edges=new_edges,
                        num_edges=self.num_edges + 3,
                        orientation=new_edge_numbering,
                        )

    def bridge2(self, edge1, side1, other, edge2, side2):
        """Return a new `Fatgraph`, formed by connecting the midpoints
        of `edge1` on `self` and `edge2` on `other`.
        
          >>> g1 = Fatgraph([Vertex([0,1,2]), Vertex([0,2,1])])
          >>> g2 = Fatgraph([Vertex([0,1,2,0,1,2])])
          >>> g = Fatgraph.bridge2(g1, 0, 0, g2, 1, 1)
          >>> g is g1
          False
          >>> g is g2
          False
          >>> g == Fatgraph([Vertex([0,1,2]), Vertex([6,2,1]), Vertex([3,4,5,3,7,5]), Vertex([0,6,8]), Vertex([4,7,8])])
          True
          
        New trivalent vertices are inserted in the middle of the
        connected edges.  Arguments `side1` and `side2` control which
        side the new edge is attached to (valid values are 0 or 1),
        i.e., which of the two inequivalent cyclic orders the new
        trivalent vertices will be given (see `Fatgraph.bridge()`).
        
        It is worth noting that this procedure adds 3 edges and 2
        vertices to the edge total of `self` and `other`::
        
          >>> g.num_edges == g1.num_edges + g2.num_edges + 3
          True
          >>> g.num_vertices == g1.num_vertices + g2.num_vertices + 2
          True
          
        This function is obviously symmetric: the triplets `self, edge1, side1`
        and `other, edge2, side2` can be swapped and the result stays the
        same (up to isomorphisms)::

          >>> g_ = Fatgraph.bridge2(g2, 1, 1, g1, 0, 0)
          >>> g == g_
          True

        *Caveat:* If `self == other` then the resulting graph is made
        up of *two copies* of `self` with a new edge connecting
        `edge1` on one copy and `edge2` on the other.
         
        """
        assert side1 in [0, 1], \
            "Fatgraph.bridge2: Invalid value for `side1`: '%s' -- should be 0 or 1" % side1
        assert side2 in [0, 1], \
            "Fatgraph.bridge2: Invalid value for `side2`: '%s' -- should be 0 or 1" % side2

        ## First, build a (non-connected) graph from the disjoint
        ## union of `self` and `other`.

        # Edges of `other` are renumbered depending on whether
        # they are internal of external edges:
        #   - internal edges in `other` have numbers ranging from 0 to
        #     `other.num_edges`: they get new numbers starting from
        #     `self.num_edges` and counting upwards
        renumber_other_edges = dict((x, x + self.num_edges)
                                    for x in xrange(other.num_edges))
        # Orientation needs the same numbering:
        new_edge_numbering = self.edge_numbering \
                             + ltranslate(renumber_other_edges, other.edge_numbering)
        # Similarly, vertices of `self` retain indices `[0..v]`, while
        # vertices of `other` follow.
        new_vertices = self.vertices \
                       + [Vertex(ltranslate(renumber_other_edges, ov))
                          for ov in other.vertices]
        renumber_other_vertices = dict((x, x + self.num_vertices)
                                       for x in xrange(other.num_vertices))
        # build new edges: vertex indices need to be shifted for
        # endpoints, but attachment indices are the same; three new
        # edges are added at the tail of the list
        new_edges = self.edges \
                    + [Edge((renumber_other_vertices[x.endpoints[0][0]], x.endpoints[0][1]),
                            (renumber_other_vertices[x.endpoints[1][0]], x.endpoints[1][1]))
                       for x in other.edges] \
                    + [None,  # connecting_edge
                       None,  # other_half1
                       None]  # other_half2

        edge2 = renumber_other_edges[edge2]  # index changed in `new_edges`

        opposite_side1 = 0 if side1 else 1
        opposite_side2 = 0 if side2 else 1

        ## assign edge indices
        connecting_edge = self.num_edges + other.num_edges
        ## break `edge1` in two halves: if `v1a` and `v1b` are the
        ## endpoints of `edge1`, then the "one_half" edge extends from
        ## the `v1a` endpoint of `edge1` to the new vertex
        ## `midpoint1`; the "other_half" edge extends from the
        ## `midpoint1` new vertex to `v1b`.
        one_half1 = edge1
        other_half1 = connecting_edge + 1
        ## break `edge2` in two halves; same as above.
        one_half2 = edge2
        other_half2 = connecting_edge + 2

        ## assign new vertex indices
        midpoint1_index = len(new_vertices)
        midpoint2_index = midpoint1_index + 1

        if side1:
            midpoint1 = Vertex([other_half1, one_half1, connecting_edge])
        else:
            midpoint1 = Vertex([one_half1, other_half1, connecting_edge])

        if side2:
            midpoint2 = Vertex([other_half2, one_half2, connecting_edge])
        else:
            midpoint2 = Vertex([one_half2, other_half2, connecting_edge])

        ## two new vertices are added: the mid-points of the connected edges.
        new_vertices += [midpoint1, midpoint2]
        ## the connecting edge has endpoints in the mid-points of
        ## `edge1` and `edge2`, and is *always* in third position.
        new_edges[connecting_edge] = Edge((midpoint1_index, 2), (midpoint2_index, 2))

        ((v1a, pos1a), (v1b, pos1b)) = new_edges[edge1].endpoints
        new_edges[one_half1] = Edge((v1a, pos1a), (midpoint1_index, side1))
        # replace `edge1` with new `other_half1` in the second endpoint
        new_vertices[v1b] = Vertex(new_vertices[v1b][:pos1b]
                                   + [other_half1]
                                   + new_vertices[v1b][pos1b + 1:])
        new_edges[other_half1] = Edge((midpoint1_index, opposite_side1), (v1b, pos1b))

        # replace `edge2` with new `other_half2` in the second
        # endpoint; again we need to distinguish the special case when
        # `edge1` and `edge2` are the same edge.
        ((v2a, pos2a), (v2b, pos2b)) = new_edges[edge2].endpoints
        new_edges[one_half2] = Edge((v2a, pos2a), (midpoint2_index, side2))
        # "other half" of second edge *always* ends at the previous
        # edge endpoint, so replace `edge2` in `v2b`.
        new_vertices[v2b] = Vertex(new_vertices[v2b][:pos2b]
                                   + [other_half2]
                                   + new_vertices[v2b][pos2b + 1:])
        new_edges[other_half2] = Edge((midpoint2_index, opposite_side2), (v2b, pos2b))

        ## inherit orientation, and add the three new edges in the order they were created
        new_edge_numbering += [connecting_edge, other_half1, other_half2]

        # build new graph 
        return Fatgraph(new_vertices,
                        edges=new_edges,
                        num_edges=self.num_edges + other.num_edges + 3,
                        orientation=new_edge_numbering,
                        )

    @ocache_contract
    def contract(self, edge):
        """Return new `Fatgraph` obtained by contracting the specified edge.

        Examples::

          >>> Fatgraph([Vertex([2,2,0]), Vertex([0,1,1])]).contract(0)
          Fatgraph([Vertex([1, 1, 0, 0])])
          >>> Fatgraph([Vertex([2,1,0]), Vertex([2,0,1])]).contract(1)
          Fatgraph([Vertex([0, 1, 1, 0])])

        The M_{1,1} trivalent graph yield the same result no matter
        what edge is contracted::

          >>> Fatgraph([Vertex([2,1,0]), Vertex([2,1,0])]).contract(0)
          Fatgraph([Vertex([1, 0, 1, 0])])
          >>> Fatgraph([Vertex([2,1,0]), Vertex([2,1,0])]).contract(1)
          Fatgraph([Vertex([0, 1, 0, 1])])
          >>> Fatgraph([Vertex([2,1,0]), Vertex([2,1,0])]).contract(2)
          Fatgraph([Vertex([1, 0, 1, 0])])
        """
        assert not self.is_loop(edge), \
            "Fatgraph.contract: cannot contract a loop."
        assert (edge >= 0) and (edge < self.num_edges), \
            "Fatgraph.contract: invalid edge number (%d):" \
            " must be in range 0..%d" \
            % (edge, self.num_edges)

        ## Plug the higher-numbered vertex into the lower-numbered one.

        # store endpoints of the edge-to-be-contracted
        ((v1, pos1), (v2, pos2)) = self.edges[edge].endpoints
        assert v1 < v2

        # Build new list of vertices, removing the contracted edge and
        # shifting all indices above:
        #   - edges numbered 0..edge-1 are unchanged;
        #   - edges numbered `edge+1`.. are renumbered, 
        #     shifting the number down one position;
        #   - edge `edge` is kept intact, will be removed by mating
        #     operation (see below).
        renumber_edges = dict((i, i)
                              for i in xrange(0, edge + 1))
        renumber_edges.update(dict((i, i - 1)
                                   for i in xrange(edge + 1, self.num_edges)))
        # See `itranslate` in utils.py for how this prescription is
        # encoded in the `renumber_edges` mapping.
        new_vertices = [Vertex([renumber_edges[e] for e in V])
                        for V in self.vertices]

        # Mate endpoints of contracted edge:
        # 1. Rotate endpoints `v1`, `v2` so that the given edge would
        #    appear *last* in `v1` and *first* in `v2` (*Note:* since
        #    `v1`, `v2` are *cyclic*, this means that we do the same
        #    operation on `v1` and `v2` alike).
        # 2. Join vertices by concatenating the list of incident
        #    edges;
        # 3. Set new `i1` vertex in place of old first endpoint:
        new_vertices[v1] = Vertex(
            Fatgraph._rotate_and_trim(new_vertices[v1], pos1)
            +
            Fatgraph._rotate_and_trim(new_vertices[v2], pos2)
        )
        # 4. Remove second endpoint from list of new vertices:
        del new_vertices[v2]

        ## Orientation of the contracted graph.
        cut = self.edge_numbering[edge]
        # edges with index below the contracted one are untouched
        renumber_edge_numbering = dict((x, x) for x in xrange(cut))
        # edges with index above the contracted one are shifted down
        # one position
        renumber_edge_numbering.update(dict((x + 1, x)
                                            for x in xrange(cut, self.num_edges)))
        new_edge_numbering = [renumber_edge_numbering[self.edge_numbering[x]]
                              for x in xrange(self.num_edges)
                              if x != edge]

        # build new graph
        return Fatgraph(
            new_vertices,
            orientation=new_edge_numbering,
        )

    # @cython.locals(L=list, p=cython.int)
    # @cython.cfunc(list)
    @staticmethod
    def _rotate_and_trim(V, p):
        """
        Return a copy of the given list `V`, with items rotated `p`
        positions leftwards; additionally, the item originally at
        position `p` is removed.

        Examples::

          >>> Fatgraph._rotate_and_trim([1,2,3,4,5], 2)
          [4, 5, 1, 2]

        """
        return V[p + 1:] + V[:p]

    def contract_boundary_cycle(self, bcy, vi1, vi2):
        """Return a new `BoundaryCycle` instance, image of `bcy` under
        the topological map that contracts the edge with endpoints
        `(v1,i1)` and `(v2,i2)` that are passed as first and second
        argument.

        XXX: return `bcy` if neither `v1` nor `v2` are contained in it.
        """
        (v1, pos1) = vi1
        (v2, pos2) = vi2
        l1 = len(self.vertices[v1])
        l2 = len(self.vertices[v2])
        new_bcy = []
        for corner in bcy:
            if corner[0] == v1:
                if pos1 == corner[1]:
                    # skip this corner, keep only one of the
                    # corners limited by the contracted edge
                    continue
                else:
                    i1 = (corner[1] - pos1 - 1) % l1
                    i2 = (corner[2] - pos1 - 1) % l1
                    assert (i1 + 1 - i2) % l1 == 0  # i1,i2 denote successive indices
                    assert i1 != l1 - 1  # would collide with contracted corners from `v2`
                    new_bcy.append((v1, i1, i2))
            elif corner[0] == v2:
                if pos2 == corner[1]:
                    # skip this corner, keep only one of the
                    # corners limited by the contracted edge
                    continue
                if pos2 == corner[2]:
                    new_bcy.append((v1, l1 + l2 - 3, 0))
                else:
                    i1 = l1 - 1 + ((corner[1] - pos2 - 1) % l2)
                    i2 = l1 - 1 + ((corner[2] - pos2 - 1) % l2)
                    assert (i1 + 1 - i2) % l1 == 0  # i1,i2 denote successive indices
                    new_bcy.append((v1, i1, i2))
            elif corner[0] > v2:
                # shift vertices after `v2` one position down
                new_bcy.append((corner[0] - 1, corner[1], corner[2]))
            else:
                # pass corner unchanged
                new_bcy.append(corner)
        if __debug__:
            cnt = {}
            for corner in new_bcy:
                try:
                    cnt[corner] += 1
                except KeyError:
                    cnt[corner] = 1
            for (corner, count) in cnt.iteritems():
                assert count == 1, \
                    "BoundaryCycle.contract():" \
                    " Corner %s appears %d times in contracted boundary cycle %s" \
                    % (corner, count, new_bcy)
        return BoundaryCycle(new_bcy)

    @ocache0
    def edge_orbits(self):
        """Compute orbits of the edges under the action of graph
        automorphism group, and a representative for each orbit.
        
        Returns a dictionary, whose keys are the representatives, and
        whose values are the orbits.  Orbits are represented as Python
        `set` objects.

        Examples::

          >>> Fatgraph([Vertex([0,1,2]), Vertex([0,1,2])]).edge_orbits()
          {0: set([0, 1, 2])}

          >>> Fatgraph([Vertex([1, 0, 2]), Vertex([2, 1, 0])]).edge_orbits()
          {0: set([0, 1, 2])}
          
        """
        orbits = dict((x, set([x])) for x in xrange(self.num_edges))
        for a in self.automorphisms():
            edge_permutation = a.pe
            for x in xrange(self.num_edges):
                if x not in orbits:
                    continue
                y = edge_permutation[x]
                if y not in orbits:
                    continue
                # `x` and `y` are in the same orbit, only keep the one
                # with lower abs. value, and remove the other.
                if y > x:
                    orbits[x].update(orbits[y])
                    del orbits[y]
        # check that all elements lie in some orbit
        assert sum(len(set(o)) for o in orbits.itervalues()) == self.num_edges, \
            "Fatgraph.edge_orbits():" \
            " Computed orbits `%s` do not exhaust edge set `%s`" \
            " [%s.edge_orbits() -> %s]" % (orbits, range(self.num_edges), self, orbits)
        return orbits

    @ocache0
    def edge_pair_orbits(self):
        """Compute orbits of pairs `(edge1, edge2)` under the action
        of graph automorphism group, and a representative for each
        orbit.
        
        Returns a dictionary, whose keys are the representatives, and
        whose values are the orbits.  Orbits are represented as Python
        `set` objects.

        Examples::

          >>> Fatgraph([Vertex([0,1,2]), Vertex([0,1,2])]).edge_pair_orbits()
          {(0, 1): set([(0, 1), (1, 2), (2, 0)]),
           (0, 0): set([(0, 0), (1, 1), (2, 2)]),
           (0, 2): set([(1, 0), (0, 2), (2, 1)])}
          
        """
        edge_pairs = [(x, y)
                      for x in xrange(self.num_edges)
                      for y in xrange(self.num_edges)]
        orbits = dict((p, set([p])) for p in edge_pairs)
        for a in self.automorphisms():
            edge_permutation = a.pe
            for p in edge_pairs:
                if p not in orbits:
                    continue
                q = (edge_permutation[p[0]], edge_permutation[p[1]])
                if q not in orbits:
                    continue
                # `p` and `q` are in the same orbit, only keep the one
                # with lower abs. value, and remove the other.
                if p < q:
                    orbits[p].update(orbits[q])
                    del orbits[q]
        # check that all elements lie in some orbit
        assert sum(len(set(o)) for o in orbits.itervalues()) == len(edge_pairs), \
            "Fatgraph.edge_pair_orbits():" \
            " Computed orbits `%s` do not exhaust edge pairs set `%s`" \
            " [%s.edge_pair_orbits() -> %s]" % (orbits, edge_pairs, self, orbits)
        return orbits

    def endpoints(self, edgeno):
        """Return the endpoints of `edge`, as a pair of `(v, pos)`
        where `v` is the endpoint vertex index, and `pos` is the
        attachment index of `edge` into the `Vertex` object
        `self.vertices[v]`.

        The pair `((v1, pos1), (v2, pos2))` is ordered such that `v1 < v2`.
        """
        return self.edges[edgeno].endpoints

    def hangcircle(self, edge, side):
        """Return a new `Fatgraph`, formed by attaching a circle with
        a new edge to a new trivalent vertex in the middle of `edge`.

          >>> g = Fatgraph([Vertex([0,1,2]), Vertex([0,2,1])])
          >>> g1 = g.hangcircle(0, 0)
          >>> g1 is g
          False
          
        Argument `side` controls which side of `edge` the circle is
        hung to (valid values are 0 or 1), i.e., which of the two
        inequivalent cyclic orders the new trivalent vertices will be
        given::
        
          >>> g = Fatgraph([Vertex([0,1,2]), Vertex([0,2,1])])
          >>> g1 = g.hangcircle(0, 0)
          >>> g1 == Fatgraph([Vertex([0,1,2]), Vertex([3,2,1]), Vertex([0,3,4]), Vertex([5,5,4])])
          True
          >>> g2 = g.hangcircle(0, 1)
          >>> g2 == Fatgraph([Vertex([0,1,2]), Vertex([3,2,1]), Vertex([3,0,4]), Vertex([5,5,4])])
          True

        It is worth noting that the new graph will have 3 edges more
        than the original one::

          >>> g1.num_edges == g.num_edges + 3
          True
          
        """
        assert side in [0, 1], \
            "Fatgraph.hangcircle: Invalid value for `side`: '%s' - should be 0 or 1" % side

        opposite_side = 0 if side else 1

        ## assign edge indices

        ## break `edge` in two halves: if `v1` and `v2` are the
        ## endpoints of `edge`, then the "one_half" edge extends from
        ## the `v1` endpoint of `edge1` to the new vertex
        ## `midpoint`; the "other_half" edge extends from the
        ## `midpoint` new vertex to `v2`.
        one_half = edge
        other_half = self.num_edges
        connecting_edge = self.num_edges + 1
        circling_edge = self.num_edges + 2

        ## assign new indices to new vertices
        midpoint_index = self.num_vertices
        T_index = self.num_vertices + 1

        ## two new vertices are added: the mid-point of `edge`, and
        ## the vertex `T` lying on the circle.
        if side == 1:
            midpoint = Vertex([other_half, one_half, connecting_edge])
        else:  # side == 0
            midpoint = Vertex([one_half, other_half, connecting_edge])
        T = Vertex([circling_edge, circling_edge, connecting_edge])
        new_vertices = self.vertices + [midpoint, T]

        ## new edges:
        ## - inherit edges from parent, and add place for three new edges
        new_edges = self.edges + [
            None,  # new edge: other_half
            None,  # new edge: connecting_edge
            None,  # new edge: circling_edge
        ]

        ## - break `edge` into two edges `one_half` and `other_half`:
        ((v1, pos1), (v2, pos2)) = self.edges[edge].endpoints
        new_edges[one_half] = Edge((v1, pos1), (midpoint_index, side))
        new_edges[other_half] = Edge((midpoint_index, opposite_side), (v2, pos2))
        new_vertices[v2] = Vertex(new_vertices[v2][:pos2]
                                  + [other_half]
                                  + new_vertices[v2][pos2 + 1:])

        ## - the connecting edge has endpoints in the mid-point of
        ## `edge` and in `T`, and is *always* in third position:
        new_edges[connecting_edge] = Edge((midpoint_index, 2), (T_index, 2))

        ## - the circling edge is a loop with vertex `T`
        new_edges[circling_edge] = Edge((T_index, 0), (T_index, 1))

        ## Inherit edge numbering from parent and extend as identity
        ## on the newly-added edges.
        new_edge_numbering = self.edge_numbering + \
                             [other_half, connecting_edge, circling_edge]

        # finally, build new graph 
        return Fatgraph(new_vertices,
                        edges=new_edges,
                        num_edges=self.num_edges + 3,
                        orientation=new_edge_numbering,
                        )

    def is_loop(self, edge):
        """Return `True` if `edge` is a loop (i.e., the two endpoint coincide).
        """
        return self.edges[edge].is_loop()

    def is_oriented(self):
        """Return `True` if `Fatgraph` is orientable.

        A `Fatgraph` is orientable iff it has no orientation-reversing
        automorphism.

        Enumerate all automorphisms, end exit with `False` result as
        soon as one orientation-reversing one is found.

        Examples::

          >>> Fatgraph([Vertex([1,0,1,0])]).is_oriented()
          False

          >>> Fatgraph([Vertex([2, 0, 1]), Vertex([2, 0, 1])]).is_oriented()
          True
          
          >>> Fatgraph([Vertex([2, 1, 0]), Vertex([2, 0, 1])]).is_oriented()
          False
          
          >>> Fatgraph([Vertex([2, 1, 1]), Vertex([2, 0, 0])]).is_oriented()
          False

          >>> Fatgraph([Vertex([3, 2, 2, 0, 1]), Vertex([3, 1, 0])], \
                    numbering=[(0, CyclicTuple((2,))),  \
                               (1, CyclicTuple((0, 1))),  \
                               (2, CyclicTuple((3, 1))),  \
                               (3, CyclicTuple((0, 3, 2))) ]) \
                               .is_oriented()
          True
          >>> Fatgraph([Vertex([2, 3, 1]), Vertex([2, 1, 3, 0, 0])], \
                       numbering=[(0, CyclicTuple((0,))),  \
                                  (2, CyclicTuple((1, 3))),  \
                                  (3, CyclicTuple((3, 0, 2))),  \
                                  (1, CyclicTuple((2, 1))) ]) \
                               .is_oriented()
          True
        """
        ## Try to find an orientation-reversing automorphism the hard way
        for a in self.automorphisms():
            if a.is_orientation_reversing():
                return False
        # no orientation reversing automorphism found
        return True

    @ocache_isomorphisms
    def isomorphisms(G1, G2):
        """Iterate over `Fatgraph` isomorphisms from `G1` to `G2`.

        An isomorphism is represented by a tuple `(pv, rot, pe)` where:

          - `pv` is a permutation of ther vertices: the `i`-th vertex
            of `G1` is sent to the `pv[i]`-th vertex of `G2`, rotated
            by `rot[i]` places leftwards;

          - `pe` is a permutation of the edges: edge `i` in `G1` is
            mapped to edge `pe[i]` in `G2`.

        This method can iterate over the automorphism group of a
        graph::

          >>> G1 = Fatgraph([Vertex([2, 1, 1]), Vertex([2, 0, 0])])
          >>> for f in G1.isomorphisms(G1): print f
          ({0: 0, 1: 1}, [0, 0], {0: 0, 1: 1, 2: 2})
          ({0: 1, 1: 0}, [0, 0], {0: 1, 1: 0, 2: 2})

        Or it can find the isomorphisms between two given graphs::

          >>> G2 = Fatgraph([Vertex([2, 2, 0]), Vertex([1, 1, 0])])
          >>> for f in G1.isomorphisms(G2): print f
          ({0: 0, 1: 1}, [2, 2], {0: 1, 1: 2, 2: 0})
          ({0: 1, 1: 0}, [2, 2], {0: 2, 1: 1, 2: 0})

        If there are no isomorphisms connecting the two graphs, then no
        item is returned by the iterator::

          >>> g3 = Fatgraph([Vertex([2, 1, 0]), Vertex([2, 0, 1])])
          >>> list(G1.isomorphisms(g3))
          []

        A concrete example taken from `M_{1,4}`:latex: ::

          >>> g1 = Fatgraph([Vertex([1, 0, 2]), Vertex([2, 1, 5]), Vertex([0, 4, 3]), Vertex([8, 5, 6]), Vertex([3, 6, 7, 7]), Vertex([8, 4, 9, 9])])
          >>> g2 = Fatgraph([Vertex([1, 0, 5, 6]), Vertex([1, 0, 2]), Vertex([5, 2, 3]), Vertex([8, 4, 3]), Vertex([7, 7, 6]), Vertex([4, 8, 9, 9])])
          >>> len(list(Fatgraph.isomorphisms(g1, g2)))
          0
        """
        # if graphs differ in vertex valences, no isomorphisms
        vs1 = G1._valence_spectrum()
        vs2 = G2._valence_spectrum()
        if not set(vs1.keys()) == set(vs2.keys()):
            return  # StopIteration
        # if graphs have unequal vertex distribution by valence, no isomorphisms
        for val in G1.vertex_valences():
            if len(vs1[val]) != len(vs2[val]):
                return  # StopIteration

        (valence, indexes) = G2._starting_vertices()
        v1_index = vs1[valence][0]
        v1 = G1.vertices[v1_index]
        for v2_index in Fatgraph._compatible_vertices(v1, G2, indexes):
            for rot in xrange(valence):
                try:
                    # pass 0: init new (pv, rots, pe) triple
                    pv = Permutation()
                    rots = [None for x in xrange(G1.num_vertices)]
                    pe = Permutation()

                    # pass 1: map `v1` to `v2` and build map
                    # of neighboring vertices for next pass
                    if not pe.extend(v1, G2.vertices[v2_index][rot:rot + valence]):
                        continue  # to next `rot`
                    pv[v1_index] = v2_index
                    rots[v1_index] = rot
                    if __debug__:
                        for x in v1:
                            assert x in pe, "Edge `%d` of vertex `%s` (in graph `%s`) not mapped to any edge of graph `%s` (at line 1740, `pe=%s`)" % (
                            x, v1, G1, G2, pe)

                    # pass 2: extend map to neighboring vertices
                    nexts = Fatgraph._neighbors(pv, pe, G1, v1_index, G2, v2_index)
                    while len(pv) < G1.num_vertices:
                        neighborhood = []
                        for (i1, i2, r) in nexts:
                            (pv, rots, pe) = Fatgraph._extend_map(pv, rots, pe, G1, i1, r, G2, i2)
                            if __debug__:
                                for x in G1.vertices[i1]:
                                    assert x in pe, "Edge `%d` of vertex `%s` (in graph `%s`) not mapped to any edge of graph `%s` (at line 1751, `pe=%s`)" % (
                                    x, G1.vertices[i1], G1, G2, pe)
                            neighborhood += Fatgraph._neighbors(pv, pe, G1, i1, G2, i2)
                        nexts = neighborhood

                # extension failed in the above block, continue with next candidate
                except Fatgraph._CannotExtendMap:
                    continue  # to next `rot`

                # finally
                yield Isomorphism(G1, G2, pv, rots, pe)

    ## auxiliary functions for `Fatgraph.isomorphism`

    @ocache0
    def _valence_spectrum(self):
        """Return a dictionary mapping valences into vertex indices.

        Examples::

           >>> Fatgraph([Vertex([1,1,0,0])])._valence_spectrum()
           defaultdict(<type 'list'>, {4: [0]})

           >>> Fatgraph([Vertex([1,1,0]), Vertex([2,2,0])])._valence_spectrum()
           defaultdict(<type 'list'>, {3: [0, 1]})

           >>> Fatgraph([Vertex([3, 1, 0, 1]), \
                      Vertex([4, 4, 0]), Vertex([3, 2, 2])])._valence_spectrum()
           defaultdict(<type 'list'>, {3: [1, 2], 4: [0]})
        """
        result = defaultdict(list)
        for (index, vertex) in enumerate(self.vertices):
            l = len(vertex)
            result[l].append(index)
        # consistency checks
        assert set(result.keys()) == set(self.vertex_valences()), \
            "Fatgraph.valence_spectrum:" \
            "Computed valence spectrum `%s` does not exhaust all " \
            " vertex valences %s" \
            % (result, self.vertex_valences())
        assert set(concat(result.values())) \
               == set(range(self.num_vertices)), \
            "Fatgraph.valence_spectrum:" \
            "Computed valence spectrum `%s` does not exhaust all " \
            " %d vertex indices" % (result, self.num_vertices)
        return result

    @ocache0
    def _starting_vertices(self):
        """
        Return the pair `(valence, vertices)`, which minimizes the
        product of valence with the number of vertices of that
        valence.  The `vertices` part of the pair is the list of
        indices in `self.vertices` corresponding to vertices with the
        chosen valence.

        Examples::

          >>> g = Fatgraph([Vertex([0,1,2]), Vertex([0,2,3]),
          ...              Vertex([3,4,4,5,5]), Vertex([1,6,6])])
          >>> g._starting_vertices()
          (5, [2])
        """
        val = max(self.vertex_valences())
        vs = None
        n = len(self.vertices) + 1
        for (val_, vs_) in self._valence_spectrum().iteritems():
            n_ = len(vs_)
            if (n_ * val_ < n * val) \
                    or (n_ * val_ == n * val and val_ < val):
                val = val_
                vs = vs_
                n = n_
        return (val, vs)

    @staticmethod
    def _compatible(v1, v2):
        """Return `True` if vertices `v1` and `v2` are compatible.
        (i.e., same valence and number of loops - one *could* be
        mapped onto the other.)
        """
        if len(v1) == len(v2) and v1.num_loops == v2.num_loops:
            return True
        else:
            return False

    @staticmethod
    def _compatible_vertices(v, g, ixs):
        """Iterate over all (indices of) vertices in `g`, which
        `v` *could* be mapped to (that is, the destination vertex
        matches `v` in valence and number of loops.

        Third argument `ixs` restricts the search to the given
        subset of vertex indices in `g`.
        """
        for i in ixs:
            if Fatgraph._compatible(v, g.vertices[i]):
                yield i

    class _CannotExtendMap(Exception):
        """Exception raised by `extend_map` on failure to extend a
        partial map.
        """
        pass

    @staticmethod
    def _extend_map(pv, rots, pe, G1, i1, r, G2, i2):
        """Extend map `(pv, rots, pe)` by mapping the `i1`-th
        vertex in `G1` to the `i2`-th vertex in `G2` (and rotating
        the source vertex by `r` places leftwards).  Return the
        extended map `(pv, rot, pe)`.

        The partial map is a triple `(pv, rot, pe)` as in
        `Fatgraph.isomorphism` (which see), with the additional
        proviso that unassigned items in `rot` are represented by
        `None`.
        """
        v1 = G1.vertices[i1]
        v2 = G2.vertices[i2]
        if not Fatgraph._compatible(v1, v2):
            raise Fatgraph._CannotExtendMap

        # XXX: rotation has to be >=0 for the [r:r+..] shift below to work
        if r < 0:
            r += len(v2)

        if i1 in pv:
            if pv[i1] != i2 or (rots[i1] - r) % len(v2) != 0:
                raise Fatgraph._CannotExtendMap
            else:
                # this pair has already been added
                return (pv, rots, pe)

        i = pv.find(i2)
        if i is not None:  # i.e., i in pv.values()
            if i != i1 or (rots[i] - r) % len(v2) != 0:
                raise Fatgraph._CannotExtendMap
            else:
                # this pair has already been added
                return (pv, rots, pe)

        # rotating `v1` leftwards is equivalent to rotating `v2` rightwards...
        v2 = v2[r:r + len(v2)]
        if not pe.extend(v1, v2):
            raise Fatgraph._CannotExtendMap

        pv[i1] = i2
        rots[i1] = r
        return (pv, rots, pe)

    @staticmethod
    def _neighbors(pv, pe, G1, v1, G2, v2):
        """List of vertex-to-vertex mappings that extend map `pv`
        in the neighborhood of vertices `v1` (in the domain) and
        `v2` (in the codomain).

        Return a list of triplets `(src, dst, rot)`, where:
           * `src` is the index of a vertex in `G1`,
             connected to `v1` by an edge `x`;
           * `dst` is the index of a vertex in `G2`,
             connected to `v2` by the image (according to `pe`)
             of edge `x`;
           * `rot` is the rotation to be applied to `G1[src]`
             so that edge `x` and its image appear
             at the same index position;
        """
        assert v2 == pv[v1]
        result = []
        for x in G1.vertices[v1]:
            if G1.edges[x].is_loop():
                continue  # with next edge `x`
            ((s1, a1), (s2, a2)) = G1.edges[x].endpoints
            src_v, src_i = (s2, a2) if (s1 == v1) else (s1, a1)
            # ignore vertices that are already in the domain of `m`
            if src_v in pv.keys():
                continue  # to next `x`
            ((d1, b1), (d2, b2)) = G2.edges[pe[x]].endpoints
            dst_v, dst_i = (d2, b2) if (d1 == v2) else (d1, b1)
            if dst_v in pv.values():
                continue  # to next `x`
            # array of (source vertex index, dest vertex index, rotation)
            result.append((src_v, dst_v, dst_i - src_i))
        return result

    def num_automorphisms(self):
        """Return the cardinality of the automorphism group of this
        `Fatgraph` object.

        Examples::

          >>> Fatgraph([Vertex([0,1,2]), Vertex([0,2,1])]).num_automorphisms()
          6
          >>> Fatgraph([Vertex([0,1,1]), Vertex([0,2,2])]).num_automorphisms()
          2
        """
        return len(list(self.automorphisms()))

    @ocache0
    def vertex_valences(self):
        return frozenset(len(v) for v in self.vertices)


def MgnTrivalentGraphsRecursiveGenerator(g, n):
    """Return a list of all connected trivalent fatgraphs having the
    prescribed genus `g` and number of boundary cycles `n`.
    
    Examples::

      >>> for g in MgnTrivalentGraphsRecursiveGenerator(0,3): print g
      Fatgraph([Vertex([1, 2, 1]), Vertex([2, 0, 0])]) 
      Fatgraph([Vertex([1, 0, 2]), Vertex([2, 0, 1])])

      >>> for g in MgnTrivalentGraphsRecursiveGenerator(1,1): print g
      Fatgraph([Vertex([1, 0, 2]), Vertex([2, 1, 0])])

    """
    # avoid infinite recursion in later statements
    if n == 0 or (g, n) < (0, 3):
        raise StopIteration

    # sanity check
    assert n > 0, \
        "MgnTrivalentGraphsRecursiveGenerator: " \
        " number of boundary cycles `n` must be positive," \
        " but got `%s` instead" % n

    logging.debug("Starting MgnTrivalentGraphsRecursiveGenerator(%d,%d) ..." % (g, n))

    ## M_{0,3} - induction base
    if (g, n) == (0, 3):
        yield Fatgraph([Vertex([1, 2, 1]), Vertex([2, 0, 0])])
        yield Fatgraph([Vertex([1, 0, 2]), Vertex([2, 0, 1])])
        logging.debug("  MgnTrivalentGraphsRecursiveGenerator(0,3) done.")

    ## M_{1,1} - induction base
    elif (g, n) == (1, 1):
        yield Fatgraph([Vertex([1, 0, 2]), Vertex([2, 1, 0])])
        logging.debug("  MgnTrivalentGraphsRecursiveGenerator(1,1) done.")

    ## General case
    else:
        try:
            checkpoint = os.path.join(runtime.options.checkpoint_dir,
                                      ("MgnTrivalentGraphsRecursiveGenerator%d,%d.list"
                                       % (g, n)))
        except AttributeError:
            # test run, `runtime.options` not defined
            checkpoint = None

        # try loading from file
        unique = None
        if checkpoint and runtime.options.restart:
            try:
                unique = load(checkpoint)
                if unique is not None:
                    logging.debug("  Loaded %d trivalent graphs from file '%s'",
                                  len(unique), checkpoint)
            except Exception, error:
                logging.debug("  Could not load saved state from file '%s': %s",
                              checkpoint, error.message)
                # ignore error and proceed

        if unique is None:
            # could not restore from saved state, have to compute
            unique = []
            discarded = 0
            timing.start("MgnTrivalentGraphsRecursiveGenerator(%d,%d)" % (g, n))
            for G in _MgnTrivalentGraphsRecursiveGenerator_main(g, n):
                # XXX: should this check be done in  *_main(g,n)?
                if (G.genus, G.num_boundary_cycles) != (g, n) or (G in unique):
                    discarded += 1
                    continue
                unique.append(G)
                yield G
            timing.stop("MgnTrivalentGraphsRecursiveGenerator(%d,%d)" % (g, n))
            logging.debug(
                "  MgnTrivalentGraphsRecursiveGenerator(%d,%d) done: %d unique graphs, discarded %d duplicates. (Elapsed: %.3fs)",
                g, n, len(unique), discarded,
                timing.get("MgnTrivalentGraphsRecursiveGenerator(%d,%d)" % (g, n)))
            # clear cached isomorphisms, they will not be needed in the sequel
            for G in unique:
                try:
                    G._cache_isomorphisms.clear()  # XXX: private impl. detail!
                except AttributeError:
                    pass
            # save to checkpoint file (if defined)
            if checkpoint:
                save(unique, checkpoint)
        else:
            # re-use loaded graphs
            for G in unique:
                yield G


def _MgnTrivalentGraphsRecursiveGenerator_main(g, n):
    logging.debug("  MgnTrivalentGraphsRecursiveGenerator(%d,%d): "
                  "pass 1: hang a circle to all edges of graphs in M_{%d,%d} ..." % (g, n, g, n - 1))
    for G in MgnTrivalentGraphsRecursiveGenerator(g, n - 1):
        for x in G.edge_orbits():
            yield G.hangcircle(x, 0)
            yield G.hangcircle(x, 1)

    logging.debug("  MgnTrivalentGraphsRecursiveGenerator(%d,%d): "
                  "pass 2: bridge all edges of a single graph in M_{%d,%d} ..." % (g, n, g, n - 1))
    for G in MgnTrivalentGraphsRecursiveGenerator(g, n - 1):
        for (x, y) in G.edge_pair_orbits():
            yield G.bridge(x, 0, y, 0)
            yield G.bridge(x, 0, y, 1)
            yield G.bridge(x, 1, y, 0)
            yield G.bridge(x, 1, y, 1)

    logging.debug("  MgnTrivalentGraphsRecursiveGenerator(%d,%d): "
                  "pass 3: bridge all edges of a single graph in M_{%d,%d} ..." % (g, n, g - 1, n + 1))
    for G in MgnTrivalentGraphsRecursiveGenerator(g - 1, n + 1):
        for (x, y) in G.edge_pair_orbits():
            yield G.bridge(x, 0, y, 0)
            yield G.bridge(x, 0, y, 1)
            yield G.bridge(x, 1, y, 0)
            yield G.bridge(x, 1, y, 1)

    ## logging.debug("  MgnTrivalentGraphsRecursiveGenerator(%d,%d): "
    ##               "pass 4: bridge two graphs of such that g_1+g_2=%d, n_1+n_2=%d ..." % (g,n, g,n+1)) 
    ## def add_up_to(x, min=0):
    ##     if x == 0 and min == 0:
    ##         yield (0,0)
    ##     elif x-min >= 0:
    ##         for y in xrange(min, x-min+1):
    ##             yield (y, x-y)
    ## for (g1, g2) in add_up_to(g, min=0):
    ##     for (n1, n2) in add_up_to(n+1, min=1):
    ##         if (g1, n1) < (0, 3) or (g2, n2) < (0,3):
    ##             continue
    ##         for G1 in MgnTrivalentGraphsRecursiveGenerator(g1,n1):
    ##             for G2 in MgnTrivalentGraphsRecursiveGenerator(g2,n2):
    ##                 for x1 in G1.edge_orbits():
    ##                     for x2 in G2.edge_orbits():
    ##                         yield Fatgraph.bridge2(G1, x1, 0, G2, x2, 0)
    ##                         yield Fatgraph.bridge2(G1, x1, 0, G2, x2, 1)
    ##                         yield Fatgraph.bridge2(G1, x1, 1, G2, x2, 0)
    ##                         yield Fatgraph.bridge2(G1, x1, 1, G2, x2, 1)


class MgnGraphsIterator(BufferingIterator):
    """Iterate over all connected fatgraphs having the
    prescribed genus `g` and number of boundary cycles `n`.
    
    Examples::

      >>> for g in MgnGraphsIterator(0,3): print g
      Fatgraph([Vertex([1, 2, 1]), Vertex([2, 0, 0])]) 
      Fatgraph([Vertex([1, 0, 2]), Vertex([2, 0, 1])])
      Fatgraph([Vertex([1, 1, 0, 0])])

      >>> for g in MgnGraphsIterator(1,1): print g
      Fatgraph([Vertex([1, 0, 2]), Vertex([2, 1, 0])])
      Fatgraph([Vertex([1, 0, 1, 0])])

    """

    def __init__(self, g, n):
        assert n > 0, \
            "MgnGraphsIterator: " \
            " number of boundary cycles `n` must be positive," \
            " but got `%s` instead" % n
        assert (g > 0) or (g == 0 and n >= 3), \
            "MgnGraphsIterator: " \
            " Invalid (g,n) pair (%d,%d): " \
            " need either g>0 or g==0 and n>2" \
            % (g, n)

        #: Prescribed genus of returned graphs
        self.g = g

        #: Prescribed number of boundary components
        self.n = n

        # Gather all 3-valent graphs.
        timing.start("MgnGraphsIterator: trivalent graphs")
        trivalent = list(MgnTrivalentGraphsRecursiveGenerator(g, n))
        timing.stop("MgnGraphsIterator: trivalent graphs")

        #: Fatgraphs to be contracted at next `.refill()` invocation
        self._batch = trivalent

        #: Graphs returned by next `.refill()` call will have this
        #  number of vertices.
        self._num_vertices = 4 * g + 2 * n - 5

        # initialize superclass with list of trivalent graphs
        BufferingIterator.__init__(self, trivalent)
        logging.info("  Found %d distinct unique trivalent fatgraphs. (Elapsed: %.3fs)",
                     len(trivalent), timing.get("MgnGraphsIterator: trivalent graphs"))

    def refill(self):
        if self._num_vertices == 0:
            raise StopIteration

        try:
            checkpoint = os.path.join(runtime.options.checkpoint_dir,
                                      ("M%d,%d-MgnGraphsIterator%d.list"
                                       % (runtime.g, runtime.n, self._num_vertices)))
        except AttributeError:
            # running a test, so no `runtime.options` defined
            checkpoint = None

        # try loading `next_batch` from old persisted state
        next_batch = None
        if checkpoint and runtime.options.restart:
            try:
                next_batch = load(checkpoint)
                if next_batch is not None:
                    logging.info("  Loaded graphs with %d vertices from file '%s'",
                                 self._num_vertices, checkpoint)
            except Exception, error:
                logging.debug("  Could not retrieve state from file '%s': %s",
                              checkpoint, error.message)
                # ignore error, continue generating graphs

        if next_batch is None:
            # really compute `next_batch`
            logging.debug("Generating graphs with %d vertices ...",
                          self._num_vertices)
            discarded = 0
            next_batch = []
            timing.start("MgnGraphsIterator: %d vertices" % self._num_vertices)
            for graph in self._batch:
                # contract all edges
                for edge in graph.edge_orbits():
                    if not graph.is_loop(edge):
                        dg = graph.contract(edge)
                        if dg not in next_batch:
                            # put graph back into next batch for processing
                            next_batch.append(dg)
                        else:
                            discarded += 1
            timing.stop("MgnGraphsIterator: %d vertices" % self._num_vertices)
            logging.info(
                "  Found %d distinct unique fatgraphs with %d vertices, discarded %d duplicates. (Elapsed: %.3fs)",
                len(self._batch), self._num_vertices, discarded,
                timing.get("MgnGraphsIterator: %d vertices" % self._num_vertices))
            if checkpoint is not None:
                save(next_batch, checkpoint)

        self._batch = next_batch
        self._num_vertices -= 1
        return next_batch


## main: run tests

if "__main__" == __name__:
    import doctest

    doctest.testmod(name='rg',
                    optionflags=doctest.NORMALIZE_WHITESPACE)
