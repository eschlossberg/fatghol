#! /usr/bin/env python
#
"""Classes and functions to deal with fatgraphs.
"""
__docformat__ = 'reStructuredText'

## logging subsystem

import logging

## stdlib imports

import debug, sys
from copy import copy
import operator
from itertools import chain,count,izip
import types


## application-local imports

from cache import (
    cache,
    cache1,
    cache_iterator,
    cache_symmetric,
    Cacheable,
    )
from combinatorics import (
    InplacePermutationIterator,
    SetProductIterator,
    Permutation,
    PermutationIterator,
    )
from cyclicseq import CyclicList,CyclicTuple
from iterators import (
    BufferingIterator,
    Iterator,
    itranslate,
    )
import persist
from utils import (
    concat,
    sign,
    )


## main

class VertexCache(object):
    """A caching factory of `Vertex` objects.
    """
    def __init__(self):
        self.cache = {}
    def __call__(self, edge_seq):
        key = tuple(edge_seq)
        if key not in self.cache:
            self.cache[key] = Vertex(key)
        return self.cache[key]
    def __str__(self):
        # needed to form readable persistent iterators cache
        return "rg.VertexCache"


class Vertex(Cacheable, CyclicList):
    """A (representative of) a vertex of a ribbon graph.

    A vertex is represented by the cyclically ordered list of its
    (decorated) edges.  The edge colorings may be accessed through a
    (read-only) sequence interface.
    """
    # *Note:* `Vertex` cannot be a `tuple` subclass because:
    #   1) `tuple` has no `index()` method and re-implementing one in
    #      pure Python would be less efficient;
    #   2) we could not implement `rotate()` and friends: tuples are
    #      immutable.

    def __init__(self, seq=None):
        Cacheable.__init__(self)
        CyclicList.__init__(self, seq)
        
    def __cmp__(self, other):
        """Return negative if x<y, zero if x==y, positive if x>y.
        Unlike standard Python sequence comparison, vertices with
        lower valence come first, and two vertices are only compared
        lexicographically if they have the same valence::

          >>> cmp(Vertex([0,1,2]), Vertex([0,1,2,3,4]))
          -1
          >>> cmp(Vertex([0,1,2,3,4]), Vertex([0,1,2]))
          1
          
          >>> cmp(Vertex([0,1,2]), Vertex([0,1,2]))
          0
          
          >>> cmp(Vertex([0,1,2,3]), Vertex([0,0,1,1]))
          1
          >>> cmp(Vertex([0,0,1,1]), Vertex([0,1,2,3]))
          -1
          
        """
        result = cmp(len(self), len(other))
        if 0 == result:
            if super(Vertex, self).__eq__(other):
                return 0
            else:
                if super(Vertex, self).__lt__(other):
                    return -1
                else:
                    return +1
        return result

    def __str__(self):
        return repr(self)
    
    @cache1
    def num_loops(self):
        """Return the number of loops attached to this vertex."""
        seen = {}
        loops = 0
        for x in xrange(len(self)):
            if self[x] in seen:
                loops += 1
            else:
                seen[self[x]] = True
        return loops


class EqualIfIsomorphic(object):
    """Instances of this class will compare equal if there is an
    isomorphism mapping one to the other.
    """
    
    @cache_symmetric
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



class Fatgraph(EqualIfIsomorphic, Cacheable):
    """A fully-decorated ribbon graph.

    Exports a (read-only) sequence interface, through which vertices
    can be accessed.

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

    def __init__(self, g_or_vs, vertextype=Vertex, **kwargs):
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
              >>> g2.endpoints_v is g1.endpoints_v
              True
              >>> g2.endpoints_i is g1.endpoints_i
              True

        """
        # set this instance's _persistent_id
        Cacheable.__init__(self)

        # dispatch based on type of arguments passed
        if isinstance(g_or_vs, Fatgraph):
            # copy-constructor
            self._vertextype = vertextype
            self.edge_numbering = g_or_vs.edge_numbering
            self.endpoints_i = g_or_vs.endpoints_i
            self.endpoints_v = g_or_vs.endpoints_v
            self.num_edges = g_or_vs.num_edges
            self.num_external_edges = g_or_vs.num_external_edges
            self.num_vertices = g_or_vs.num_vertices
            self.vertices = g_or_vs.vertices

        else:
            # initialize new instance
            assert debug.is_sequence_of_type(Vertex, g_or_vs), \
                   "Fatgraph.__init__: parameter `g_or_vs` must be" \
                   " sequence of `%s` instances;" \
                   " got `%s` instead, which has type `%s`." \
                   % (Vertex, g_or_vs, [type(x) for x in g_or_vs])

            #: Factory method to make a `Vertex` instance from a linear
            #  list of incident edge colorings.
            self._vertextype = vertextype

            #: list of vertices
            self.vertices = g_or_vs

            #: Number of edge colors
            self.num_edges = kwargs.get('num_edges',
                                        sum(len(v) for v in self.vertices) / 2)

            #: Number of external (loose-end) edges
            self.num_external_edges = kwargs.get('num_external_edges', 0)

            #: Number of vertices  XXX: why is this settable with kwarg???
            self.num_vertices = kwargs.get('num_vertices', len(g_or_vs))

            if 'endpoints' in kwargs:
                (self.endpoints_v, self.endpoints_i) = kwargs.get('endpoints')
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
                self.endpoints_v = [ [] for dummy in xrange(self.num_edges) ]
                self.endpoints_i = [ [] for dummy in xrange(self.num_edges) ]
                for current_vertex_index in xrange(self.num_vertices):
                    for (edge_index_in_vertex, edge) \
                            in enumerate(self.vertices[current_vertex_index]):
                        assert edge in range(self.num_edges), \
                                   "Fatgraph.__init__:"\
                                   " edge number %d not in range 0..%d" \
                                   % (edge, self.num_edges)
                        self.endpoints_v[edge].append(current_vertex_index)
                        self.endpoints_i[edge].append(edge_index_in_vertex)

            ## Orientation is given by an ordering of the edges,
            ## which directly translates into an orientation of the
            ## associated cell.  
            if 'orientation' in kwargs:
                self.edge_numbering = kwargs.get('orientation')
            else:
                self.edge_numbering = [ x for x in xrange(self.num_edges) ]

        # before computing invariants, check that internal data
        # structures are in a consistent state
        assert self._ok()

        # used for isomorphism testing
        self.invariants = (
            self.num_vertices,
            self.num_edges,
            self.num_external_edges if self.num_external_edges > 0
                                    else self.num_boundary_cycles(),
            #self.vertex_valences(),
            #self.vertex_valence_distribution(),
            )


    def _ok(self):
        """Perform coherency checks on internal state variables of
        `Fatgraph` instance and return `True` if they all pass.
        """
        assert self.num_edges > 0, \
               "Fatgraph `%s` has 0 edges." % (self)
        # check regular edges endpoints
        for (edge, ep_v, ep_i) in izip(count(),
                                       self.endpoints_v[:self.num_edges],
                                       self.endpoints_i[:self.num_edges]):
            assert isinstance(ep_v, list)  # Fatgraph.contract() updates this destructively
            assert isinstance(ep_i, list)  # Fatgraph.contract() updates this destructively
            assert len(ep_v) == 2
            assert len(ep_i) == 2
            assert isinstance(ep_v[0], int)
            assert isinstance(ep_v[1], int)
            assert isinstance(ep_i[0], int)
            assert isinstance(ep_i[1], int)
            assert (0 <= ep_v[0] < self.num_vertices)
            assert (0 <= ep_v[1] < self.num_vertices)
            assert (0 <= ep_i[0] < len(self.vertices[ep_v[0]]))
            assert (0 <= ep_i[1] < len(self.vertices[ep_v[1]])), \
                   "Fatgraph `%s`:"\
                   " invalid attachment indices `%s`" \
                   " for endpoints %s of regular edge %d"\
                   % (self, ep_i, ep_v, edge)
##                    "Fatgraph `%s` has invalid regular endpoints array `%s/%s`" \
##                    " invalid endpoints pair %s/%s for edge %d" \
##                    % (self, self.endpoints_v, self.endpoints_i,
##                       ep_v, ep_i, edge)
            assert (edge in self.vertices[ep_v[0]]), \
                    "Invalid endpoints %s for edge %d of graph `%s`" \
                    % (ep_v, edge, self)
            assert (edge in self.vertices[ep_v[1]]), \
                    "Invalid endpoints %s for edge %d of graph `%s`" \
                    % (ep_v, edge, self)
        # check external edges endpoints
        for (edge, ep_v, ep_i) in izip(count(),
                                       self.endpoints_v[self.num_edges:],
                                       self.endpoints_i[self.num_edges:]):
            assert isinstance(ep_v, list)  # Fatgraph.contract() updates this destructively
            assert isinstance(ep_i, list)  # Fatgraph.contract() updates this destructively
            xedge = -self.num_external_edges + edge
            assert (ep_v[1] is None)
            assert isinstance(ep_v[0], int)
            assert (ep_i[1] is None)
            assert isinstance(ep_i[0], int)
            assert (0 <= ep_v[0] < self.num_vertices)
            assert (0 <= ep_i[0] < len(self.vertices[ep_v[0]]))
##                    "Fatgraph `%s` has invalid external endpoints array: `%s/%s`" \
##                    % (self, self.endpoints_v, self.endpoints_i)
            assert (xedge in self.vertices[ep_v[0]]), \
                   "Invalid endpoints %s for external edge %d of graph `%s`" \
                   % (ep, xedge, self)
        # check that each edge occurs exactly two times in vertices
        cnt = [ 0 for x in xrange(self.num_edges + self.num_external_edges) ]
        for v in self.vertices:
            for edge in v:
                cnt[edge] += 1
        for edge, cnt in enumerate(cnt):
            if edge < self.num_edges:
                assert cnt == 2, \
                       "Regular edge %d appears in %d vertices" \
                       % (edge, cnt)
            else:
                assert cnt == 1, \
                       "External edge %d appears in %d vertices" \
                       % (edge, cnt)

        assert self.edge_numbering is not None
        
        return True


    def __getitem__(self, index):
        return self.vertices[index]


    def __iter__(self):
        """Return iterator over vertices."""
        return iter(self.vertices)


    def __repr__(self):
        if hasattr(self, 'num_external_edges') and self.num_external_edges > 0:
            return "Fatgraph(%s, num_external_edges=%d)" \
                   % (repr(self.vertices), self.num_external_edges)
        elif hasattr(self, 'vertices'):
            return "Fatgraph(%s)" % repr(self.vertices)
        else:
            return "Fatgraph(<Initializing...>)"

    
    def __str__(self):
        return repr(self)


    def automorphisms(self):
        """Enumerate automorphisms of this `Fatgraph` object.

        See `.isomorphisms()` for details of how a `Fatgraph`
        isomorphism is represented.
        """
        return self.isomorphisms(self)


    class BoundaryCycle(frozenset):
        """A boundary cycle of a Fatgraph.

        Boundary cycles are a cyclic sequence of 'corners': a corner
        consists of a vertex `v` and (an unordered pair of) two
        consecutive indices (in the cyclic order at `v`, so, either `j
        == i+1` or `i` and `j` are the starting and ending indices).

        Two boundary cycles are equal if they comprise the same
        corners.
        """
        def __init__(self, triples, graph=None):
            """Construct a `BoundaryCycle` instance from a sequence of
            triples `(v, i, j)`, where `i` and `j` are consecutive (in
            the cyclic order sense) indices at a vertex `v`.
            """
            self.origin = graph
            frozenset.__init__(self, triples)
            if __debug__:
                if graph is not None:
                    for (v, i, j) in self:
                        l = len(graph.vertices[v])-1
                        assert (abs(i-j) == 1) or (set((i,j)) == set((0,l))), \
                               "Fatgraph.BoundaryCycle():" \
                               " Non-consecutive indices in triple `%s`" \
                               % ((v,i,j),)

        def transform(self, iso):
            """Return a new `BoundaryCycle` instance, obtained by
            transforming each corner according to a graph isomorphism.
            """
            assert self.origin is not None
            (pv, rots, pe) = iso
            triples = []
            for (v, i, j) in self:
                l = len(self.origin.vertices[v])
                # create transformed triple 
                v_ = pv[v]
                i_ = (i + rots[v]) % l # XXX: is it `-` or `+`?
                j_ = (j + rots[v]) % l
                # ensure the contract is honored, that `j` is the
                # index _following_ `i` in the cyclic order
                if i_ == 0 and j_ == l:
                    i_, j_ = j_, i_
                triples.append((v_, i_, j_))
            return Fatgraph.BoundaryCycle(triples)
                

    @cache
    def boundary_cycles(self):
        """Return a list of boundary cycles of this `Fatgraph` object.

        Boundary cycles are represented as a cyclic list of 'corners':
        a corner is a triple `(v, i, j)` consisting of a vertex and
        two consecutive indices (in the cyclic order, so, either `j ==
        i+1` or `i` and `j` are the starting and ending indices)::
        
          >>> Fatgraph([Vertex([2,1,0]),Vertex([2,0,1])]).boundary_cycles()
          [BoundaryCycle([(1, 0, 1), (0, 2, 0)]),
           BoundaryCycle([(1, 1, 2), (0, 1, 2)]),
           BoundaryCycle([(1, 2, 0), (0, 0, 1)])]

        This verbose representation allows one to distinguish the
        boundary cycles made from the same set of edges::

          >>> Fatgraph([Vertex([0,1,2,0,1,2])]).boundary_cycles()
          [BoundaryCycle([(0, 0, 1), (0, 2, 3), (0, 4, 5)]),
           BoundaryCycle([(0, 0, 5), (0, 3, 4), (0, 1, 2)])]
        """
        assert self.num_external_edges == 0, \
               "Fatgraph.boundary_cycles: "\
               " cannot compute boundary cycles for" \
               " a graph with nonzero external edges: %s" % self

        # auxiliary function
        def other_end(graph, edge, vertex, attachment):
            """Return the other endpoint of `edge` as a pair `(v, i)`.
            """
            ends = zip(graph.endpoints_v[edge], graph.endpoints_i[edge])
            if ends[0] == (vertex, attachment):
                return ends[1]
            else:
                return ends[0]

        # Build the collection of "corners" of `graph`,
        # structured just like the set of vertices.
        # By construction, `corners[v][i]` has the the
        # form `(v,i,j)` where `j` is the index following
        # `i` in the cyclic order.
        corners = [ [ (v, i, (i+1)%len(self.vertices[v])) for i in xrange(len(self.vertices[v])) ]
                    for v in xrange(self.num_vertices) ]

        result = []
        while True:
            # fast-forward to the first unused corner
            for v in xrange(self.num_vertices):
                for i in xrange(len(self.vertices[v])):
                    if corners[v][i] is not None:
                        break
                if corners[v][i] is not None:
                    break
            if corners[v][i] is None:
                break

            # build a list of corners comprising the same boundary
            # cycle: start with one corner, follow the edge starting
            # at the second delimiter of the corner to its other
            # endpoint, and repeat until we come back to the starting
            # point.  
            corner = None
            start = (v,i)
            triples = []
            while (v,i) != start or len(triples) == 0:
                corner = corners[v][i]
                corners[v][i] = None
                triples.append(corner)
                assert v == corner[0]
                assert i == corner[1]
                j = corner[2]
                edge = self.vertices[v][j]
                (v,i) = other_end(self, edge, v, j)
            result.append(Fatgraph.BoundaryCycle(triples, graph=self))

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
        assert side1 in [0,1], \
               "Fatgraph.bridge: Invalid value for `side1`: '%s' - should be 0 or 1" % side1
        assert side2 in [0,1], \
               "Fatgraph.bridge: Invalid value for `side2`: '%s' - should be 0 or 1" % side2
        
        opposite_side1 = 0 if side1 else 1
        opposite_side2 = 0 if side2 else 1

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

        if side1:
            midpoint1 = self._vertextype([other_half1, one_half1, connecting_edge])
        else:
            midpoint1 = self._vertextype([one_half1, other_half1, connecting_edge])

        if side2:
            midpoint2 = self._vertextype([other_half2, one_half2, connecting_edge])
        else:
            midpoint2 = self._vertextype([one_half2, other_half2, connecting_edge])

        ## two new vertices are added: the mid-points of the connected edges.
        new_vertices = self.vertices + [midpoint1, midpoint2]
        ## the connecting edge has endpoints in the mid-points of
        ## `edge1` and `edge2`, and is *always* in third position.
        new_endpoints_v = self.endpoints_v + [ [midpoint1_index, midpoint2_index] ]
        new_endpoints_i = self.endpoints_i + [ [2,2] ]
        
        (v1a, v1b) = self.endpoints_v[edge1]
        (pos1a, pos1b) = self.endpoints_i[edge1]
        new_endpoints_v[one_half1] = [v1a, midpoint1_index]
        new_endpoints_i[one_half1] = [pos1a, side1]
        if edge1 != edge2:
            # replace `edge1` with new `other_half1` in the second endpoint
            new_vertices[v1b] = self._vertextype(new_vertices[v1b][:pos1b]
                                                 + [other_half1]
                                                 + new_vertices[v1b][pos1b+1:])
            new_endpoints_v.append([midpoint1_index, v1b])  # other_half1
            new_endpoints_i.append([opposite_side1, pos1b]) # other_half1
        else:
            # same edge, "other half" ends at the second endpoint
            new_endpoints_v.append([midpoint1_index, midpoint2_index])
            new_endpoints_i.append([opposite_side1, side2]) 

        # replace `edge2` with new `other_half2` in the second
        # endpoint; again we need to distinguish the special case when
        # `edge1` and `edge2` are the same edge.
        (v2a, v2b) = self.endpoints_v[edge2]
        (pos2a, pos2b) = self.endpoints_i[edge2]
        if edge1 != edge2:
            new_endpoints_v[one_half2] = [v2a, midpoint2_index]
            new_endpoints_i[one_half2] = [pos2a, side2]
        else:
            # `edge1 == edge2`, so `one_half2 == other_half1`
            new_endpoints_v[one_half2] = [midpoint1_index, midpoint2_index]
            new_endpoints_i[one_half2] = [opposite_side1, side2]
        # "other half" of second edge *always* ends at the previous
        # edge endpoint, so replace `edge2` in `v2b`.
        new_vertices[v2b] = self._vertextype(new_vertices[v2b][:pos2b]
                                             + [other_half2]
                                             + new_vertices[v2b][pos2b+1:])
        new_endpoints_v.append([midpoint2_index, v2b])  # other_half2
        new_endpoints_i.append([opposite_side2, pos2b]) # other_half2

        # build new graph 
        new_edge_numbering = self.edge_numbering + \
                             [other_half1, other_half2, connecting_edge]
        return Fatgraph(new_vertices,
                     vertextype = self._vertextype,
                     endpoints = (new_endpoints_v, new_endpoints_i),
                     num_edges = self.num_edges + 3,
                     num_external_edges = self.num_external_edges,
                     orientation = new_edge_numbering,
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
        assert side1 in [0,1], \
               "Fatgraph.bridge2: Invalid value for `side1`: '%s' -- should be 0 or 1" % side1
        assert side2 in [0,1], \
               "Fatgraph.bridge2: Invalid value for `side2`: '%s' -- should be 0 or 1" % side2

        ## First, build a (non-connected) graph from the disjoint
        ## union of `self` and `other`.

        # Edges of `other` are renumbered depending on whether
        # they are internal of external edges:
        #   - internal edges in `other` have numbers ranging from 0 to
        #     `other.num_edges`: they get new numbers starting from
        #     `self.num_edges` and counting upwards
        renumber_other_edges = dict((x,x+self.num_edges)
                                    for x in xrange(other.num_edges))
        #   - external edges in `other` have *negative* indices: they
        #     are renumbered starting from `self.num_external_edges-1`
        #     and counting downwards.
        renumber_other_edges.update((x,-self.num_external_edges+x)
                                    for x in xrange(other.num_external_edges))
        # Orientation needs the same numbering:
        new_edge_numbering = self.edge_numbering \
                             + list(itranslate(renumber_other_edges, other.edge_numbering))
        # Similarly, vertices of `self` retain indices `[0..v]`, while
        # vertices of `other` follow.
        new_vertices = self.vertices \
                       + [ self._vertextype(itranslate(renumber_other_edges, ov))
                           for ov in other.vertices ]
        renumber_other_vertices = dict((x, x+self.num_vertices)
                                       for x in xrange(other.num_vertices))
        # vertex indices need to be shifted for endpoints
        new_endpoints_v = self.endpoints_v \
                          + [ [renumber_other_vertices[x], renumber_other_vertices[y]]
                              for (x,y) in other.endpoints_v ]
        # but vertex positions are the same
        new_endpoints_i = self.endpoints_i + other.endpoints_i

        # FIXME: From this point onwards, the code basically is the
        # same as in `Fatgraph.bridge`, copied and edited here for
        # efficiency reasons.

        edge2 = renumber_other_edges[edge2] # need new number
        
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
            midpoint1 = self._vertextype([other_half1, one_half1, connecting_edge])
        else:
            midpoint1 = self._vertextype([one_half1, other_half1, connecting_edge])

        if side2:
            midpoint2 = self._vertextype([other_half2, one_half2, connecting_edge])
        else:
            midpoint2 = self._vertextype([one_half2, other_half2, connecting_edge])

        ## two new vertices are added: the mid-points of the connected edges.
        new_vertices += [midpoint1, midpoint2]
        ## the connecting edge has endpoints in the mid-points of
        ## `edge1` and `edge2`, and is *always* in third position.
        new_endpoints_v += [[midpoint1_index, midpoint2_index]]
        new_endpoints_i += [[2,2]]
        
        (v1a, v1b) = new_endpoints_v[edge1]
        (pos1a, pos1b) = new_endpoints_i[edge1]
        new_endpoints_v[one_half1] = [v1a, midpoint1_index]
        new_endpoints_i[one_half1] = [pos1a, side1]
        # replace `edge1` with new `other_half1` in the second endpoint
        new_vertices[v1b] = self._vertextype(new_vertices[v1b][:pos1b]
                                             + [other_half1]
                                             + new_vertices[v1b][pos1b+1:])
        new_endpoints_v.append([midpoint1_index, v1b])  # other_half1
        new_endpoints_i.append([opposite_side1, pos1b]) # other_half1

        # replace `edge2` with new `other_half2` in the second
        # endpoint; again we need to distinguish the special case when
        # `edge1` and `edge2` are the same edge.
        (v2a, v2b) = new_endpoints_v[edge2]
        (pos2a, pos2b) = new_endpoints_i[edge2]
        new_endpoints_v[one_half2] = [v2a, midpoint2_index]
        new_endpoints_i[one_half2] = [pos2a, side2]
        # "other half" of second edge *always* ends at the previous
        # edge endpoint, so replace `edge2` in `v2b`.
        new_vertices[v2b] = self._vertextype(new_vertices[v2b][:pos2b]
                                             + [other_half2]
                                             + new_vertices[v2b][pos2b+1:])
        new_endpoints_v.append([midpoint2_index, v2b])  # other_half2
        new_endpoints_i.append([opposite_side2, pos2b]) # other_half2

        # build new graph 
        new_edge_numbering +=  [other_half1, other_half2, connecting_edge]
        return Fatgraph(new_vertices,
                     vertextype = self._vertextype,
                     endpoints = (new_endpoints_v, new_endpoints_i),
                     num_edges = self.num_edges + other.num_edges + 3,
                     num_external_edges = self.num_external_edges + other.num_external_edges,
                     orientation = new_edge_numbering,
                     )


    def _cmp_orient(self, other, iso):
        pe = iso[2]
        image_edge_numbering = Permutation(dict((self.edge_numbering[x],
                                                 other.edge_numbering[pe[x]])
                                                for x in xrange(self.num_edges)))
        return image_edge_numbering.sign()


    @cache
    def contract(self, edgeno):
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
        # check that we are not contracting a loop or an external edge
        assert not self.is_loop(edgeno), \
               "Fatgraph.contract: cannot contract a loop."
        assert (self.endpoints_v[edgeno][0] is not None) \
               and (self.endpoints_v[edgeno][1] is not None), \
               "Fatgraph.contract: cannot contract an external edge."
        assert (edgeno >= 0) and (edgeno < self.num_edges), \
               "Fatgraph.contract: invalid edge number (%d):"\
               " must be in range 0..%d" \
               % (edgeno, self.num_edges)

        ## Plug the higher-numbered vertex into the lower-numbered one.
        
        # store endpoints of the edge-to-be-contracted
        (v1, v2) = self.endpoints_v[edgeno]
        (pos1, pos2) = self.endpoints_i[edgeno]
        if v1 > v2:
            # swap endpoints so that `v1 < v2`
            v1, v2 = v2, v1
            pos1, pos2 = pos2, pos1

        # save highest-numbered index of vertices to be contracted
        l1 = len(self.vertices[v1]) - 1
        l2 = len(self.vertices[v2]) - 1

        # Build new list of vertices, removing the contracted edge and
        # shifting all indices above:
        #   - edges numbered 0..edgeno-1 are unchanged;
        #   - edges numbered `edgeno+1`.. are renumbered, 
        #     shifting the number down one position;
        #   - edge `edgeno` is kept intact, will be removed by mating
        #     operation (see below).
        renumber_edges = dict((i+1,i)
                              for i in xrange(edgeno, self.num_edges))
        # See `itranslate` in utils.py for how this prescription is
        # encoded in the `renumber_edges` mapping.
        new_vertices = [ self._vertextype(itranslate(renumber_edges, v))
                         for v in self.vertices ]

        # Mate endpoints of contracted edge:
        # 1. Rotate endpoints `v1`, `v2` so that the given edge would
        #    appear *last* in `v1` and *first* in `v2` (*Note:* since
        #    the contracted edge has already been deleted and `v1`,
        #    `v2` are *cyclic*, this means that we do the same
        #    operation on `v1` and `v2` alike).
        # 2. Join vertices by concatenating the list of incident
        #    edges;
        # 3. Set new `i1` vertex in place of old first endpoint:
        new_vertices[v1] = self._vertextype(
            new_vertices[v1][pos1+1:] + new_vertices[v1][:pos1]
            +
            new_vertices[v2][pos2+1:] + new_vertices[v2][:pos2]
            )
        # 4. Remove second endpoint from list of new vertices:
        del new_vertices[v2]

        # vertices with index above `v2` are now shifted down one place
        renumber_vertices = dict((i+1,i)
                                 for i in xrange(v2, self.num_vertices))
        # vertex `v2` is mapped to vertex `v1`
        renumber_vertices[v2] = v1
        new_endpoints_v = [ list(itranslate(renumber_vertices, ep))
                          for ep in  self.endpoints_v ]
        del new_endpoints_v[edgeno]

        new_endpoints_i = [ ep_i[:] for ep_i in self.endpoints_i ]
        del new_endpoints_i[edgeno]
        # renumber attachment indices, according to the mating of
        # vertices `v1` and `v2`:
        # - on former vertex `v1`:
        #   * indices (pos1+1)..l1 are mapped to 0..(l1-pos1-1)
        #     in the mated vertex;
        #   * index pos1 is deleted;
        #   * indices 0..pos1-1 are mapped to (l1-pos1)..l1-1;
        renumber_pos1 = { pos1:None }
        renumber_pos1.update(dict((pos1+1+x, x)
                             for x in xrange(l1-pos1)))
        renumber_pos1.update(dict((x, l1-pos1+x)
                             for x in xrange(pos1)))
        for edge in self.vertices[v1]:
            if edge == edgeno:
                continue # skip contracted edge
            if v1 == self.endpoints_v[edge][0]:
                new_endpoints_i[renumber_edges.get(edge, edge)][0] = renumber_pos1[self.endpoints_i[edge][0]]
            if v1 == self.endpoints_v[edge][1]:
                new_endpoints_i[renumber_edges.get(edge, edge)][1] = renumber_pos1[self.endpoints_i[edge][1]]
        # - on former vertex `v2`:
        #   * indices (pos2+1)..l2 are mapped to l1..(l1+l2-pos2-1);
        #   * index pos2 is deleted;
        #   * indices 0..pos2-1 are mapped to (l1+l2-pos2)..l1+l2-1:
        renumber_pos2 = { pos2:None }
        renumber_pos2.update(dict((pos2+1+x, l1+x)
                                  for x in xrange(l2-pos2)))
        renumber_pos2.update(dict((x, l1+l2-pos2+x)
                                  for x in xrange(pos2)))
        for edge in self.vertices[v2]:
            if edge == edgeno:
                continue # skip contracted edge
            if v2 == self.endpoints_v[edge][0]:
                new_endpoints_i[renumber_edges.get(edge, edge)][0] = renumber_pos2[self.endpoints_i[edge][0]]
            if v2 == self.endpoints_v[edge][1]:
                new_endpoints_i[renumber_edges.get(edge, edge)][1] = renumber_pos2[self.endpoints_i[edge][1]]

        ## Orientation of the contracted graph.

        cut = self.edge_numbering[edgeno]
        renumber_edge_numbering = { cut:None }
        renumber_edge_numbering.update(dict((x,x) for x in xrange(cut)))
        renumber_edge_numbering.update(dict((x+1,x)
                                      for x in xrange(cut,self.num_edges)))
        new_edge_numbering = [ renumber_edge_numbering[self.edge_numbering[x]]
                         for x in xrange(self.num_edges)
                         if x != edgeno ]
        
        # consistency check
        if __debug__:
            assert len(new_endpoints_v) == self.num_edges - 1
            assert len(new_endpoints_i) == len(new_endpoints_v)
            assert len(new_edge_numbering) == self.num_edges - 1
            for x in xrange(self.num_edges - 1):
                assert 0 <= new_edge_numbering[x] < self.num_edges - 1
            g = Fatgraph(new_vertices,
                         vertextype = self._vertextype,
                         endpoints = (new_endpoints_v, new_endpoints_i),
                         num_edges = self.num_edges - 1,
                         num_external_edges = self.num_external_edges,
                         orientation = new_edge_numbering,
                         )
            assert g.num_boundary_cycles() == self.num_boundary_cycles(), \
                   "Fatgraph.contract(%s, %d):" \
                   " Contracted graph `%s` does not have the same number" \
                   " of boundary cycles of parent graph." \
                   % (self, edgeno, g)

        # build new graph 
        return Fatgraph(new_vertices,
                     vertextype = self._vertextype,
                     endpoints = (new_endpoints_v, new_endpoints_i),
                     num_edges = self.num_edges - 1,
                     num_external_edges = self.num_external_edges,
                     orientation = new_edge_numbering,
                     )


    def edges(self):
        return xrange(self.num_edges)

    
    @cache
    def edge_orbits(self):
        """Compute orbits of the edges under the action of graph
        automorphism group, and a representative for each orbit.
        
        Returns a dictionary, whose keys are the representatives, and
        whose values are the orbits.

        Orbits are represented as Python `list` objects; the order the
        items appear in a list `L` is the order the orbit is swept by
        applying graph automorphisms to `L[0]`.

        Examples::

          >>> g = Fatgraph([Vertex([0,1,2]), Vertex([0,1,2])])
          >>> g.edge_orbits()
          { 0:[0,1,2] }
          
        """
        orbits = dict( (x,[x]) for x in xrange(self.num_edges) )
        seen = set()
        for a in self.automorphisms():
            edge_permutation = a[2]
            for x in xrange(self.num_edges):
                if x in seen:
                    continue
                y = edge_permutation[x]
                # `x` and `y` are in the same orbit, only keep the one
                # with lower abs. value, and remove the other.
                if x < y:
                    orbits[x].append(y)
                    if y not in seen:
                        del orbits[y]
                        seen.add(y)
                else: # x > y
                    continue
        return orbits


    @cache
    def edge_pair_orbits(self):
        """Compute orbits of pairs `(edge1, edge2)` under the action
        of graph automorphism group, and a representative for each
        orbit.
        
        Returns a dictionary, whose keys are the representatives, and
        whose values are the orbits.

        Orbits are represented as Python `list` objects; the order the
        items appear in a list `L` is the order the orbit is swept by
        applying graph automorphisms to `L[0]`.

        Examples::

          >>> g = Fatgraph([Vertex([0,1,2]), Vertex([0,1,2])])
          >>> g.edge_pair_orbits()
          { 0:[0,1,2] }
          
        """
        edge_pairs = [ (x,y) 
                       for x in xrange(self.num_edges)
                       for y in xrange(self.num_edges) ]
        orbits = dict( (p,[p]) for p in edge_pairs )
        seen = set()
        for a in self.automorphisms():
            edge_permutation = a[2]
            for p in edge_pairs:
                if p in seen:
                    continue
                q = (edge_permutation[p[0]], edge_permutation[p[1]])
                # `p` and `q` are in the same orbit, only keep the one
                # with lower abs. value, and remove the other.
                if p < q:
                    orbits[p].append(q)
                    if q not in seen:
                        del orbits[q]
                        seen.add(q)
                else: # x > y
                    continue
        return orbits


    def genus(self):
        """Return the genus g of this `Fatgraph` object."""
        n = self.num_boundary_cycles()
        K = self.num_vertices
        L = self.num_edges
        # by Euler, K-L+n=2-2*g
        return (L - K - n + 2) / 2


    def graft(self, G, v):
        """Return new `Fatgraph` formed by grafting graph `G` into vertex
        with index `v`.  The number of"external" edges in `G` must match the
        valence of `v`.
        """
        assert G.num_external_edges == len(self.vertices[v]), \
               "Fatgraph.graft:" \
               " attempt to graft %d-legged graph `%s`"\
               " into %d-valent vertex `%s`" \
               % (G.num_external_edges, G,
                  len(self.vertices[v]), self.vertices[v])
        vertextype = self._vertextype # micro-optimization

        # edges of `G` are renumbered depending on whether
        # they are internal of external edges:
        #   - internal edges in `G` have numbers ranging from 0 to
        #     `G.num_edges`: they get new numbers starting from
        #     `self.num_edges` and counting upwards
        renumber_g_edges = dict((x,x+self.num_edges)
                                for x in xrange(G.num_edges))
        #   - external edges in `G` are mated with edges incoming to
        #     vertex `v`: the first external edge (labeled -1)
        #     corresponds to the first edge in `v`, the second
        #     external edge (labeled -2) to the second edge in `v`,
        #     and so on.
        renumber_g_edges.update((-n-1,l)
                                for (n,l) in enumerate(self.vertices[v]))

        # the first `v-1` vertices of the new graph are the first
        # `v-1` vertices of `self`; then come vertices `v+1`,... of
        # `self`; vertices from `G` come last in the new graph
        new_vertices = (self.vertices[:v] 
                        + self.vertices[v+1:] 
                        + [ vertextype(itranslate(renumber_g_edges, gv))
                            for gv in G.vertices ])

        return Fatgraph(new_vertices, vertextype=vertextype,
                     num_edges = self.num_edges + G.num_edges,
                     num_external_edges = self.num_external_edges)


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
        assert side in [0,1], \
               "Fatgraph.hangcircle: Invalid value for `side`: '%s' - should be 0 or 1" % side1
        
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
        if side:
            midpoint = self._vertextype([other_half, one_half, connecting_edge])
        else:
            midpoint = self._vertextype([one_half, other_half, connecting_edge])
        T = self._vertextype([circling_edge, circling_edge, connecting_edge])
        new_vertices = self.vertices + [midpoint, T]

        ## new edge endpoints:
        ## - start with a copy of the original ednpoints:
        new_endpoints_v = copy(self.endpoints_v)
        new_endpoints_i = copy(self.endpoints_i)

        ## - replace `edge` with new `other_half` in the second endpoint:
        (v1, v2) = self.endpoints_v[edge]
        (pos1, pos2) = self.endpoints_i[edge]
        new_endpoints_v[one_half] = [v1, midpoint_index]
        new_endpoints_i[one_half] = [pos1, side]
        new_vertices[v2] = self._vertextype(new_vertices[v2][:pos2]
                                             + [other_half]
                                             + new_vertices[v2][pos2+1:])
        new_endpoints_v.append([midpoint_index, v2])  # other_half1
        new_endpoints_i.append([opposite_side, pos2]) # other_half1

        ## - the connecting edge has endpoints in the mid-point of
        ## `edge` and in `T`, and is *always* in third position:
        new_endpoints_v.append([midpoint_index, T_index])
        new_endpoints_i.append([2,2])

        ## - the circling edge is a loop with vertex `T`
        new_endpoints_v.append([T_index, T_index])
        new_endpoints_i.append([0,1])
        
        # finally, build new graph 
        new_edge_numbering = self.edge_numbering + \
                             [other_half, connecting_edge, circling_edge]
        return Fatgraph(new_vertices,
                     vertextype = self._vertextype,
                     endpoints = (new_endpoints_v, new_endpoints_i),
                     num_edges = self.num_edges + 3,
                     num_external_edges = self.num_external_edges,
                     orientation = new_edge_numbering,
                     )
    
        
    def is_connected(self):
        """Return `True` if graph is connected.

        Count all vertices that we can reach from the 0th vertex,
        using a breadth-first algorithm; the graph is connected iff
        this count equals the number of vertices.

        See:
          http://brpreiss.com/books/opus4/html/page554.html#SECTION0017320000000000000000
          http://brpreiss.com/books/opus4/html/page561.html#SECTION0017341000000000000000
          
        Examples::
          >>> Fatgraph([Vertex([3, 3, 0, 0]), Vertex([2, 2, 1, 1])]).is_connected()
          False
          >>> Fatgraph([Vertex([3, 1, 2, 0]), Vertex([3, 0, 2, 1])]).is_connected()
          True
        """
        endpoints_v = self.endpoints_v
        endpoints_i = self.endpoints_i
        visited_edges = set()
        visited_vertices = set()
        vertices_to_visit = [0]
        for vi in vertices_to_visit:
            # enqueue neighboring vertices that are not connected by
            # an already-visited edge
            for l in self.vertices[vi]:
                if l not in visited_edges:
                    # add other endpoint of this edge to the to-visit list
                    if endpoints_v[l][0] == vi:
                        other = endpoints_v[l][1]
                    else:
                        other = endpoints_v[l][0]
                    if other not in visited_vertices:
                        vertices_to_visit.append(other)
                    visited_edges.add(l)
                visited_vertices.add(vi)
        return (len(visited_vertices) == len(self.vertices))


    def is_loop(self, edge):
        """Return `True` if `edge` is a loop (i.e., the two endpoint coincide).
        """
        return self.endpoints_v[edge][0] == self.endpoints_v[edge][1]
        

    def is_orientation_reversing(self, automorphism):
        """Return `True` if `automorphism` reverses orientation of
        this `Fatgraph` instance."""
        return (-1 == Fatgraph._cmp_orient(self, self, automorphism))


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
            if self.is_orientation_reversing(a):
                return False
        # no orientation reversing automorphism found
        return True


    @cache_iterator
    def isomorphisms(g1, g2):
        """Iterate over `Fatgraph` isomorphisms from `g1` to `g2`.

        An isomorphism is represented by a tuple `(pv, rot, pe)` where:

          - `pv` is a permutation of ther vertices: the `i`-th vertex
            of `g1` is sent to the `pv[i]`-th vertex of `g2`, rotated
            by `rot[i]` places leftwards;

          - `pe` is a permutation of the edges: edge `i` in `g1` is
            mapped to edge `pe[i]` in `g2`.

        This method can iterate over the automorphism group of a
        graph::

          >>> g1 = Fatgraph([Vertex([2, 1, 1]), Vertex([2, 0, 0])])
          >>> for f in g1.isomorphisms(g1): print f
          ({0: 0, 1: 1}, [0, 0], {0: 0, 1: 1, 2: 2})
          ({0: 1, 1: 0}, [0, 0], {0: 1, 1: 0, 2: 2})

        Or it can find the isomorphisms between two given graphs::

          >>> g2 = Fatgraph([Vertex([2, 2, 0]), Vertex([1, 1, 0])])
          >>> for f in g1.isomorphisms(g2): print f
          ({0: 0, 1: 1}, [2, 2], {0: 1, 1: 2, 2: 0})
          ({0: 1, 1: 0}, [2, 2], {0: 2, 1: 1, 2: 0})

        If there are no isomorphisms connecting the two graphs, then no
        item is returned by the iterator::

          >>> g3 = Fatgraph([Vertex([2, 1, 0]), Vertex([2, 0, 1])])
          >>> list(g1.isomorphisms(g3))
          []
        """
        
        # As this procedure is quite complex, we break it into a
        # number of auxiliary functions.
        
        def starting_vertices(graph):
            """
            Return the pair `(valence, vertices)`, which minimizes the
            product of valence with the number of vertices of that
            valence.

            Examples::

              >>> g = Fatgraph([Vertex([0,1,2]), Vertex([0,2,3]), Vertex([3,4,4,5,5])])
              >>> g.starting_vertices()
              (5, [Vertex([3, 4, 4, 5, 5])])
            """
            val = max(graph.vertex_valences())
            vs = None
            n = len(graph.vertices)+1
            for (val_, vs_) in graph.valence_spectrum().iteritems():
                n_ = len(vs_)
                if (n_*val_ < n*val) \
                       or (n_*val_ == n*val and val_<val):
                    val = val_
                    vs = vs_
                    n = n_
            return (val, vs)

        def compatible(v1, v2):
            """Return `True` if vertices `v1` and `v2` are compatible.
            (i.e., same valence and number of loops - one *could* be
            mapped onto the other.)
            """
            if len(v1) == len(v2) and v1.num_loops() == v2.num_loops():
                return True
            else:
                return False
                
        def admissible_vertex_mappings(v, g, ixs):
            """Iterate over all (indices of) vertices in `g`, which
            `v` *could* be mapped to (that is, the destination vertex
            matches `v` in valence and number of loops.

            Third argument `ixs` restricts the search to the given
            subset of vertex indices in `g`.
            """
            for i in ixs:
                if compatible(v, g.vertices[i]):
                    yield i

        class CannotExtendMap(Exception):
            """Exception raised by `extend_map` on failure to extend a
            partial map.
            """
            pass

        def extend_map(m, g1, i1, r, g2, i2):
            """Extend map `m` by mapping the `i1`-th vertex in `g1` to
            the `i2`-th vertex in `g2` (and rotating the source vertex
            by `r` places leftwards).  Return the extended map `(pv,
            rot, pe)`.

            The partial map `m` is represented as a triple `(pv, rot,
            pe)` as in `Fatgraph.isomorphism` (which see), with the
            additional proviso that unassigned items in `rot` are
            represented by `None`.
            """
            (pv, rots, pe) = m
            v1 = g1.vertices[i1]
            v2 = g2.vertices[i2]
            if not compatible(v1, v2):
                raise CannotExtendMap

            # XXX: rotation has to be >=0 for the [r:r+..] shift below to work
            if r < 0:
                r += len(v2)

            if pv.has_key(i1):
                if pv[i1] != i2 or (rots[i1] - r) % len(v2) != 0:
                    raise CannotExtendMap
                else:
                    # this pair has already been added
                    return m

            pv[i1] = i2
            rots[i1] = r

            # rotating `v1` leftwards is equivalent to rotating `v2` rightwards...
            v2 = v2[r:r+len(v2)]
            if not pe.extend(v1, v2):
                raise CannotExtendMap

            return (pv, rots, pe)

        def neighbors(m, g1, i1, g2, i2):
            """List of vertex-to-vertex mappings that extend map
            `m` in the neighborhood of of `i1` (in the domain) and
            `i2` (in the codomain).

            Return a list of triplets `(src, dst, rot)`, where:
               * `src` is the index of a vertex in `g1`,
                 connected to `i1` by an edge `x`;
               * `dst` is the index of a vertex in `g2`,
                 connected to `i2` by the image (according to `m`)
                 of edge `x`;
               * `rot` is the rotation to be applied to `g1[src]`
                 so that edge `x` and its image appear
                 at the same index position;
            
            XXX: prune vertices that are already mapped by `m`?
            """
            result = []
            for x in g1[i1]:
                    src_endpoints_v = g1.endpoints_v[x]
                    # ignore loops
                    if src_endpoints_v[0] == src_endpoints_v[1]:
                        continue # to next `x`
                    src_endpoints_i = g1.endpoints_i[x]
                    src_v = src_endpoints_v[0] if (src_endpoints_v[1] == i1) else src_endpoints_v[1]
                    # ignore vertices that are already in the domain of `m`
                    if src_v in m[0]:
                        continue # to next `x`
                    src_i = src_endpoints_i[0] if (src_endpoints_v[1] == i1) else src_endpoints_i[1]
                    dst_endpoints_v = g2.endpoints_v[m[2][x]]
                    dst_endpoints_i = g2.endpoints_i[m[2][x]]
                    dst_v = dst_endpoints_v[0] if (dst_endpoints_v[1] == i2) else dst_endpoints_v[1]
                    dst_i = dst_endpoints_i[0] if (dst_endpoints_v[1] == i2) else dst_endpoints_i[1]
                    # array of (source vertex index, dest vertex index, rotation)
                    result.append((src_v, dst_v, dst_i-src_i))
            return result
            
        # if graphs differ in vertex valences, no isomorphisms
        vs1 = g1.valence_spectrum()
        vs2 = g2.valence_spectrum()
        if not set(vs1.keys()) == set(vs2.keys()):
            return # StopIteration
        # if graphs have unequal vertex distribution by valence, no isomorphisms
        for val in g1.vertex_valences():
            if len(vs1[val]) != len(vs2[val]):
                return # StopIteration

        (val, vs) = starting_vertices(g1)
        src0 = vs[0]
        v1 = g1[src0]
        for dst0 in admissible_vertex_mappings(v1, g2, vs2[val]):
            for rot0 in xrange(val):
                try:
                    # pass 0: init new (pv, rot, pe) triple
                    pv = Permutation()
                    rot = [ None for x in xrange(g1.num_vertices) ]
                    pe = Permutation()

                    # pass 1: map `v1` to `v2` and build map
                    # of neighboring vertices for next pass
                    pv[src0] = dst0
                    rot[src0] = rot0
                    if not pe.extend(v1, g2[dst0][rot0:rot0+val]):
                        continue # to next `rot0`
                    if __debug__:
                        for x in v1:
                            assert x in pe, "Edge `%d` of vertex `%s` (in graph `%s`) not mapped to any edge of graph `%s` (at line 1740, `pe=%s`)" % (x, v1, g1, g2, pe)

                    # pass 2: extend map to neighboring vertices
                    m = (pv, rot, pe)
                    nexts = neighbors(m, g1, src0, g2, dst0)
                    while len(m[0]) < g1.num_vertices:
                        neighborhood = []
                        for (i1, i2, r) in nexts:
                            m = extend_map(m, g1, i1, r, g2, i2)
                            if __debug__:
                                for x in g1[i1]:
                                    assert x in m[2], "Edge `%d` of vertex `%s` (in graph `%s`) not mapped to any edge of graph `%s` (at line 1751, `pe=%s`)" % (x, g1[i1], g1, g2, m[2])
                            neighborhood += neighbors(m, g1, i1, g2, i2)
                        nexts = neighborhood

                # extension failed in the above block, continue with next candidate
                except CannotExtendMap:
                    continue # to next `rot0`

                # sanity checks
                assert len(m[0]) == g1.num_vertices
                assert len(m[2]) == g1.num_edges
                assert set(m[0].keys()) == set(range(g1.num_vertices))
                assert set(m[0].values()) == set(range(g2.num_vertices))
                assert set(m[2].keys()) == set(range(g1.num_edges))
                assert set(m[2].values()) == set(range(g2.num_edges))

                # finally
                yield m


    def num_boundary_cycles(self):
        """Return the number of boundary components of this `Fatgraph` object.

        Each boundary component is represented by the list of (colored)
        edges.

        Examples::
          >>> Fatgraph([Vertex([2,1,0]), Vertex([2,1,0])]).num_boundary_cycles()
          1
          >>> Fatgraph([Vertex([2,1,0]), Vertex([2,0,1])]).num_boundary_cycles()
          3
        """
        return len(self.boundary_cycles())


    ##@cache_symmetric -- not important, as comparisons are ever done one-way
    def projection(self, other):
        """Return the component of the projection of `self` on the
        basis vector `other`.  This can be either 0 (if `self` and
        `other` are not isomorphic), or +1/-1 depending on comparison
        of the orientation of `self` with the pull-back orientation on
        `other`.

        If the two graphs are *not* isomorphic, then the result is 0::

          >>> g1 = Fatgraph([Vertex([0,1,2,0,1,2])])
          >>> g2 = Fatgraph([Vertex([0,1,2]), Vertex([0,1,2])])
          >>> Fatgraph.projection(g1, g2)
          0

        Any graph obviously projects onto itself with coefficient `1`::

          >>> Fatgraph.projection(g1, g1)
          1
          >>> Fatgraph.projection(g2, g2)
          1

        And similarly does any graph isomorphic to a given graph::
        
          >>> g3 = Fatgraph(g2)
          >>> Fatgraph.projection(g2, g3)
          1

        Flipping the orientation on an edge reverses the coefficient
        sign::

          >>> g4 = Fatgraph(g2)
          >>> # make a copy of `g4.edge_numbering`, since it's shared with `g2`
          >>> g4 .edge_numbering = copy(g4.edge_numbering)
          >>> g4.edge_numbering[0], g4.edge_numbering[1] = \
              g4.edge_numbering[1], g4.edge_numbering[0]
          >>> Fatgraph.projection(g2, g4)
          -1
          
        """
        assert isinstance(other, Fatgraph), \
               "Fatgraph.projection:" \
               " called with non-Fatgraph argument `other`: %s" % other
        assert self.is_oriented(), \
               "Fatgraph.projection: cannot project non-orientable graph: %s" \
               % self
        assert other.is_oriented(), \
               "Fatgraph.projection: cannot project non-orientable graph: %s" \
               % other
        try:
            iso = self.isomorphisms(other).next()
            return Fatgraph._cmp_orient(self, other, iso)
        except StopIteration:
            # list of morphisms is empty, graphs are not equal.
            return 0

    @cache1
    def valence_spectrum(self):
        """Return a dictionary mapping valences into vertex indices.

        Examples::

           >>> Fatgraph([Vertex([1,1,0,0])]).valence_spectrum()
           {4: [0]}

           >>> Fatgraph([Vertex([1,1,0]), Vertex([2,2,0])]).valence_spectrum()
           {3: [0, 1]}

           >>> Fatgraph([Vertex([3, 1, 0, 1]), \
                      Vertex([4, 4, 0]), Vertex([3, 2, 2])]).valence_spectrum()
           {3: [1, 2], 4: [0]}
        """
        result = {}
        for (index, vertex) in enumerate(self.vertices):
            l = len(vertex)
            if l in result:
                result[l].append(index)
            else:
                result[l] = [index]
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

    @cache1
    def vertex_valences(self):
        return frozenset(len(v) for v in self.vertices)

    @cache1
    def vertex_valence_distribution(self):
        spec = self.valence_spectrum()
        return dict((v, len(spec[v]))
                    for v in spec.iterkeys())
    
    
def MakeNumberedGraphs(graph):
    """Return all distinct (up to isomorphism) decorations of
    `Fatgraph` instance `graph` with a numbering of the boundary
    cycles.

    Examples::
    
      >>> ug1 = Fatgraph([Vertex([2,0,0]), Vertex([2,1,1])])
      >>> for g in MakeNumberedGraphs(ug1): print g
      NumberedFatgraph(Fatgraph([Vertex([2, 0, 0]), Vertex([2, 1, 1])]),
                       numbering={BoundaryCycle([(0, 2, 0), (1, 2, 0), (0, 0, 1), (1, 0, 1)]): 0,
                                  BoundaryCycle([(1, 1, 2)]): 1,
                                  BoundaryCycle([(0, 1, 2)]): 2})
      NumberedFatgraph(Fatgraph([Vertex([2, 0, 0]), Vertex([2, 1, 1])]),
                       numbering={BoundaryCycle([(0, 1, 2)]): 0,
                                  BoundaryCycle([(0, 2, 0), (1, 2, 0), (0, 0, 1), (1, 0, 1)]): 1,
                                  BoundaryCycle([(1, 1, 2)]): 2})
      NumberedFatgraph(Fatgraph([Vertex([2, 0, 0]), Vertex([2, 1, 1])]),
                       numbering={BoundaryCycle([(0, 1, 2)]): 0,
                                  BoundaryCycle([(1, 1, 2)]): 1,
                                  BoundaryCycle([(0, 2, 0), (1, 2, 0), (0, 0, 1), (1, 0, 1)]): 2})
       
    Note that, when only one numbering out of many possible ones is
    returned because of isomorphism, the returned numbering may not be
    the trivial one (it is infact the first permutation of 0..n
    returned by `InplacePermutationIterator`)::
      
      >>> ug2 = Fatgraph([Vertex([2,1,0]), Vertex([2,0,1])])
      >>> MakeNumberedGraphs(ug2)
      [NumberedFatgraph(Fatgraph([Vertex([2, 1, 0]), Vertex([2, 0, 1])]),
                        numbering={BoundaryCycle([(1, 2, 0), (0, 0, 1)]): 0,
                                   BoundaryCycle([(0, 2, 0), (1, 0, 1)]): 1,
                                   BoundaryCycle([(0, 1, 2), (1, 1, 2)]): 2})]

    When the graph has only one boundary component, there is only one
    possible numbering, which is actually returned::
    
      >>> ug3 = Fatgraph([Vertex([1,0,1,0])])
      >>> MakeNumberedGraphs(ug3)
      [NumberedFatgraph(Fatgraph([Vertex([1, 0, 1, 0])]),
                        numbering={BoundaryCycle([(0, 3, 0), (0, 2, 3), (0, 1, 2), (0, 0, 1)]): 0})]
      
    """
    bc = graph.boundary_cycles()
    n = len(bc) # == graph.num_boundary_cycles()

    ## Find out which automorphisms permute the boundary cycles among
    ## themselves.
    # XXX: implement `Permutation.__hash__()` and turn `K` into a `set`.
    K = []
    for a in graph.automorphisms():
        k = Permutation()
        k_is_ok = True
        for src in xrange(n):
            dst_cy = bc[src].transform(a)
            try:
                dst = bc.index(dst_cy)
            except ValueError:
                # `dst_cy` not in `bc`
                k_is_ok = False
                break # continue with next `a`
            k[src] = dst
        if k_is_ok and (k not in K):
            # `a` induces permutation `k` on the set `bc`
            K.append(k)

    ## There will be as many distinct numberings as there are cosets
    ## of `K` in `Sym(n)`.
    if len(K) > 1:
        def unseen(candidate, K, already):
            """Return `False` iff any of the images of `candidate` by an
            element of group `K` is contained in set `already`.
            """
            for k in K:
                if k.rearrange(candidate) in already:
                    return False
            return True
        numberings = [ ]
        for candidate in InplacePermutationIterator(range(n)):
            if unseen(candidate, K, numberings):
                numberings.append(copy(candidate))
    else:
        # if `K` is the one-element group, then all orbits are trivial
        numberings = PermutationIterator(range(n))
    
    result = [ NumberedFatgraph(graph, zip(bc, numbering))
               for numbering in numberings ]
    return result



class NumberedFatgraph(Fatgraph):
    """A `Fatgraph` decorated with a numbering of the boundary components.

    A numbered fatgraph is constructed from a `Fatgraph` instance
    (called the *underlying graph*) and a numbering (that is, a
    bijective map assigning an integer to each boundary components of
    the underlying graph).

    Examples::

      >>> BoundaryCycle = Fatgraph.BoundaryCycle # shortcut
      >>> ug = Fatgraph([Vertex([1, 0, 1, 0])])  # underlying graph
      >>> ng = NumberedFatgraph(ug, \
                 numbering=[(BoundaryCycle([(0,3,0), (0,2,3), (0,1,2), (0,0,1)]), 0)])

    Since `NumberedFatgraphs` are just decorated `Fatgraphs`, they
    only differ in the way two `NumberedFatgraph` instances are deemed
    isomorphic:
      - they must have isomorphic underlying graphs;
      - the numberings must match under the isomorphism map.
    For example::
    
      >>> ug = Fatgraph([Vertex([1, 1, 0, 0])])
      >>> ng1 = NumberedFatgraph(ug, numbering=[(BoundaryCycle([(0,3,0), (0,1,2)]), 0), \
                                                (BoundaryCycle([(0,0,1)]), 1),          \
                                                (BoundaryCycle([(0,2,3)]), 2)])
      >>> ng2 = NumberedFatgraph(ug, numbering=[(BoundaryCycle([(0,3,0), (0,1,2)]), 1), \
                                                (BoundaryCycle([(0,0,1)]), 0),          \
                                                (BoundaryCycle([(0,2,3)]), 2)])
      >>> ng1 == ng2
      False

      Fatgraph instances equipped with a numbering are compared as
      numbered graphs (that is, the isomorphism should transform the
      numbering on the source graph onto the numbering of the
      destination)::

        >>> NumberedFatgraph.__eq__(
        ...     NumberedFatgraph(Fatgraph([Vertex([2,0,1]), Vertex([2,1,0])]), 
        ...                      numbering=[(BoundaryCycle([(0,2,0), (1,0,1)]), 0), 
        ...                                 (BoundaryCycle([(0,0,1), (1,2,0)]), 1), 
        ...                                 (BoundaryCycle([(0,1,2), (1,1,2)]), 2) ] ),
        ...     NumberedFatgraph(Fatgraph([Vertex([2,0,1]), Vertex([2,1,0])]), 
        ...                      numbering=[(BoundaryCycle([(0,2,0), (1,0,1)]), 0), 
        ...                                 (BoundaryCycle([(0,0,1), (1,2,0)]), 2), 
        ...                                 (BoundaryCycle([(0,1,2), (1,1,2)]), 1) ] ) )
        True

        >>> NumberedFatgraph.__eq__(
        ...     NumberedFatgraph(Fatgraph([Vertex([1, 0, 0, 2, 2, 1])]), 
        ...                      numbering=[(BoundaryCycle([(0,5,0)]), 0), 
        ...                                 (BoundaryCycle([(0,0,1), (0,2,3), (0,4,5)]), 1), 
        ...                                 (BoundaryCycle([(0,1,2)]), 3), 
        ...                                 (BoundaryCycle([(0,3,4)]), 2) ]),
        ...     NumberedFatgraph( 
        ...                      Fatgraph([Vertex([2, 2, 1, 1, 0, 0])]), 
        ...                      numbering=[(BoundaryCycle([(0,2,3)]), 0), 
        ...                                 (BoundaryCycle([(0,3,4), (0,5,0), (0,1,2)]), 3), 
        ...                                 (BoundaryCycle([(0,4,5)]), 1), 
        ...                                 (BoundaryCycle([(0,0,1)]), 2) ]) )
        False

        >>> NumberedFatgraph.__eq__(
        ...     NumberedFatgraph(Fatgraph([Vertex([1, 0, 0, 2, 2, 1])]), 
        ...                      numbering=[(BoundaryCycle([(0,5,0)]), 0), 
        ...                                 (BoundaryCycle([(0,0,1), (0,2,3), (0,4,5)]), 1),
        ...                                 (BoundaryCycle([(0,1,2)]), 3), 
        ...                                 (BoundaryCycle([(0,3,4)]), 2) ]),
        ...     NumberedFatgraph(Fatgraph([Vertex([2, 2, 1, 1, 0, 0])]), 
        ...                      numbering=[(BoundaryCycle([(0,2,3)]), 3), 
        ...                                 (BoundaryCycle([(0,3,4), (0,5,0), (0,1,2)]), 2), 
        ...                                 (BoundaryCycle([(0,4,5)]), 0), 
        ...                                 (BoundaryCycle([(0,0,1)]), 1) ]) )
        False

        >>> NumberedFatgraph.__eq__(
        ...     NumberedFatgraph(Fatgraph([Vertex([3, 2, 2, 0, 1]), Vertex([3, 1, 0])]), 
        ...                      numbering=[(BoundaryCycle([(0,4,0), (1,0,1)]), 0),
        ...                                 (BoundaryCycle([(0,0,1), (0,2,3), (1,2,0)]), 1),
        ...                                 (BoundaryCycle([(0,1,2)]), 2),
        ...                                 (BoundaryCycle([(0,3,4), (1,1,2)]), 3) ]),
        ...     NumberedFatgraph(Fatgraph([Vertex([2, 3, 1]), Vertex([2, 1, 3, 0, 0])]), 
        ...                      numbering=[(BoundaryCycle([(0,2,0), (1,0,1)]), 3),
        ...                                 (BoundaryCycle([(0,0,1), (1,2,3), (1,4,0)]), 1),
        ...                                 (BoundaryCycle([(0,1,2), (1,1,2)]), 0) ,
        ...                                 (BoundaryCycle([(1,3,4)]), 2) ]) )
        True

        >>> NumberedFatgraph.__eq__(
        ...     NumberedFatgraph(Fatgraph([Vertex([0, 1, 2, 0, 2, 1])]),
        ...                      numbering=[(BoundaryCycle([(0,1,2), (0,4,5)]), 0),
        ...                                 (BoundaryCycle([(0,5,0), (0,3,4), (0,2,3), (0,0,1)]), 1)]),
        ...     NumberedFatgraph(Fatgraph([Vertex([0, 1, 2, 0, 2, 1])]),
        ...                      numbering=[(BoundaryCycle([(0,1,2), (0,4,5)]), 1),
        ...                                 (BoundaryCycle([(0,5,0), (0,3,4), (0,2,3), (0,0,1)]), 0)]) )
        False
      """

    def __init__(self, underlying, numbering, vertextype=Vertex):
        Fatgraph.__init__(self, underlying)
        self.underlying = underlying
        self.numbering = numbering


    def __repr__(self):
        """Output a printed representation, such that `eval(repr(x)) == x`.

        The `numbering` attribute of the `NumberedFatgraph` differs from
        the standard Python printing of dictionaries, in that its printed
        form is sorted by values (i.e., by boundary component index)
        to make doctests more stable.
        """
        def sort_by_values_then_by_keys(kv1,kv2):
            (k1, v1) = kv1
            (k2, v2) = kv2
            q = cmp(v1, v2)
            if 0 == q:
                return cmp(k1, k2)
            else:
                return q
        parts = []
        for (k, v) in sorted(self.numbering.iteritems(),
                             cmp=sort_by_values_then_by_keys):
            parts.append("%s: %s" % (repr(k), repr(v)))
        canonical = "{" + str.join(", ", parts) + "}"
        return "NumberedFatgraph(%s, numbering=%s)" \
               % (repr(self.underlying), canonical)
    

    @cache
    def contract(self, edgeno):
        """Return a new `NumberedFatgraph` instance, obtained by
        contracting the specified edge.

        Examples::

          >>> g0 = NumberedFatgraph(Fatgraph([Vertex([1, 2, 1]), Vertex([2, 0, 0])]),
          ...                       numbering={BoundaryCycle([(0, 1, 2), (1, 2, 0),
          ...                                                 (0, 0, 1), (1, 0, 1)]): 0,
          ...                                  BoundaryCycle([(1, 1, 2)]): 1,
          ...                                  BoundaryCycle([(0, 2, 0)]): 2})
          >>> g0.contract(2)
          NumberedFatgraph(Fatgraph([Vertex([1, 1, 0, 0])]),
                           numbering={BoundaryCycle([(0, 3, 0), (0, 1, 2)]): 0,
                                      BoundaryCycle([(0, 2, 3)]): 1,
                                      BoundaryCycle([(0, 0, 1)]): 2})

          >>> g1 = NumberedFatgraph(Fatgraph([Vertex([1, 0, 2]), Vertex([2, 1, 5]),
          ...                                 Vertex([0, 4, 3]), Vertex([4, 5, 3])]),
          ...                       numbering={BoundaryCycle([(0, 2, 0), (0, 1, 2), (0, 0, 1),
          ...                                                 (1, 1, 2), (2, 0, 1), (3, 0, 1),
          ...                                                 (1, 2, 0), (2, 2, 0), (1, 0, 1),
          ...                                                                       (3, 1, 2)]): 0,
          ...                                  BoundaryCycle([(3, 2, 0), (2, 1, 2)]): 1})
          >>> g1.contract(0)
          NumberedFatgraph(Fatgraph([Vertex([1, 0, 3, 2]), Vertex([1, 0, 4]), Vertex([3, 4, 2])]),
                           numbering={BoundaryCycle([(2, 1, 2), (0, 1, 2), (0, 0, 1), (1, 1, 2),
                                                     (0, 3, 0), (2, 0, 1), (1, 2, 0), (1, 0, 1)]): 0,
                                      BoundaryCycle([(0, 2, 3), (2, 2, 0)]): 1})

        """
        # check that we are not contracting a loop or an external edge
        assert not self.is_loop(edgeno), \
               "NumberedFatgraph.contract: cannot contract a loop."
        assert (edgeno >= 0) and (edgeno < self.num_edges), \
               "NumberedFatgraph.contract: invalid edge number (%d):"\
               " must be in range 0..%d" \
               % (edgeno, self.num_edges)

        # store endpoints of the edge-to-be-contracted
        (v1, v2) = self.endpoints_v[edgeno]
        (pos1, pos2) = self.endpoints_i[edgeno]
        if v1 > v2:
            # swap endpoints so that `v1 < v2`
            v1, v2 = v2, v1
            pos1, pos2 = pos2, pos1
        l1 = len(self.vertices[v1])
        l2 = len(self.vertices[v2])
        
        # transform corners according to contraction; see
        # `Fatgraph.contract()` for an explanation of how the
        # underlying graph is altered during contraction.
        new_numbering = {}
        for (bcy, n) in self.numbering.iteritems():
            new_cy = []
            for corner in bcy:
                if corner[0] == v1:
                    if pos1 == corner[1]:
                        # skip this corner, keep only one of the
                        # corners limited by the contracted edge
                        continue
                    else:
                        i1 = (corner[1] - pos1 - 1) % l1
                        i2 = (corner[2] - pos1 - 1) % l1
                        new_cy.append((v1, i1, i2))
                elif corner[0] == v2:
                    if pos2 == corner[1]:
                        # skip this corner, keep only one of the
                        # corners limited by the contracted edge
                        continue
                    if pos2 == corner[2]:
                        new_cy.append((v1, l1+l2-3, 0))
                    else:
                        i1 = l1-1 + ((corner[1] - pos2 - 1) % l2)
                        i2 = l1-1 + ((corner[2] - pos2 - 1) % l2)
                        new_cy.append((v1, i1, i2))
                elif corner[0] > v2:
                    # shift vertices after `v2` one position down
                    new_cy.append((corner[0]-1, corner[1], corner[2]))
                else:
                    # pass corner unchanged
                    new_cy.append(corner)
            new_numbering[Fatgraph.BoundaryCycle(new_cy)] = n

        return NumberedFatgraph(self.underlying.contract(edgeno),
                                numbering=new_numbering)
        

    @cache_iterator
    def isomorphisms(self, other):
        """Iterate over isomorphisms from `self` to `other`.

        See `Fatgraph.isomrphisms` for a discussion of the
        representation of isomorphisms and example usage.
        """
        for (pv, rot, pe) in Fatgraph.isomorphisms(self.underlying, other.underlying):
            pe_does_not_preserve_bc = False
            for bc1 in self.underlying.boundary_cycles():
                bc2 = bc1.transform((pv, rot, pe))
                # there are cases (see examples in the
                # `Fatgraph.__eq__` docstring, in which the
                # above algorithm may find a valid
                # mapping, changing from `g1` to an
                # *alternate* representation of `g2` -
                # these should fail as they don't preserve
                # the boundary cycles, so we catch them
                # here.
                if (bc2 not in other.numbering) \
                       or (self.numbering[bc1] != other.numbering[bc2]):
                    pe_does_not_preserve_bc = True
                    break
            if pe_does_not_preserve_bc:
                continue # to next underlying graph isomorphism
            yield (pv, rot, pe)


    def _numbering_get(self):
        """Return the numbering previously set on this instance via
        `._numbering_set()`.
        """
        return self._numbering

    def _numbering_set(self, initializer):
        """Set the `.numbering` attribute from a valid `dict`
        initializer.

        The initializer can be:
          - either a sequence of tuples `(bcy, n)`, where each `n` is
            a non-negative integer, and each `bcy` is a
            `BoundaryCycle` instance,
          - or a `dict` instance mapping `BoundaryCycle` instances to
            `int`s.
        In either case, an assertion is raised if:
          - the number of pairs in the initializer does not match
            the number of boundary cycles;
          - the set of integer keys is not `[0 .. n]`;
          - there are duplicate boundary cycles or integers
            in the initializer;

        The `numbering` attribute is set to a dictionary mapping the
        boundary cycle `bcy` to the integer `n`::

          >>> ug0 = Fatgraph([Vertex([1,2,0]), Vertex([1,0,2])])
          >>> bc = ug0.boundary_cycles()  # three b.c.'s
          >>> ng0 = NumberedFatgraph(ug0, [ (bcy,n) for (n,bcy) in enumerate(bc)])
          >>> ng0.numbering == {
          ...    Fatgraph.BoundaryCycle([(0,0,1), (1,2,0)]): 0, 
          ...    Fatgraph.BoundaryCycle([(0,1,2), (1,1,2)]): 1, 
          ...    Fatgraph.BoundaryCycle([(0,2,0), (1,0,1)]): 2,
          ... }
          True
        """
        assert len(initializer) == self.num_boundary_cycles()
        self._numbering = dict(initializer)
        if __debug__:
            count = [ 0 for x in xrange(self.num_boundary_cycles()) ]
            for (bcy,n) in self._numbering.iteritems():
                assert type(n) is types.IntType, \
                       "NumberedFatgraph._numbering_set: 2nd argument has wrong type:" \
                       " expecting (BoundaryCycle, Int) pair, got `(%s, %s)`." \
                       " Reversed-order arguments?" \
                       % (bcy, n)
                assert isinstance(bcy, Fatgraph.BoundaryCycle), \
                       "NumberedFatgraph._numbering_set: 1st argument has wrong type:" \
                       " expecting (BoundaryCycle, Int) pair, got `(%s, %s)`." \
                       " Reversed-order arguments?" \
                       % (bcy, n)
                assert bcy in self.boundary_cycles(), \
                       "NumberedFatgraph._numbering_set():" \
                       " Cycle `%s` is no boundary cycle of graph `%s` " \
                       % (bcy, self.underlying)
                count[n] += 1
                if count[n] > 1:
                    raise AssertionError("NumberedFatgraph._numbering_set():" \
                                         " Duplicate key %d" % n)
            assert sum(count) != self.num_boundary_cycles()-1, \
                   "NumberedFatgraph._numbering_set():" \
                   " Initializer does not exhaust range `0..%d`: %s" \
                   % (self.num_boundary_cycles()-1, initializer)


    numbering = property(_numbering_get, _numbering_set)
    
    

@cache_iterator
def MgnTrivalentGraphsRecursiveGenerator(g, n):
    """Iterate over all connected trivalent fatgraphs having the
    prescribed genus `g` and number of boundary cycles `n`.
    
    Examples::

      >>> for g in MgnTrivalentGraphsRecursiveGenerator(0,3): print g
      Fatgraph([Vertex([1, 2, 1]), Vertex([2, 0, 0])]) 
      Fatgraph([Vertex([1, 0, 2]), Vertex([2, 0, 1])])

      >>> for g in MgnTrivalentGraphsRecursiveGenerator(1,1): print g
      Fatgraph([Vertex([1, 0, 2]), Vertex([2, 1, 0])])

    """
    # avoid infinite recursion in later statements
    if n==0 or (g,n)<(0,3):
        raise StopIteration

    # sanity check
    assert n > 0, \
           "MgnTrivalentGraphsRecursiveGenerator: " \
           " number of boundary cycles `n` must be positive,"\
           " but got `%s` instead" % n

    logging.debug("Starting MgnTrivalentGraphsRecursiveGenerator(%d,%d) ..." % (g,n))

    ## M_{0,3} - induction base
    if (g,n) == (0,3):
        yield Fatgraph([Vertex([1, 2, 1]), Vertex([2, 0, 0])])
        yield Fatgraph([Vertex([1, 0, 2]), Vertex([2, 0, 1])])
        logging.debug("  MgnTrivalentGraphsRecursiveGenerator(0,3) done.")

    ## M_{1,1} - induction base
    elif (g,n) == (1,1):
        yield Fatgraph([Vertex([1, 0, 2]), Vertex([2, 1, 0])])
        logging.debug("  MgnTrivalentGraphsRecursiveGenerator(1,1) done.")

    ## General case
    else:
        def graphs(g,n):
            logging.debug("  MgnTrivalentGraphsRecursiveGenerator(%d,%d): "
                          "pass 1: hang a circle to all edges of graphs in M_{%d,%d} ..." % (g,n, g,n-1))
            for G in MgnTrivalentGraphsRecursiveGenerator(g,n-1):
                for x in G.edge_orbits():
                    yield G.hangcircle(x,0)
                    yield G.hangcircle(x,1)

            logging.debug("  MgnTrivalentGraphsRecursiveGenerator(%d,%d): "
                          "pass 2: bridge all edges of a single graph in M_{%d,%d} ..." % (g,n, g,n-1))
            for G in MgnTrivalentGraphsRecursiveGenerator(g,n-1):
                for (x,y) in G.edge_pair_orbits():
                        yield G.bridge(x,0, y,0)
                        yield G.bridge(x,0, y,1)
                        yield G.bridge(x,1, y,0)
                        yield G.bridge(x,1, y,1)

            logging.debug("  MgnTrivalentGraphsRecursiveGenerator(%d,%d): "
                          "pass 3: bridge all edges of a single graph in M_{%d,%d} ..." % (g,n, g-1,n+1))
            for G in MgnTrivalentGraphsRecursiveGenerator(g-1,n+1):
                for (x,y) in G.edge_pair_orbits():
                        yield G.bridge(x,0, y,0)
                        yield G.bridge(x,0, y,1)
                        yield G.bridge(x,1, y,0)
                        yield G.bridge(x,1, y,1)

            logging.debug("  MgnTrivalentGraphsRecursiveGenerator(%d,%d): "
                          "pass 4: bridge two graphs of such that g_1+g_2=%d, n_1+n_2=%d ..." % (g,n, g,n+1)) 
            def add_up_to(x, min=0):
                if x == 0 and min == 0:
                    yield (0,0)
                elif x-min >= 0:
                    for y in xrange(min, x-min+1):
                        yield (y, x-y)
            for (g1, g2) in add_up_to(g, min=0):
                for (n1, n2) in add_up_to(n+1, min=1):
                    if (g1, n1) < (0, 3) or (g2, n2) < (0,3):
                        continue
                    for G1 in MgnTrivalentGraphsRecursiveGenerator(g1,n1):
                        for G2 in MgnTrivalentGraphsRecursiveGenerator(g2,n2):
                            for x1 in G1.edge_orbits():
                                for x2 in G2.edge_orbits():
                                    yield Fatgraph.bridge2(G1, x1, 0, G2, x2, 0)
                                    yield Fatgraph.bridge2(G1, x1, 0, G2, x2, 1)
                                    yield Fatgraph.bridge2(G1, x1, 1, G2, x2, 0)
                                    yield Fatgraph.bridge2(G1, x1, 1, G2, x2, 1)

        unique = []
        for G in graphs(g,n):
            # XXX: should this check be done in graphs(g,n)?
            if (G.genus(), G.num_boundary_cycles()) != (g,n) \
                   or (G in unique):
                G.release()
                continue
            unique.append(G)
            yield G

        logging.debug("  MgnTrivalentGraphsRecursiveGenerator(%d,%d) done, found %d unique graphs." % (g,n, len(unique)))



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

    def __init__(self, g, n, trivalent_graphs_generator=MgnTrivalentGraphsRecursiveGenerator):
        assert n > 0, \
               "MgnGraphsIterator: " \
               " number of boundary cycles `n` must be positive,"\
               " but got `%s` instead" % n
        assert (g > 0) or (g == 0 and n >= 3), \
               "MgnGraphsIterator: " \
               " Invalid (g,n) pair (%d,%d): "\
               " need either g>0 or g==0 and n>2" \
               % (g,n)

        #: Prescribed genus of returned graphs
        self.g = g

        #: Prescribed number of boundary components
        self.n = n
        
        # Gather all 3-valent graphs.
        trivalent = list(trivalent_graphs_generator(g,n))

        #: Fatgraphs to be contracted at next `.refill()` invocation
        self._batch = trivalent

        #: Number of edges of graphs that will be returned by next
        #  `.refill()` call.  Starts with `6*g + 3*n - 7`, which is the
        #  highest-numbered edge in trivalent graphs.
        self._current_edge = 6*g + 3*n - 7

        self._num_vertices = 4*g + 2*n - 4
        
        # initialize superclass with list of trivalent graphs
        BufferingIterator.__init__(self, trivalent)

    def refill(self):
        if self._num_vertices == 0:
            raise StopIteration

        logging.debug("Generating graphs with %d vertices ...",
                     self._num_vertices)
        discarded = 0
        next_batch = []
        for graph in self._batch:
            # contract all edges
            for edge in xrange(graph.num_edges):
                if not graph.is_loop(edge):
                    dg = graph.contract(edge)
                    if dg not in next_batch:
                        # put graph back into next batch for processing
                        next_batch.append(dg)
                    else:
                        dg.release()
                        discarded += 1
        self._batch = next_batch
        self._num_vertices -= 1

        logging.debug("  Found %d distinct unique graphs with %d vertices, discarded %d.",
                     len(next_batch), 1+self._num_vertices, discarded)
        return next_batch

#MgnGraphsIterator = persist.PersistedIterator(_MgnGraphsIterator)



class MgnNumberedGraphsIterator(BufferingIterator):
    """Iterate over all connected numbered fatgraphs having the
    prescribed genus `g` and number of boundary cycles `n`.
    
    Examples::

      >>> for g in MgnNumberedGraphsIterator(0,3): print g
      NumberedFatgraph(Fatgraph([Vertex([1, 2, 1]), Vertex([2, 0, 0])]),
                       numbering={BoundaryCycle([(0, 1, 2), (1, 2, 0), (0, 0, 1), (1, 0, 1)]): 0,
                                  BoundaryCycle([(1, 1, 2)]): 1,
                                  BoundaryCycle([(0, 2, 0)]): 2})
      NumberedFatgraph(Fatgraph([Vertex([1, 2, 1]), Vertex([2, 0, 0])]),
                       numbering={BoundaryCycle([(0, 2, 0)]): 0,
                                  BoundaryCycle([(0, 1, 2), (1, 2, 0), (0, 0, 1), (1, 0, 1)]): 1,
                                  BoundaryCycle([(1, 1, 2)]): 2})
      NumberedFatgraph(Fatgraph([Vertex([1, 2, 1]), Vertex([2, 0, 0])]),
                       numbering={BoundaryCycle([(0, 2, 0)]): 0,
                                  BoundaryCycle([(1, 1, 2)]): 1,
                                  BoundaryCycle([(0, 1, 2), (1, 2, 0), (0, 0, 1), (1, 0, 1)]): 2})
      NumberedFatgraph(Fatgraph([Vertex([1, 0, 2]), Vertex([2, 0, 1])]),
                       numbering={BoundaryCycle([(0, 0, 1), (1, 1, 2)]): 0,
                                  BoundaryCycle([(0, 2, 0), (1, 2, 0)]): 1,
                                  BoundaryCycle([(0, 1, 2), (1, 0, 1)]): 2})
      NumberedFatgraph(Fatgraph([Vertex([1, 1, 0, 0])]),
                       numbering={BoundaryCycle([(0, 0, 1)]): 0,
                                  BoundaryCycle([(0, 2, 3)]): 1,
                                  BoundaryCycle([(0, 1, 2), (0, 3, 0)]): 2})
      NumberedFatgraph(Fatgraph([Vertex([1, 1, 0, 0])]),
                       numbering={BoundaryCycle([(0, 1, 2), (0, 3, 0)]): 0,
                                  BoundaryCycle([(0, 0, 1)]): 1,
                                  BoundaryCycle([(0, 2, 3)]): 2})
      NumberedFatgraph(Fatgraph([Vertex([1, 1, 0, 0])]),
                       numbering={BoundaryCycle([(0, 2, 3)]): 0,
                                  BoundaryCycle([(0, 1, 2), (0, 3, 0)]): 1,
                                  BoundaryCycle([(0, 0, 1)]): 2})
     

      >>> list(MgnNumberedGraphsIterator(1,1)) == [
      ...     NumberedFatgraph(Fatgraph([Vertex([1, 0, 2]), Vertex([2, 1, 0])]),
      ...                      numbering=[(Fatgraph.BoundaryCycle([(0, 2, 0), (1, 1, 2), (0, 1, 2),
      ...                                                          (1, 0, 1), (0, 0, 1), (1, 2, 0)]), 0) ]),
      ...     NumberedFatgraph(Fatgraph([Vertex([1, 0, 1, 0])]),
      ...                      numbering=[(Fatgraph.BoundaryCycle([(0, 3, 0), (0, 2, 3),
      ...                                                          (0, 1, 2), (0, 0, 1)]), 0) ])
      ...  ]
      True
    """

    def __init__(self, g, n):
        self.__naked_graphs_iterator = MgnGraphsIterator(g, n)
        BufferingIterator.__init__(self)

    def refill(self):
        return MakeNumberedGraphs(self.__naked_graphs_iterator.next())

#MgnNumberedGraphsIterator = persist.PersistedIterator(_MgnNumberedGraphsIterator)



## main: run tests

#import pydb
#pydb.debugger()

if "__main__" == __name__:
    import doctest
    doctest.testmod(name='rg',
                    optionflags=doctest.NORMALIZE_WHITESPACE)
