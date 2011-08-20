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
from itertools import chain,count,izip
import types


## application-local imports

from cache import (
    cache,
    cache_iterator,
    cache_symmetric,
    )
from combinatorics import (
    InplacePermutationIterator,
    SetProductIterator,
    Permutation,
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
##     __slots__ = [
##         'cache',
##         ]
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


class Vertex(CyclicList):
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

##     __slots__ = [ '_num_loops' ]

    def __init__(self, seq=None):
        self._num_loops = None
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
    
    def is_canonical_representative(self):
        """Return `True` if this `Vertex` object is maximal
        (w.r.t. lexicographic order) among representatives of the same
        cyclic sequence.
        
        Examples::
        
          >>> Vertex([3,2,1]).is_canonical_representative()
          True
          >>> Vertex([2,1,3]).is_canonical_representative()
          False
          >>> Vertex([1,1]).is_canonical_representative()
          True
          >>> Vertex([1]).is_canonical_representative()
          True
        """
        L = len(self)
        for i in xrange(1,L):
            for j in xrange(0,L):
                # k := (i+j) mod L
                k = i+j
                if k >= L:
                    k -= L
                if self[k] < self[j]:
                    # continue with next i
                    break
                elif self[k] > self[j]:
                    return False
                # else, continue comparing
        return True

    def make_canonical(self):
        """Alter `Vertex` *in place* so that it is represented by a
        canonical sequence.  Return modified sequence for convenience.
        
        Examples::
          >>> Vertex([3,2,1]).make_canonical()
          Vertex([3, 2, 1])
          >>> Vertex([2,1,3]).make_canonical()
          Vertex([3, 2, 1])
        """
        L = len(self)
        r = 0
        for i in xrange(1,L):
            for j in xrange(0,L):
                # k := (i+j) mod L
                k = i+j
                if k >= L:
                    k -= L
                if self[k] < self[j]:
                    # continue with next i
                    break
                elif self[k] > self[j]:
                    r = i
                # else, continue comparing
        if r > 0:
            self.rotate(r)
        return self

    def num_loops(self):
        """Return the number of loops attached to this vertex."""
        if self._num_loops is None:
            seen = {}
            loops = 0
            for x in xrange(len(self)):
                if self[x] in seen:
                    loops += 1
                else:
                    seen[self[x]] = True
            self._num_loops = loops
        return self._num_loops


class Fatgraph(object):
    """A fully-decorated ribbon graph.

    Exports a (read-only) sequence interface, through which vertices
    can be accessed.
    """
    # the only reason to use `__slots__` here is to keep a record of
    # all instance attribute names.
##     __slots__ = [
##         '_fasteq_cache',
##         '_id',
##         '_id_factory',
##         '_numbering',
##         '_vertextype',
##         '_vertex_valences',
##         'endpoints_v',
##         'endpoints_i',
##         'numbering',
##         'num_edges',
##         'num_external_edges',
##         'num_vertices',
##         'edge_numbering',
##         'vertices',
##         ]

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
        if isinstance(g_or_vs, Fatgraph):
            # copy-constructor
            self._vertex_valences = g_or_vs._vertex_valences
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
                   " sequence of `Vertex` instances."

            #: Factory method to make a `Vertex` instance from a linear
            #  list of incident edge colorings.
            self._vertextype = vertextype

            #: list of vertices
            self.vertices = g_or_vs

            #: list of vertex valences 
            self._vertex_valences = tuple(sorted(kwargs.get('_vertex_valences',
                                                            (len(v) for v in g_or_vs))))

            #: Number of edge colors
            self.num_edges = kwargs.get('num_edges',
                                        sum(self._vertex_valences) / 2)

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

        assert self._ok()

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

    @cache_symmetric
    def __eq__(self, other):
        """Return `True` if Fatgraphs `self` and `other` are isomorphic.

        Examples::

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
        assert isinstance(other, Fatgraph), \
               "Fatgraph.__eq__:" \
               " called with non-Fatgraph argument `other`: %s" % other
        # shortcuts
        if self is other:
            return True
        elif ((self.num_edges != other.num_edges)
            or (self.num_vertices != other.num_vertices)
            or (self._vertex_valences != other._vertex_valences)):
            return False
        elif (self.vertices == other.vertices) \
                 and (self.endpoints_v == other.endpoints_v) \
                 and (self.endpoints_i == other.endpoints_i):
            return True
        else:
            # go the long way: try to find an explicit isomorphims
            # between graphs `self` and `other`
            try:
                # if there is any morphism, then return `True`
                self.isomorphisms(other).next()
                return True
            except StopIteration:
                # list of morphisms is empty, graphs are not equal.
                return False

    def __getitem__(self, index):
        return self.vertices[index]

    def __hash__(self):
        return self._persistent_id

    def __iter__(self):
        """Return iterator over vertices."""
        return iter(self.vertices)

    # both `__eq__` and `__ne__` are needed for testing equality of objects;
    # see `<http://www.voidspace.org.uk/python/articles/comparison.shtml>`
    def __ne__(self, other):
        """The opposite of `__eq__` (which see)."""
        return not self.__eq__(other)

    def __repr__(self):
        if hasattr(self, 'num_external_edges') and self.num_external_edges > 0:
            return "Fatgraph(%s, num_external_edges=%d)" \
                   % (repr(self.vertices), self.num_external_edges)
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

    @cache
    def boundary_components(self):
        """Return the number of boundary components of this `Fatgraph` object.

        Each boundary component is represented by the list of (colored)
        edges::

          >>> Fatgraph([Vertex([2,1,0]),Vertex([2,0,1])]).boundary_components()
          [CyclicTuple((2, 0)), CyclicTuple((1, 2)), CyclicTuple((0, 1))]

        If both sides of an edge belong to the same boundary
        component, that edge appears twice in the list::

          >>> Fatgraph([Vertex([2,1,1]),Vertex([2,0,0])]).boundary_components()
          [CyclicTuple((2, 0, 2, 1)), CyclicTuple((1,)), CyclicTuple((0,))]
          
          >>> Fatgraph([Vertex([2,1,0]),Vertex([2,1,0])]).boundary_components()
          [CyclicTuple((2, 1, 0, 2, 1, 0))]
          
        """
        assert self.num_external_edges == 0, \
               "Fatgraph.boundary_components: "\
               " cannot compute boundary components for" \
               " a graph with nonzero external edges: %s" % self
        
        # micro-optimizations
        L = self.num_edges
        endpoints_v = self.endpoints_v
        endpoints_i = self.endpoints_i

        # pass1: build a "copy" of `graph`, replacing each edge
        # coloring with a triplet `(other, index, edge)` pointing to
        # the other endpoint of that same edge: the element at
        # position `index` in vertex `other`.
        pass1 = []
        for (index_of_vertex_in_graph, vertex) in enumerate(self.vertices):
            replacement = []
            for (index_of_edge_in_vertex, edge) in enumerate(vertex):
                (v1, v2) = endpoints_v[edge]
                (i1, i2) = endpoints_i[edge]
                if v1 != v2:
                    if v1 == index_of_vertex_in_graph:
                        other_end = v2
                        other_index = i2
                    else:
                        other_end = v1
                        other_index = i1
                else:
                    other_end = v1 # == v2, that is *this* vertex
                    if index_of_edge_in_vertex == i1:
                        other_index = i2
                    else:
                        other_index = i1
                # replace other_index with index of *next* edge
                # (in the vertex cyclic order)
                if other_index == len(self.vertices[other_end])-1:
                    other_index = 0
                else:
                    other_index += 1
                replacement.append((other_end, other_index, edge))
            pass1.append(replacement)

        # pass2: now build a linear list, each element of the list
        # corresponding to an half-edge, of triples `(pos, seen,
        # edge)` where `pos` is the index in this list where the other
        # endpoint of that edge is located, `seen` is a flag, set to
        # `False` for half-edges that have not yet been walked
        # through, and `edge` is the corresponding edge.
        pass2 = []
        # build indices to the where each vertex begins in the linear list
        vi=[0]
        for vertex in self.vertices:
            vi.append(vi[-1]+len(vertex))
        # build list from collapsing the 2-level structure
        for vertex in pass1:
            for triplet in vertex:
                pass2.append([vi[triplet[0]]+triplet[1], False, triplet[2]])

        # pass3: pick up each element of the linear list, and follow it
        # until we come to an already marked one.
        result = []
        pos = 0
        while pos < len(pass2):
            # fast forward to an element that we've not yet seen
            while (pos < len(pass2)) and (pass2[pos][1] == True):
                pos += 1
            if pos >= len(pass2):
                break
            # walk whole chain of edges
            i = pos
            result.append([]) # new boundary component
            while pass2[i][1] == False:
                result[-1].append(pass2[i][2])
                pass2[i][1] = True
                i = pass2[i][0]
            pos += 1

        # consistency check: each edge must occur two times
        # (either twice in the same b.c., or in two distinct
        # components)
        if __debug__:
            cnt = [0] * self.num_edges
            for bc in result:
                for edge in bc:
                    cnt[edge] += 1
            for x in xrange(len(cnt)):
                assert cnt[x] == 2, \
                       "Fatgraph.boundary_components:"\
                       " edge %d occurs %d times "\
                       " in boundary components `%s`"\
                       " of graph `%s`"\
                       % (x, cnt[x], result, self)

        # that's all, folks!
        return list(CyclicTuple(bc) for bc in result)
        

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

        If boundary components have already been computed, they are
        adapted and set in the contracted graph too::

          >>> g1 = Fatgraph([Vertex([2,1,1]), Vertex([2,0,0])])
          >>> g1.boundary_components() # compute b.c.'s
          [CyclicTuple((2, 0, 2, 1)), CyclicTuple((1,)), CyclicTuple((0,))]
          >>> g2 = g1.contract(2)
          >>> g2.boundary_components()
          [CyclicTuple((1, 0)), CyclicTuple((1,)), CyclicTuple((0,))]

          >>> g1 = Fatgraph([Vertex([2,1,0]), Vertex([2,0,1])])
          >>> g1.boundary_components() # compute b.c.'s
          [CyclicTuple((2, 0)), CyclicTuple((1, 2)), CyclicTuple((0, 1))]
          >>> g2 = g1.contract(2)
          >>> g2.boundary_components()
          [CyclicTuple((1,)), CyclicTuple((0, 1)), CyclicTuple((0,))]

        In the above examples, notice that any reference to edge `2`
        has been removed from the boundary cycles after contraction.
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

        ## To keep code simpler we plug the higher-numbered vertex
        ## into the lower-numbered one, and possibly reverse the
        ## orientation on the resulting graph.
        
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
                 )
            assert g.num_boundary_components() == self.num_boundary_components()

        # build new graph 
        return Fatgraph(new_vertices,
                     vertextype = self._vertextype,
                     endpoints = (new_endpoints_v, new_endpoints_i),
                     num_edges = self.num_edges - 1,
                     num_external_edges = self.num_external_edges,
                     orientation = new_edge_numbering,
                     )
            

    @cache
    def genus(self):
        """Return the genus g of this `Fatgraph` object."""
        n = self.num_boundary_components()
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
          >>> g2 = g.hangcircle(0, 1)
          >>> g1 == Fatgraph([Vertex([0,1,2]), Vertex([3,2,1]), Vertex([0,3,4]), Vertex([5,5,4])])
          True
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
    
        
    def is_canonical(self):
        """Return `True` if this `Fatgraph` object is canonical.

        A graph is canonical iff:
        1) Each vertex is represented by the maximal sequence, among all
           sequences representing the same cyclic order.
        2) Vertices are sorted in lexicographic order.

        Examples::
          >>> Fatgraph([Vertex([2,1,0]), Vertex([2,1,0])]).is_canonical()
          True             
          >>> Fatgraph([Vertex([2,1,0]), Vertex([2,0,1])]).is_canonical()
          True             
          >>> Fatgraph([Vertex([2,0,1]), Vertex([2,1,0])]).is_canonical()
          False
          >>> Fatgraph([Vertex([0,1,2]), Vertex([2,1,0])]).is_canonical()
          False 
        """
        previous_vertex = None
        for vertex in self.vertices:
            if not vertex.is_canonical_representative():
                return False
            if previous_vertex and (previous_vertex < vertex):
                return False
            previous_vertex = vertex
        return True


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
    def isomorphisms(self, other):
        """Iterate over isomorphisms from `self` to `other`.

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

        Examples::

        """
        ## Compute all permutations of vertices that preserve valence.
        ## (A permutation `p` preserves vertex valence if vertices `v`
        ## and `p[v]` have the same valence.)
        vs1 = self.valence_spectrum()
        vs2 = other.valence_spectrum()
        
        # save valences as we have no guarantees that keys() method
        # will always return them in the same order
        vsk = vs1.keys()

        # graphs differ in vertex valences, no isomorphisms
        if not set(vsk) == set(vs2.keys()):
            return
        # graphs have unequal vertex distribution by valence, no isomorphisms
        if not  dict((val, len(vs1[val])) for val in vsk) \
               == dict((val, len(vs2[val])) for val in vsk):
            return

        src_indices = concat([ vs1[val] for val in vsk ])
        permutations_of_vertices_of_same_valence = [
            list(reversed([ p[:] for p in InplacePermutationIterator(vs2[val]) ]))
            for val in vsk
            ]
        dst_index_permutations = [
            concat(ps)
            for ps
            in SetProductIterator(permutations_of_vertices_of_same_valence)
            ]

        ## Build a list of vertex-to-vertex mappings; each of these
        ## mappings is a list of pairs `(i1, i2)`, where `i1` is a
        ## vertex index in `g1` and `i2` is a vertex index in `g2`.
        ## Vertices `i1` and `i2` have the same valence, and they must
        ## also agree on secondary invariants, like the number of
        ## loops.

        # use number of attached loops as a secondary invariant
        loops1 = [ v.num_loops() for v in self.vertices ]
        loops2 = [ v.num_loops() for v in other.vertices ]

        n = self.num_vertices
        candidate_pvs = []
        for dst_indices in dst_index_permutations:
            vertex_mapping_is_good_candidate = True
            # each `dsts` is a permutation of the vertex indices of `g2`
            for i in xrange(n):
                if loops1[src_indices[i]] != loops2[dst_indices[i]]:
                    # cannot map `src` into `dst`, proceed to next `i`
                    vertex_mapping_is_good_candidate = False
                    break
            if vertex_mapping_is_good_candidate:
                candidate_pvs.append(Permutation(dict((src_indices[i],
                                                       dst_indices[i])
                                                      for i in xrange(n))))

        ## Browse the list of vertex-to-vertex mappings; for each one, check:
        ##   1. that it induces a permutation on the edge colors;
        ##   2. that this permutation preserves the adjacency list;
        ##   3. if there is a numbering on the boundary cycles, that it is
        ##      preserved by the induced mapping.

        dst_valences = [ len(v) for v in other.vertices ]
        rots = [ range(val) for val in dst_valences ]
        for pv in candidate_pvs:
            # try mapping vertices with every possible rotation
            # vertex at index `i1` in `self` should be mapped to
            # vertex at index `i2` in `other` with shift `rot[i2]`
            # (that is, `self.vertices[i1][0:len]` should be mapped
            # linearly onto `other.vertices[i2][rot:rot+len]`).
            for rot in SetProductIterator(rots):
                pe = Permutation()
                pe_is_ok = True  # optimistic start
                for (src, dst) in pv.iteritems():
                    shift = rot[dst]
                    v1 = self.vertices[src]
                    v2 = other.vertices[dst][shift:shift+dst_valences[dst]]
                    if not pe.extend(v1, v2):
                        # cannot extend, proceed to next `rot`
                        pe_is_ok = False
                        break
                if pe_is_ok and (len(pe) > 0):
                    # Check that the combined action of `pv` and `pe`
                    # preserves the adjacency relation.  Note:
                    #   - we make list comprehensions of both adjacency lists
                    #     to avoid inverting `pe`: that is, we compare the
                    #     the adjacency lists in the order they have in `self`,
                    #     but with the vertex numbering from `other`;
                    #   - elements of the adjacency lists are made into
                    #     `set`s for unordered comparison;
                    if 0 != cmp([ set(other.endpoints_v[pe[x]])
                                  for x in xrange(other.num_edges) ],
                                [ set(pv.itranslate(self.endpoints_v[x]))
                                  for x in xrange(self.num_edges) ]):
                        continue # to next `rot`
                    yield (pv, rot, pe)


    @cache
    def num_boundary_components(self):
        """Return the number of boundary components of this `Fatgraph` object.

        Each boundary component is represented by the list of (colored)
        edges.

        Examples::
          >>> Fatgraph([Vertex([2,1,0]), Vertex([2,1,0])]).num_boundary_components()
          1
          >>> Fatgraph([Vertex([2,1,0]), Vertex([2,0,1])]).num_boundary_components()
          3
        """
        return len(self.boundary_components())


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

    @cache
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
        assert set(result.keys()) == set(self._vertex_valences), \
               "Fatgraph.valence_spectrum:" \
               "Computed valence spectrum `%s` does not exhaust all " \
               " vertex valences %s" \
               % (result, self._vertex_valences)
        assert set(concat(result.values())) \
               == set(range(self.num_vertices)), \
               "Fatgraph.valence_spectrum:" \
               "Computed valence spectrum `%s` does not exhaust all " \
               " %d vertex indices" % (result, self.num_vertices)
        return result


def MakeNumberedGraphs(graph):
    """Return all distinct (up to isomorphism) decorations of
    `Fatgraph` instance `graph` with a numbering of the boundary
    cycles.

    Examples::
    
      >>> ug1 = Fatgraph([Vertex([2,0,0]), Vertex([2,1,1])])
      >>> for g in MakeNumberedGraphs(ug1): print g
      NumberedFatgraph(Fatgraph([Vertex([2, 0, 0]), Vertex([2, 1, 1])]),    
                       numbering={CyclicTuple((2, 1, 2, 0)): 0,   
                                  CyclicTuple((1,)): 1,           
                                  CyclicTuple((0,)): 2})
      NumberedFatgraph(Fatgraph([Vertex([2, 0, 0]), Vertex([2, 1, 1])]),    
                       numbering={CyclicTuple((0,)): 0,           
                                  CyclicTuple((2, 1, 2, 0)): 1,   
                                  CyclicTuple((1,)): 2})        
      NumberedFatgraph(Fatgraph([Vertex([2, 0, 0]), Vertex([2, 1, 1])]),    
                       numbering={CyclicTuple((0,)): 0,          
                                  CyclicTuple((1,)): 1,
                                  CyclicTuple((2, 1, 2, 0)): 2})
       
    Note that, when only one numbering out of many possible ones is
    returned because of isomorphism, the returned numbering may not be
    the trivial one (it is infact the first permutation of 0..n
    returned by `InplacePermutationIterator`)::
      
      >>> ug2 = Fatgraph([Vertex([2,1,0]), Vertex([2,0,1])])
      >>> MakeNumberedGraphs(ug2)
      [NumberedFatgraph(Fatgraph([Vertex([2, 1, 0]), Vertex([2, 0, 1])]), 
                        numbering={CyclicTuple((2, 0)): 0,
                                   CyclicTuple((0, 1)): 1,     
                                   CyclicTuple((1, 2)): 2})]

    When the graph has only one boundary component, there is only one
    possible numbering, which is actually returned::
    
      >>> ug3 = Fatgraph([Vertex([1,0,1,0])])
      >>> MakeNumberedGraphs(ug3)
      [NumberedFatgraph(Fatgraph([Vertex([1, 0, 1, 0])]), 
                        numbering={CyclicTuple((1, 0, 1, 0)): 0})]
      
    """
    graphs = []
    bc = graph.boundary_components()
    n = len(bc) # == graph.num_boundary_components()

    for numbering in InplacePermutationIterator(range(n)):
        # make a copy of `graph` and add the given numbering
        g = NumberedFatgraph(graph,
                             [(numbering[x], bc[x]) for x in xrange(n)])

        # only add `g` to list if it is *not* isomorphic to a graph
        # already in the list
        if g not in graphs:
            graphs.append(g)

    return graphs


class NumberedFatgraph(Fatgraph):
    """A `Fatgraph` decorated with a numbering of the boundary components.

    A numbered fatgraph is constructed from a `Fatgraph` instance
    (called the *underlying graph*) and a numbering (that is, a
    bijective map assigning an integer to each boundary components of
    the underlying graph).

    Examples::

      >>> ug = Fatgraph([Vertex([1, 0, 1, 0])])  # underlying graph
      >>> ng = NumberedFatgraph(ug, numbering={CyclicTuple((1, 0, 1, 0)): 0})

    Since `NumberedFatgraphs` are just decorated `Fatgraphs`, they
    only differ in the way two `NumberedFatgraph` instances are deemed
    isomorphic:
      - they must have isomorphic underlying graphs;
      - the numberings must match under the isomorphism map.
    For example::
    
      >>> ug = Fatgraph([Vertex([1, 1, 0, 0])])
      >>> ng1 = NumberedFatgraph(ug, numbering={CyclicTuple((0,)): 1,   \
                                                CyclicTuple((1, 0)): 0, \
                                                CyclicTuple((1,)): 2})
      >>> ng2 = NumberedFatgraph(ug, numbering={CyclicTuple((0,)): 2,   \
                                                CyclicTuple((1, 0)): 1, \
                                                CyclicTuple((1,)): 0})
      >>> ng1 == ng2
      False
      """

##     __slots__ = [
##         '_numbering',
##         '_persistent_id',
##         'numbering',
##         'underlying',
##         ]


    def __init__(self, underlying, numbering, vertextype=Vertex):
        Fatgraph.__init__(self, underlying)
        self.underlying = underlying # XXX: = self (?)
        self.numbering = numbering

    @cache_symmetric
    def __eq__(self, other):
        """Return `True` if NumberedFatgraphs `self` and `other` are isomorphic.

        Fatgraph instances equipped with a numbering are compared as
        numbered graphs (that is, the isomorphism should transform the
        numbering on the source graph onto the numbering of the
        destination)::

          >>> NumberedFatgraph(Fatgraph([Vertex([2,0,1]), Vertex([2,1,0])]), 
          ...                  numbering=[(0, CyclicTuple((0,1))), 
          ...                             (1, CyclicTuple((0,2))), 
          ...                             (2, CyclicTuple((2,1))) ] ) 
          ... == NumberedFatgraph( 
          ...                  Fatgraph([Vertex([2,0,1]), Vertex([2,1,0])]), 
          ...                  numbering=[(0, CyclicTuple((1,0))), 
          ...                             (2, CyclicTuple((0,2))), 
          ...                             (1, CyclicTuple((2,1))) ])
          True

          >>> NumberedFatgraph(Fatgraph([Vertex([1, 0, 0, 2, 2, 1])]), 
          ...                  numbering=[(0, CyclicTuple((2,))), 
          ...                             (1, CyclicTuple((0,2,1))), 
          ...                             (3, CyclicTuple((0,))), 
          ...                             (2, CyclicTuple((1,))) ]) 
          ... == NumberedFatgraph( 
          ...                  Fatgraph([Vertex([2, 2, 1, 1, 0, 0])]), 
          ...                  numbering=[(0, CyclicTuple((2,))), 
          ...                             (1, CyclicTuple((0,))), 
          ...                             (3, CyclicTuple((2,1,0))), 
          ...                             (2, CyclicTuple((1,))) ])
          False
        
          >>> NumberedFatgraph(Fatgraph([Vertex([1, 0, 0, 2, 2, 1])]), 
          ...                  numbering=[(0, CyclicTuple((2,))), 
          ...                             (1, CyclicTuple((0,2,1))), 
          ...                             (3, CyclicTuple((0,))), 
          ...                             (2, CyclicTuple((1,))) ]) 
          ... == NumberedFatgraph( 
          ...                  Fatgraph([Vertex([2, 2, 1, 1, 0, 0])]), 
          ...                  numbering=[(3, CyclicTuple((2,))), 
          ...                             (0, CyclicTuple((0,))), 
          ...                             (2, CyclicTuple((2,1,0))), 
          ...                             (1, CyclicTuple((1,))) ])
          False

          >>> NumberedFatgraph(Fatgraph([Vertex([3, 2, 2, 0, 1]), Vertex([3, 1, 0])]), 
          ...                  numbering=[(0, CyclicTuple((2,))),  
          ...                             (1, CyclicTuple((0, 1))),  
          ...                             (2, CyclicTuple((3, 1))),  
          ...                             (3, CyclicTuple((0, 3, 2))) ]) 
          ... == NumberedFatgraph( 
          ...                  Fatgraph([Vertex([2, 3, 1]), Vertex([2, 1, 3, 0, 0])]), 
          ...                  numbering=[(0, CyclicTuple((0,))), 
          ...                             (2, CyclicTuple((1, 3))), 
          ...                             (3, CyclicTuple((3, 0, 2))), 
          ...                             (1, CyclicTuple((2, 1))) ])
          True

        Examples::

          >>> NumberedFatgraph.__eq__(NumberedFatgraph(
          ...                             Fatgraph([Vertex([0, 1, 2, 0, 2, 1])]),
          ...                             numbering={CyclicTuple((0,)): 1,
          ...                                        CyclicTuple((1, 2, 0, 1, 2)): 0}),
          ...                         NumberedFatgraph(
          ...                             Fatgraph([Vertex([0, 1, 2, 0, 2, 1])]),
          ...                             numbering={CyclicTuple((0,)): 0,
          ...                                        CyclicTuple((1, 2, 0, 1, 2)): 1}))
          False


        """
        assert isinstance(other, NumberedFatgraph), \
               "NumberedFatgraph.__eq__:" \
               " `other` argument is no `NumberedFatgraph` instance: %s" % other
        # shortcuts
        if self is other:
            return True
        elif ((self.underlying.num_edges != other.underlying.num_edges)
            or (self.underlying.num_vertices != other.underlying.num_vertices)
            or (self.underlying._vertex_valences != other.underlying._vertex_valences)):
            return False
        elif (self.underlying.vertices == other.underlying.vertices) \
                 and (self.underlying.endpoints_v == other.underlying.endpoints_v) \
                 and (self.underlying.endpoints_i == other.underlying.endpoints_i) \
                 and (self.numbering == other.numbering):
            return True
        else:
            # go the long way: try to find an explicit isomorphims
            # between graphs `self` and `other`
            try:
                # if there is any morphism, then return `True`
                self.isomorphisms(other).next()
                return True
            except StopIteration:
                # list of morphisms is empty, graphs are not equal.
                return False


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
        """Return new `NumberedFatgraph` obtained by contracting the specified edge."""
        # check that we are not contracting a loop or an external edge
        assert not self.is_loop(edgeno), \
               "NumberedFatgraph.contract: cannot contract a loop."
        assert (edgeno >= 0) and (edgeno < self.num_edges), \
               "NumberedFatgraph.contract: invalid edge number (%d):"\
               " must be in range 0..%d" \
               % (edgeno, self.num_edges)

        # Build new list of vertices, removing the contracted edge and
        # shifting all indices above:
        #   - edges numbered 0..edgeno-1 are unchanged;
        #   - edges numbered `edgeno+1`.. are renumbered, 
        #     shifting the number down one position;
        #   - edge `edgeno` is kept intact, will be removed by mating
        #     operation (see below).
        renumber_edges = dict((i+1,i)
                              for i in xrange(edgeno, self.num_edges))
        #   - edge `edgeno` is removed (subst with `None`)
        renumber_edges[edgeno] = None  
        numbering = [ (n, CyclicTuple(itranslate(renumber_edges, bcy)))
                      for (bcy, n) in self.numbering.iteritems() ]

        return NumberedFatgraph(self.underlying.contract(edgeno), numbering)
        

    @cache_iterator
    def isomorphisms(self, other):
        """Iterate over isomorphisms from `self` to `other`.

        See `Fatgraph.isomrphisms` for a discussion of the
        representation of isomorphisms and example usage.
        """
        for (pv, rot, pe) in Fatgraph.isomorphisms(self.underlying, other.underlying):
            pe_does_not_preserve_bc = False
            for bc1 in self.underlying.boundary_components():
                bc2 = CyclicTuple(pe.itranslate(bc1))
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


    def numbering_get(self):
        """Return the numbering previously set on this instance via
        `.set_numbering()`.
        """
        return self._numbering

    def _numbering_set(self, tuples):
        """Set the `.numbering` attribute from a sequence of tuples
        `(n, bcy)`.  Each `n` is a non-negative integer, and each
        `bcy` is an edge cycle representing a boundary component.

        The numbering attribute is set to a dictionary mapping the
        boundary cycle `bcy` to the integer `n`::

          >>> ug0 = Fatgraph([Vertex([1,2,0]), Vertex([1,0,2])])
          >>> bc = ug0.boundary_components()  # three b.c.'s
          >>> ng0 = NumberedFatgraph(ug0, enumerate(bc))
          >>> ng0.numbering_get()             \
              == { CyclicTuple((0, 2)): 2, \
                   CyclicTuple((1, 0)): 0, \
                   CyclicTuple((2, 1)): 1  }
          True

        When two boundary components are represented by the same edge
        cycle, that edge cycle is mapped to a `frozenset` instance
        containing the (distinct) indices assigned to it::
        
          >>> ug1 = Fatgraph([Vertex([1,2,0,1,2,0])])
          >>> bc = ug1.boundary_components()  # two b.c.'s
          >>> ng1 = NumberedFatgraph(ug1, [(0, bc[1]), (1, bc[0])])
          >>> ng1.numbering_get() \
              == {CyclicTuple((1, 2, 0)): frozenset([0, 1])}
          True
          
        By definition of a fatgraph, *at most* two boundary components
        may be represented by the same boundary cycle.

        """
        if tuples is None:
            self._numbering = None
            return
        # allow also `.numbering_set({ bcy: n, ...})`
        if isinstance(tuples, dict):
            tuples = tuples.iteritems()
        numbering = {}
        for (n,bcy) in tuples:
            if bcy in numbering:
                # a single boundary cycle can represent at most *two*
                # boundary components
                assert type(numbering[bcy]) is types.IntType, \
                       "Inconsistent numbering: "\
                       " Boundary cycle `%s` belongs to boundary"\
                       " component numbered %d, but already tagged %s" \
                       % (bcy, n, str.join(",",numbering[bcy]))
                numbering[bcy] = frozenset((n, numbering[bcy]))
            else:
                numbering[bcy] = n
        self._numbering = numbering

    numbering = property(numbering_get, _numbering_set)
    
    

class _ConnectedGraphsIterator(BufferingIterator):
    """Iterate over all connected numbered graphs having vertices of
    the prescribed valences.
    
    Examples::

      >>> for g in ConnectedGraphsIterator([4]): print g
      NumberedFatgraph(Fatgraph([Vertex([1, 0, 1, 0])]), 
                       numbering={CyclicTuple((1, 0, 1, 0)): 0})
      NumberedFatgraph(Fatgraph([Vertex([1, 1, 0, 0])]),
                       numbering={CyclicTuple((1, 0)): 0,
                                  CyclicTuple((0,)): 1,
                                  CyclicTuple((1,)): 2})    
      NumberedFatgraph(Fatgraph([Vertex([1, 1, 0, 0])]),
                       numbering={CyclicTuple((1,)): 0,
                                  CyclicTuple((1, 0)): 1,            
                                  CyclicTuple((0,)): 2})
      NumberedFatgraph(Fatgraph([Vertex([1, 1, 0, 0])]),                       
                       numbering={CyclicTuple((1,)): 0,
                                  CyclicTuple((0,)): 1,
                                  CyclicTuple((1, 0)): 2})
      
      >>> for g in ConnectedGraphsIterator([3,3]): print g
      NumberedFatgraph(Fatgraph([Vertex([2, 0, 1]), Vertex([2, 0, 1])]),          
                       numbering={CyclicTuple((2, 0, 1, 2, 0, 1)): 0}) 
      NumberedFatgraph(Fatgraph([Vertex([2, 1, 0]), Vertex([2, 0, 1])]), 
                       numbering={CyclicTuple((2, 0)): 0,
                                  CyclicTuple((0, 1)): 1,               
                                  CyclicTuple((1, 2)): 2})  
      NumberedFatgraph(Fatgraph([Vertex([2, 1, 1]), Vertex([2, 0, 0])]),          
                       numbering={CyclicTuple((2, 0, 2, 1)): 0,         
                                  CyclicTuple((0,)): 1,                 
                                  CyclicTuple((1,)): 2})              
      NumberedFatgraph(Fatgraph([Vertex([2, 1, 1]), Vertex([2, 0, 0])]), 
                       numbering={CyclicTuple((1,)): 0,
                                  CyclicTuple((2, 0, 2, 1)): 1,         
                                  CyclicTuple((0,)): 2})               
      NumberedFatgraph(Fatgraph([Vertex([2, 1, 1]), Vertex([2, 0, 0])]),
                       numbering={CyclicTuple((1,)): 0,
                                  CyclicTuple((0,)): 1,
                                  CyclicTuple((2, 0, 2, 1)): 2})

    Generation of all graphs with prescribed vertex valences `(v_1,
    v_2, ..., v_n)` goes this way:
    
      1) Generate all lists `L` of length `2*n` comprising the symbols
         `{0,...,n-1}`, each of which is repeated exactly twice;

      2) Pick such a list `L` and break it into smaller pieces of
         length `v_1`, ..., `v_n`, each one corresponding to a vertex
         (this is actually done in the `Fatgraph` class constructor),
         effectively building a graph `G`.

      3) Test the graph `G` for connectedness: if it's not connected,
         then go back to step 2).

      4) Compare `G` with all graphs previously found: if there is a
         permutation of the edge labels that transforms `G` into an
         already-found graph, then go back to step 2).

    """

##     __slots__ = [
##         '_graphs',
##         ]

    def __init__(self, vertex_valences, vertextype=VertexCache()):
        assert debug.is_sequence_of_integers(vertex_valences), \
               "ConnectedGraphsIterator: " \
               " argument `vertex_valences` must be a sequence of integers,"\
               " but got %s" % vertex_valences
        assert 0 == sum(vertex_valences) % 2, \
               "ConnectedGraphsIterator: " \
               " sum of vertex valences must be divisible by 2"

        self._graphs = GivenValenceGraphsIterator(vertex_valences,
                                                  vertextype=vertextype)

        # initialize superclass
        BufferingIterator.__init__(self)

    def refill(self):
        return MakeNumberedGraphs(self._graphs.next())
ConnectedGraphsIterator = persist.PersistedIterator(_ConnectedGraphsIterator)


class _GivenValenceGraphsIterator(object):
    """Iterate over all connected (un-numbered) ribbon graphs having
    vertices of the prescribed valences.
    
    Examples::

      >>> for g in GivenValenceGraphsIterator([4]): print g
      Fatgraph([Vertex([1, 0, 1, 0])])
      Fatgraph([Vertex([1, 1, 0, 0])])

      >>> for g in GivenValenceGraphsIterator([3,3]): print g
      Fatgraph([Vertex([2, 0, 1]), Vertex([2, 0, 1])])
      Fatgraph([Vertex([2, 1, 0]), Vertex([2, 0, 1])])
      Fatgraph([Vertex([2, 1, 1]), Vertex([2, 0, 0])])

    Generation of all graphs with prescribed vertex valences `(v_1,
    v_2, ..., v_n)` proceeds this way:
    
      1) Generate all lists `L` of length `2*n` comprising the symbols
         `{0,...,n-1}`, each of which is repeated exactly twice;

      2) Pick such a list `L` and break it into smaller pieces of
         length `v_1`, ..., `v_n`, each one corresponding to a vertex
         (this is actually done in the `Fatgraph` class constructor),
         effectively building a graph `G`.

      3) Test the graph `G` for connectedness: if it's not connected,
         then go back to step 2).

      4) Compare `G` with all graphs previously found: if there is a
         permutation of the edge labels that transforms `G` into an
         already-found graph, then go back to step 2).

    """

##     __slots__ = [
##         'graphs',
##         'vertextype',
##         '_edge_seq_iterator',
##         '_morphism_factory',
##         '_vertex_valences',
##         ]

    def __init__(self, vertex_valences, vertextype=VertexCache()):
        assert debug.is_sequence_of_integers(vertex_valences), \
               "GivenValenceGraphsIterator: " \
               " argument `vertex_valences` must be a sequence of integers,"\
               " but got %s" % vertex_valences
        assert 0 == sum(vertex_valences) % 2, \
               "GivenValenceGraphsIterator: " \
               " sum of vertex valences must be divisible by 2"

        self.vertextype = vertextype
        self.graphs = []
        self._morphism_factory = None
        self._vertex_valences = vertex_valences

        # build list [0,0,1,1,...,n-1,n-1]
        starting_edge_seq=[]
        for l in xrange(0, sum(vertex_valences)/2):
            starting_edge_seq += [l,l]
        self._edge_seq_iterator = InplacePermutationIterator(starting_edge_seq)

    def __iter__(self):
        return self
    
    def next(self):
        for edge_seq in self._edge_seq_iterator:
            # Break up `edge_seq` into smaller sequences corresponding
            # to vertices.
            vertices = []
            base = 0
            for current_vertex_index in xrange(len(self._vertex_valences)):
                VLEN = self._vertex_valences[current_vertex_index]
                vertices.append(self.vertextype(edge_seq[base:base+VLEN]))
                base += VLEN

            current = Fatgraph(vertices,
                            vertextype=self.vertextype,)
            if not (current.is_canonical() and current.is_connected()):
                continue
            
            if not current in self.graphs:
                self.graphs.append(current)
                return current
            # otherwise, continue with next `current` graph

        # no more graphs to generate
        raise StopIteration
GivenValenceGraphsIterator = persist.PersistedIterator(_GivenValenceGraphsIterator)


def AlgorithmB(n):
    """Iterate over all binary trees with `n+1` internal nodes in
    pre-order.  Or, equivalently, iterate over all full binary trees
    with `n+2` leaves.

    Returns a pair `(l,r)` of list, where `l[j]` and `r[j]` are the
    left and right child nodes of node `j`.  A `None` in `l[j]`
    (resp. `r[j]`) means that node `j` has no left (resp. right) child.

    The number of such trees is equal to the n-th Catalan number::

      >>> [ len(list(AlgorithmB(n))) for n in xrange(6) ]
      [1, 2, 5, 14, 42, 132]

    This is "Algorithm B" in Knuth's Volume 4, fasc. 4, section 7.2.1.6,
    with the only difference that node indices start from 0 here.
    """
    # B1 -- Initialize
    l = [ k+1 for k in xrange(n) ] + [None]
    r = [ None ] * (n+1)
    while True:
        # B2 -- Visit
        yield (l, r)
        # B3 -- Find `j`
        j = 0
        while l[j] == None:
            r[j] = None
            l[j] = j+1
            j += 1
            if j >= n:
                raise StopIteration
        # B4 -- Find `k` and `y`
        y = l[j]
        k = None
        while r[y] != None:
            k = y
            y = r[y]
        # B5 -- Promote `y`
        if k is not None:
            r[k] = None
        else:
            l[j] = None
        r[y] = r[j]
        r[j] = y


def Tree(nodeseq=[], vertextype=VertexCache()):
    """Construct a tree fatgraph from sequence of internal nodes.

    Items in `nodeseq` are sequences `(c[0], c[1], ..., c[n])` where
    `c[i]` are the labels (index number) of child nodes; if any
    `c[i]` is `None`, then a new terminal node is appended as child.
    The node created from the first item in `nodeseq` gets the label
    `0`, the second node gets the label `1`, and so on.
        
    Each internal node with `n` children is represented as a fatgraph
    vertex with `n+1` edges; the first one connects the node with its
    parent, and the other ones with the children, in the order they
    were given in the constructor.  Terminal nodes (i.e., leaves) are
    represented as edges with one loose end; that is, terminal nodes
    are *not* represented as vertices.

    Loose-end edges are given a negative index color, to easily
    distinguish them from regular edges, and are not counted in the
    `num_edges` attribute.

      >>> Tree([(1, 2), (None, None), (3, None), (None, None)])
      Fatgraph([Vertex([-1, 0, 1]), Vertex([0, -2, -3]),
             Vertex([1, 2, -4]), Vertex([2, -5, -6])],
             num_external_edges=6)

    """
    edge_to_parent = {}
    next_external_edge_label = -2  # grows downwards: -2,-3,...
    next_internal_edge_label = 0   # grows upwards: 1,2,...
    internal_edge_endpoints_v = []
    external_edge_endpoints_v = []
    internal_edge_endpoints_i = []
    external_edge_endpoints_i = []
    vertices = []
    next_vertex_index = 0
    for cs in nodeseq:
        # Edges incident to this vertex; for the root vertex, the
        # connection to the parent node is just the first external
        # edge (labeled `-1`)
        edgeno = edge_to_parent.get(next_vertex_index, -1)
        edges = [ edgeno ]
        if edgeno < 0:
            # external edge (only happens on first iteration)
            external_edge_endpoints_v.append([next_vertex_index, None])
            external_edge_endpoints_i.append([0, None])
        else:
            # internal edge, add this vertex as second endpoint
            internal_edge_endpoints_v[edgeno].append(next_vertex_index)
            internal_edge_endpoints_i[edgeno].append(0)
        for child in cs:
            if child is None:
                # terminal node here
                edges.append(next_external_edge_label)
                external_edge_endpoints_v.append([next_vertex_index, None])
                external_edge_endpoints_i.append([len(edges)-1, None])
                next_external_edge_label -= 1
            else:
                # internal node here
                edges.append(next_internal_edge_label)
                edge_to_parent[child] = next_internal_edge_label
                internal_edge_endpoints_v.append([next_vertex_index])
                internal_edge_endpoints_i.append([len(edges)-1])
                next_internal_edge_label += 1
        vertices.append(vertextype(edges))
        next_vertex_index += 1

    return Fatgraph(vertices,
                 endpoints = (internal_edge_endpoints_v +
                                list(reversed(external_edge_endpoints_v)),
                              internal_edge_endpoints_i +
                                list(reversed(external_edge_endpoints_i))),
                 num_edges = next_internal_edge_label,
                 num_external_edges = -next_external_edge_label-1,
                 vertextype=vertextype)


class TreeIterator(BufferingIterator):
    """Iterate over trees with a specified number of leaves.

    Internal nodes are allowed to have any number of children: the
    iterator is not restricted to binary trees.
    """

    def __init__(self, num_leaves):
        self._internal_edges = num_leaves - 2

        self._trees = [ Tree(zip(l,r))
                        for l,r in AlgorithmB(num_leaves - 2) ]

        BufferingIterator.__init__(self, self._trees)

    def refill(self):
        if self._internal_edges > 1:
            self._internal_edges -= 1
            l = self._internal_edges - 1 # label of edge to contract
            new_trees = [ t.contract(l) for t in self._trees ]
            self._trees = new_trees
            return new_trees
        else:
            raise StopIteration


def MgnGraphsInsertionGenerator(g,n):
    """Iterate over all connected trivalent fatgraphs having the
    prescribed genus `g` and number of boundary cycles `n`.
    
    Examples::

      >>> for g in MgnGraphsInsertionGenerator(0,3): print g
      Fatgraph([Vertex([1, 2, 1]), Vertex([2, 0, 0])]) 
      Fatgraph([Vertex([1, 0, 2]), Vertex([2, 0, 1])])

      >>> for g in MgnGraphsInsertionGenerator(1,1): print g
      Fatgraph([Vertex([1, 0, 2]), Vertex([2, 1, 0])])

    """

    assert n > 0, \
           "MgnGraphsInsertionGenerator: " \
           " number of boundary cycles `n` must be positive,"\
           " but got `%s` instead" % n
    assert (g > 0) or (g == 0 and n >= 3), \
           "MgnGraphsInsertionGenerator: " \
           " Invalid (g,n) pair (%d,%d): "\
           " need either g>0 or g==0 and n>2" \
           % (g,n)

    #: Unique (up to isomorphism) graphs found so far
    graphs = []

    #: Minimum number of edges of a (g,n)-graph
    max_valence = 2 * (2*g + n - 1)

    ## pass 1: Gather all roses.
    logging.debug("  MgnGraphsInsertionGenerator: Computing roses with %d leaves ...", max_valence/2)
    roses = []
    discarded = 0
    for rose in GivenValenceGraphsIterator((max_valence,)):
        if (rose.genus() != g) \
               or (rose.num_boundary_components() != n) \
               or (rose in roses):
            discarded += 1
            continue
        roses.append(rose)
        # a rose is a valid fatgraph too
        #graphs.extend(MakeNumberedGraphs(rose))
    logging.debug("    MgnGraphsInsertionGenerator: Found %d distinct unique roses; discarded %d.",
                 len(roses), discarded)

    ## pass 2: Gather all 3-valent graphs.
    trivalent = []
    #: Full binary trees
    logging.debug("  MgnGraphsInsertionGenerator: Computing full binary trees with %d leaves ...",
                 max_valence - 3)
    trees = [ Tree(zip(l,r))
              for l,r in AlgorithmB(max_valence - 3) ]

    logging.debug("  MgnGraphsInsertionGenerator: Computing trivalent fat graphs ...")
    discarded = 0
    for rose in roses:
        # now substitute the unique vertex with any possible tree
        # and any possible rotation
        for places in xrange(max_valence):
            # need to make a deep copy, because `Vertex` objects are shared
            rotated_rose = Fatgraph([copy(rose[0])])
            #rotated_rose = rose
            rotated_rose[0].rotate(places)
            for tree in trees:
                graph = rotated_rose.graft(tree, 0)
                if (graph.genus() != g) \
                       or (graph.num_boundary_components() != n) \
                       or (graph in trivalent):
                    discarded += 1
                    continue
                trivalent.append(graph)
    logging.debug("    MgnGraphsInsertionGenerator: Found %d distinct trivalent graphs, discarded %d.",
                 len(trivalent), discarded)

    return iter(trivalent)

#MgnGraphsInsertionGenerator = persist.PersistedIterator(_MgnGraphsInsertionGenerator)


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
                for x in xrange(G.num_edges):
                    yield G.hangcircle(x,0)
                    yield G.hangcircle(x,1)

            logging.debug("  MgnTrivalentGraphsRecursiveGenerator(%d,%d): "
                          "pass 2: bridge all edges of a single graph in M_{%d,%d} ..." % (g,n, g,n-1))
            for G in MgnTrivalentGraphsRecursiveGenerator(g,n-1):
                for x in xrange(G.num_edges):
                    # since G.bridge() is symmetric, we need only consider
                    # edge pairs `(x,y)` where `y <= x`.
                    for y in xrange(x+1):
                        yield G.bridge(x,0, y,0)
                        yield G.bridge(x,0, y,1)
                        yield G.bridge(x,1, y,0)
                        yield G.bridge(x,1, y,1)

            logging.debug("  MgnTrivalentGraphsRecursiveGenerator(%d,%d): "
                          "pass 3: bridge all edges of a single graph in M_{%d,%d} ..." % (g,n, g-1,n+1))
            for G in MgnTrivalentGraphsRecursiveGenerator(g-1,n+1):
                for x in xrange(G.num_edges):
                    # since G.bridge() is symmetric, we need only consider
                    # edge pairs `(x,y)` where `y <= x`.
                    for y in xrange(x+1):
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
                            for x1 in xrange(G1.num_edges):
                                for x2 in xrange(G2.num_edges):
                                    yield Fatgraph.bridge2(G1, x1, 0, G2, x2, 0)
                                    yield Fatgraph.bridge2(G1, x1, 0, G2, x2, 1)
                                    yield Fatgraph.bridge2(G1, x1, 1, G2, x2, 0)
                                    yield Fatgraph.bridge2(G1, x1, 1, G2, x2, 1)

        unique = []
        for G in graphs(g,n):
            # XXX: should this check be done in graphs(g,n)?
            if (G.genus(), G.num_boundary_components()) != (g,n):
                continue
            if G in unique:
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

##     __slots__ = [
##         '_batch',
##         '_current_edge',
##         '_num_vertices',
##         '_vertextype',
##         'g',
##         'n',
##         ]

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
        self._batch = next_batch
        self._num_vertices -= 1

        logging.debug("  Found %d distinct unique graphs with %d vertices, discarded %d.",
                     len(next_batch), self._num_vertices, discarded)
        return next_batch

#MgnGraphsIterator = persist.PersistedIterator(_MgnGraphsIterator)



class MgnNumberedGraphsIterator(BufferingIterator):
    """Iterate over all connected numbered fatgraphs having the
    prescribed genus `g` and number of boundary cycles `n`.
    
    Examples::

      >>> for g in MgnNumberedGraphsIterator(0,3): print g
      NumberedFatgraph(Fatgraph([Vertex([1, 2, 1]), Vertex([2, 0, 0])]),
                       numbering={CyclicTuple((1,)): 0,
                                  CyclicTuple((0,)): 1,
                                  CyclicTuple((2, 0, 2, 1)): 2})
      NumberedFatgraph(Fatgraph([Vertex([1, 2, 1]), Vertex([2, 0, 0])]),
                       numbering={CyclicTuple((2, 0, 2, 1)): 0,
                                  CyclicTuple((1,)): 1,
                                  CyclicTuple((0,)): 2})
      NumberedFatgraph(Fatgraph([Vertex([1, 2, 1]), Vertex([2, 0, 0])]),
                       numbering={CyclicTuple((0,)): 0,
                                  CyclicTuple((2, 0, 2, 1)): 1,
                                  CyclicTuple((1,)): 2})
      NumberedFatgraph(Fatgraph([Vertex([1, 0, 2]), Vertex([2, 0, 1])]),
                       numbering={CyclicTuple((1, 2)): 0,
                                  CyclicTuple((2, 0)): 1,
                                  CyclicTuple((0, 1)): 2})
      NumberedFatgraph(Fatgraph([Vertex([1, 1, 0, 0])]),
                       numbering={CyclicTuple((1, 0)): 0,
                                  CyclicTuple((0,)): 1,
                                  CyclicTuple((1,)): 2})
      NumberedFatgraph(Fatgraph([Vertex([1, 1, 0, 0])]),
                       numbering={CyclicTuple((1,)): 0,
                                  CyclicTuple((1, 0)): 1,
                                  CyclicTuple((0,)): 2})
      NumberedFatgraph(Fatgraph([Vertex([1, 1, 0, 0])]),
                       numbering={CyclicTuple((1,)): 0,
                                  CyclicTuple((0,)): 1,
                                  CyclicTuple((1, 0)): 2})

      >>> for g in MgnNumberedGraphsIterator(1,1): print g
      NumberedFatgraph(Fatgraph([Vertex([1, 0, 2]), Vertex([2, 1, 0])]),
                       numbering={CyclicTuple((1, 0, 2, 1, 0, 2)): 0})
      NumberedFatgraph(Fatgraph([Vertex([1, 0, 1, 0])]),
                       numbering={CyclicTuple((1, 0, 1, 0)): 0})

    """

    def __init__(self, g, n, vertextype=VertexCache()):
        #self.__naked_graphs_iterator = MgnGraphsIterator(g, n)
        self.__naked_graphs_iterator = MgnGraphsIterator(g, n)
        BufferingIterator.__init__(self)

    def refill(self):
        return MakeNumberedGraphs(self.__naked_graphs_iterator.next())

#MgnNumberedGraphsIterator = persist.PersistedIterator(_MgnNumberedGraphsIterator)



## main: run tests

if "__main__" == __name__:
    import doctest
    doctest.testmod(name='rg',
                    optionflags=doctest.NORMALIZE_WHITESPACE)
