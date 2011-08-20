#! /usr/bin/env python
#
"""Classes and functions to deal with ribbon graphs.
"""
__docformat__ = 'reStructuredText'


import debug, pydb, sys
sys.excepthook = pydb.exception_hook


from combinatorics import (
    InplacePermutationIterator,
    SetProductIterator,
    Permutation,
    )
from cyclicseq import CyclicList,CyclicTuple
from utils import (
    BufferingIterator,
    concat,
    itranslate
    )

from copy import copy
from itertools import chain,count,izip


class VertexCache(object):
    """A caching factory of `Vertex` objects.
    """
    __slots__ = [
        'cache',
        ]
    def __init__(self):
        self.cache = {}
    def __call__(self, edge_seq):
        key = tuple(edge_seq)
        if key not in self.cache:
            self.cache[key] = Vertex(key)
        return self.cache[key]


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

    __slots__ = [ '_num_loops' ]

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
        """Return `True` if this `Vertex` object is maximal among
        representatives of the same cyclic sequence.
        
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
                
class Graph(object):
    """A fully-decorated ribbon graph.

    Exports a (read-only) sequence interface, through which vertices
    can be accessed.
    """
    # the only reason to use `__slots__` here is to keep a record of
    # all instance attribute names.
    __slots__ = [
        '_boundary_components',
        '_fasteq_cache',
        '_genus',
        '_id',
        '_id_factory',
        '_numbering',
        '_num_boundary_components',
        '_valence_spectrum',
        '_vertextype',
        '_vertex_valences',
        'endpoints_v',
        'endpoints_i',
        'numbering',
        'num_edges',
        'num_external_edges',
        'num_vertices',
        'orient_v',
        'orient_a',
        'vertices',
        ]

    def __init__(self, vertices, vertextype=Vertex,
                 __fasteq_cache={}, __id_factory=count(), **kwargs):
        """Construct a `Graph` instance, taking list of vertices.

        Argument `vertices` must be a sequence of `Vertex` class
        instances::  

          >>> G1 = Graph([Vertex([2,0,1]), Vertex([2,1,0])])

        Note that the list of vertices is assigned, *not copied* into
        the instance variable.

        """
        assert debug.is_sequence_of_type(Vertex, vertices), \
               "Graph.__init__: parameter `vertices` must be" \
               " sequence of `Vertex` instances."

        #: class-wide cache for `__eq__` results
        self._fasteq_cache = __fasteq_cache

        #: unique numeric id of this `Graph` object; id's of deleted
        #: object should *not* be re-used
        self._id = __id_factory.next()
        self._id_factory = __id_factory
        
        # the following values will be computed on-demand

        #: Cached boundary cycles of this graph; see
        #  `.boundary_components()` for an explanation of the
        #  format. This is initially `None` and is actually computed on
        #  first invocation of the `.boundary_components()` method.
        self._boundary_components = kwargs.get('_boundary_components', None)

        #: Cached genus; initially `None`, and actually computed on
        #  first invocation of the `.genus()` method.
        self._genus = kwargs.get('_genus', None)

        #: cached number of boundary cycles of this graph; this is
        #  initially `None` and is actually computed on first invocation
        #  of the `.num_boundary_components()` method.
        self._num_boundary_components = kwargs.get('_num_boundary_components',
                                                   None)

        #: Valence spectrum; initially `None`, and actually computed on
        #  first invocation of the `.valence_spectrum()` method.
        self._valence_spectrum = kwargs.get('_valence_spectrum', None)

        #: Factory method to make a `Vertex` instance from a linear
        #  list of incident edge colorings.
        self._vertextype = vertextype
        
        #: list of vertices
        self.vertices = vertices
        
        #: list of vertex valences 
        self._vertex_valences = tuple(sorted(kwargs.get('_vertex_valences',
                                                        (len(v) for v in vertices))))

        #: Number of edge colors
        self.num_edges = kwargs.get('num_edges',
                                    sum(self._vertex_valences) / 2)

        #: Number of external (loose-end) edges
        self.num_external_edges = kwargs.get('num_external_edges', 0)

        #: Number of vertices
        self.num_vertices = kwargs.get('num_vertices', len(self.vertices))
        
        #: Order on the boundary cycles, or `None`.
        self.numbering = kwargs.get('numbering', None)

        if 'endpoints' in kwargs:
            (self.endpoints_v, self.endpoints_i) = kwargs.get('endpoints')
        else:
            #: Adjacency list of this graph.  For each edge, store a pair
            #: `(v1, v2)` where `v1` and `v2` are (indices of)
            #: endpoint vertices of an edge, and a corresponding pair
            #: `(i1, i2)` where `i1` and `i2` are indices of the given
            #: edge at vertices `v1` and `v2`.
            #: Each pair `(v, i)` represents a flag by the endpoint
            #: vertex and the index of the edge in the vertex.  (The
            #: vertex index alone is not enough for representing the
            #: edge arrow for loops.)
            self.endpoints_v = [ [] for dummy in xrange(self.num_edges) ]
            self.endpoints_i = [ [] for dummy in xrange(self.num_edges) ]
            for current_vertex_index in xrange(len(self.vertices)):
                for (edge_index_in_vertex, edge) \
                        in enumerate(self.vertices[current_vertex_index]):
                    assert edge in range(self.num_edges), \
                               "Graph.__init__:"\
                               " edge number %d not in range 0..%d" \
                               % (edge, self.num_edges)
                    self.endpoints_v[edge].append(current_vertex_index)
                    self.endpoints_i[edge].append(edge_index_in_vertex)

        ## Orientation is represented according to Conant+Vogtmann:
        ##   - an ordering of the vertices;
        ##   - an arrow on each edge.
        ## Unless an explicit ordering is given in the constructor arguments,
        ## use the obvious one: vertex have the order they are stored in
        ## `self.vertices` and all edges get a positive arrow.
        if 'orientation' in kwargs:
            self.orient_v, self.orient_a = kwargs.get('orientation')
        else:
            #: Order on the vertices: the mapping from (internal)
            #: vertex index to the vertex number given according to
            #: the order.
            self.orient_v = [ x for x in xrange(self.num_vertices) ]
            #: Arrow on the edges; this is `-1` iff the arrow goes
            #: opposite to the order the edge endpoints are stored in
            #: the adjacency list, and `+1` otherwise.
            self.orient_a = [ +1 for x in xrange(self.num_edges) ]

        assert self._ok()

    def _ok(self):
        """Perform coherency checks on internal state variables of
        `Graph` instance and return `True` if they all pass.
        """
        assert self.num_edges > 0, \
               "Graph `%s` has 0 edges." % (self)
        # check regular edges endpoints
        for (edge, ep_v, ep_i) in izip(count(),
                                       self.endpoints_v[:self.num_edges],
                                       self.endpoints_i[:self.num_edges]):
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
                   "Graph `%s`:"\
                   " invalid attachment indices `%s`" \
                   " for endpoints %s of regular edge %d"\
                   % (self, ep_i, ep_v, edge)
##                    "Graph `%s` has invalid regular endpoints array `%s/%s`" \
##                    " invalid endpoints pair %s/%s for edge %d" \
##                    % (self, self.endpoints_v, self.endpoints_i,
##                       ep_v, ep_i, edge)
            assert (edge in self.vertices[ep_v[0]])
            assert (edge in self.vertices[ep_v[1]])
##                    "Invalid endpoints %s for edge %d of graph `%s`" \
##                    % (ep_v, edge, self)
        # check external edges endpoints
        for (edge, ep_v, ep_i) in izip(count(),
                                       self.endpoints_v[self.num_edges:],
                                       self.endpoints_i[self.num_edges:]):
            xedge = -self.num_external_edges + edge
            assert (ep_v[1] is None)
            assert isinstance(ep_v[0], int)
            assert (ep_i[1] is None)
            assert isinstance(ep_i[0], int)
            assert (0 <= ep_v[0] < self.num_vertices)
            assert (0 <= ep_i[0] < len(self.vertices[ep_v[0]]))
##                    "Graph `%s` has invalid external endpoints array: `%s/%s`" \
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

        assert self.orient_v is not None
        assert self.orient_a is not None
        
        return True

    def __eq__(self, other):
        """Return `True` if Graphs `self` and `other` are isomorphic.

        Examples::

          >>> Graph([Vertex([1,0,0,1])]) == Graph([Vertex([1,1,0,0])])
          True

          >>> Graph([Vertex([2,0,0]), Vertex([2,1,1])]) \
                == Graph([Vertex([2,2,0]), Vertex([1,1,0])])
          True

          >>> Graph([Vertex([2,0,1]), Vertex([2,0,1])]) \
                == Graph([Vertex([2,1,0]), Vertex([2,0,1])])
          False

          >>> Graph([Vertex([2,0,1]), Vertex([2,0,1])]) \
                == Graph([Vertex([2,0,0]), Vertex([2,1,1])])
          False

          >>> Graph([Vertex([2,0,0]), Vertex([2,1,1])]) \
                == Graph([Vertex([1,1,0,0])])
          False

        Graph instances equipped with a numbering are compared as
        numbered graphs (that is, the isomorphism should transform the
        numbering on the source graph onto the numbering of the
        destination)::

          >>> Graph([Vertex([2,0,1]), Vertex([2,1,0])], \
                     numbering=[(0, CyclicTuple((0,1))), \
                                (1, CyclicTuple((0,2))), \
                                (2, CyclicTuple((2,1))) ] ) \
              == Graph([Vertex([2,0,1]), Vertex([2,1,0])], \
                        numbering=[(0, CyclicTuple((1,0))), \
                                   (2, CyclicTuple((0,2))), \
                                   (1, CyclicTuple((2,1))) ])
          True

          >>> Graph([Vertex([1, 0, 0, 2, 2, 1])], \
                     numbering=[(0, CyclicTuple((2,))), \
                                (1, CyclicTuple((0,2,1))), \
                                (3, CyclicTuple((0,))), \
                                (2, CyclicTuple((1,))) ]) \
                == Graph([Vertex([2, 2, 1, 1, 0, 0])], \
                          numbering=[(0, CyclicTuple((2,))), \
                                     (1, CyclicTuple((0,))), \
                                     (3, CyclicTuple((2,1,0))), \
                                     (2, CyclicTuple((1,))) ])
          False
        
          >>> Graph([Vertex([1, 0, 0, 2, 2, 1])], \
                     numbering=[(0, CyclicTuple((2,))), \
                                (1, CyclicTuple((0,2,1))), \
                                (3, CyclicTuple((0,))), \
                                (2, CyclicTuple((1,))) ]) \
                == Graph([Vertex([2, 2, 1, 1, 0, 0])], \
                          numbering=[(3, CyclicTuple((2,))), \
                                     (0, CyclicTuple((0,))), \
                                     (2, CyclicTuple((2,1,0))), \
                                     (1, CyclicTuple((1,))) ])
          False

          >>> Graph([Vertex([3, 2, 2, 0, 1]), Vertex([3, 1, 0])], \
                    numbering=[(0, CyclicTuple((2,))),  \
                               (1, CyclicTuple((0, 1))),  \
                               (2, CyclicTuple((3, 1))),  \
                               (3, CyclicTuple((0, 3, 2))) ]) \
              == Graph([Vertex([2, 3, 1]), Vertex([2, 1, 3, 0, 0])], \
                       numbering=[(0, CyclicTuple((0,))), \
                                  (2, CyclicTuple((1, 3))), \
                                  (3, CyclicTuple((3, 0, 2))), \
                                  (1, CyclicTuple((2, 1))) ])
          True

        Examples::

          >>> Graph.__eq__(Graph([Vertex([0, 1, 2, 0, 2, 1])],
          ...                    numbering={CyclicTuple((0,)): 1,
          ...                               CyclicTuple((1, 2, 0, 1, 2)): 0}),
          ...              Graph([Vertex([0, 1, 2, 0, 2, 1])],
          ...                    numbering={CyclicTuple((0,)): 0,
          ...                               CyclicTuple((1, 2, 0, 1, 2)): 1}))
          False


        """
        assert isinstance(other, Graph), \
               "Graph.__eq__:" \
               " called with non-Graph argument `other`: %s" % other
        # try cached result first
        args = frozenset([self._id, other._id])
        if args not in self._fasteq_cache:
            # shortcuts
            if (self is other) or (self._id == other._id):
                self._fasteq_cache[args] = True
            elif ((self.num_edges != other.num_edges)
                or (self.num_vertices != other.num_vertices)
                or (self._vertex_valences != other._vertex_valences)):
                self._fasteq_cache[args] = False
            elif (self.vertices == other.vertices) \
               and (self.endpoints_v == other.endpoints_v) \
               and (self.endpoints_i == other.endpoints_i) \
               and (self.numbering == other.numbering):
                self._fasteq_cache[args] = True
            else:
                # go the long way: try to find an explicit isomorphims
                # between graphs `self` and `other`
                try:
                    # if there is any morphism, then return `True`
                    self.isomorphisms(other).next()
                    self._fasteq_cache[args] = True
                except StopIteration:
                    # list of morphisms is empty, graphs are not equal.
                    self._fasteq_cache[args] = False
        return self._fasteq_cache[args]

    def __getitem__(self, index):
        return self.vertices[index]

    def __hash__(self):
        return self._id

    def __iter__(self):
        """Return iterator over vertices."""
        return iter(self.vertices)

    # both `__eq__` and `__ne__` are needed for testing equality of objects;
    # see `<http://www.voidspace.org.uk/python/articles/comparison.shtml>`
    def __ne__(self, other):
        """The opposite of `__eq__` (which see)."""
        return not self.__eq__(other)

    def __repr__(self):
        # the hairy if-clause down here prints an attribute iff:
        #   - it is not `None` (which may happen for both `.numbering`
        #     and `.num_external_edges`), and
        #   - if it is an integer, it is not 0 (which may happen
        #     for `.num_external_edges`.
        extra = dict((x, getattr(self, x))
                     for x in ['numbering', 'num_external_edges']
                     if ((getattr(self, x) is not None)
                         and (not isinstance(getattr(self, x), int)
                              or (getattr(self, x) > 0))))
        return "Graph(%s%s)" % (repr(self.vertices),
                                  "".join((", %s=%s" % (k,v) for k,v
                                             in extra.iteritems())))
    
    def __str__(self):
        return repr(self)

    def automorphisms(self):
        """Enumerate automorphisms of this `Graph` object.

        See `.isomorphisms()` for details of how a `Graph`
        isomorphism is represented.
        """
        return self.isomorphisms(self)

    
    def boundary_components(self):
        """Return the number of boundary components of this `Graph` object.

        Each boundary component is represented by the list of (colored)
        edges::

          >>> Graph([Vertex([2,1,0]),Vertex([2,0,1])]).boundary_components()
          [CyclicTuple((2, 0)), CyclicTuple((1, 2)), CyclicTuple((0, 1))]

        If both sides of an edge belong to the same boundary
        component, that edge appears twice in the list::

          >>> Graph([Vertex([2,1,1]),Vertex([2,0,0])]).boundary_components()
          [CyclicTuple((2, 0, 2, 1)), CyclicTuple((1,)), CyclicTuple((0,))]
          
          >>> Graph([Vertex([2,1,0]),Vertex([2,1,0])]).boundary_components()
          [CyclicTuple((2, 1, 0, 2, 1, 0))]
          
        """
        assert self.num_external_edges == 0, \
               "Graph.boundary_components: "\
               " cannot compute boundary components for" \
               " a graph with nonzero external edges: %s" % self
        
        # if no cached result, compute it now...
        if self._boundary_components is None:
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
                           "Graph.boundary_components:"\
                           " edge %d occurs %d times "\
                           " in boundary components `%s`"\
                           " of graph `%s`"\
                           % (x, cnt[x], result, self)

            # save result for later reference
            self._boundary_components = list(CyclicTuple(bc) for bc in result)
        
        # that's all, folks!
        return self._boundary_components

    def clone(self):
        """Return a new `Graph` instance, sharing all attribute values
        with this one, except for the following ones:

          `_id`
            unique for each instance

          `_id_factory`
            shared among all instances

        """
        return Graph(self.vertices,
                     endpoints = (self.endpoints_v, self.endpoints_i),
                     num_edges = self.num_edges,
                     num_vertices = self.num_vertices,
                     numbering = self.numbering,
                     orientation = (self.orient_v, self.orient_a),
                     vertextype = self._vertextype,
                     _boundary_components = self._boundary_components,
                     _genus = self._genus,
                     _num_boundary_components = self._num_boundary_components,
                     _valence_spectrum = self._valence_spectrum,
                     _vertex_valences = self._vertex_valences,
                     )

    def _cmp_orient(self, other, iso):
        pv = iso[0]
        pe = iso[2]
        
        # compute permutation on vertex order
        image_orient_v = pv.rearrange(self.orient_v[:])
        result = Permutation(dict((image_orient_v[x], other.orient_v[x])
                                  for x in xrange(self.num_vertices))).sign()

        # for positively-oriented edges, return +1 if endpoint
        # indices are in decreasing lexicographic order, -1
        # otherwise; sign is reversed if edge has negative
        # orientation.
        def arrow(edge, endpoints_v, endpoints_i, orient_a):
            e0 = endpoints_v[edge][0]
            e1 = endpoints_v[edge][1]
            if e0 == e1:
                # use attachment indices to determine arrow
                e0 = endpoints_i[edge][0]
                e1 = endpoints_i[edge][1]
            assert e0 != e1
            if e0 > e1:
                return orient_a[edge]
            else:
                return -orient_a[edge]

        orig_arrows = [ arrow(x, other.endpoints_v, other.endpoints_i,
                              other.orient_a)
                        for x in xrange(other.num_edges) ]
        image_arrows = [ arrow(x,
                               [ (pv[self.endpoints_v[pe[y]][0]],
                                  pv[self.endpoints_v[pe[y]][1]])
                                 for y in xrange(self.num_edges)],
                               pe.rearrange(self.endpoints_i[:]),
                               pe.rearrange(self.orient_a[:]))
                         for x in xrange(self.num_edges) ]
        for x in xrange(self.num_edges):
            result *= orig_arrows[x] * image_arrows[x]

        return result


    def contract(self, edgeno):
        """Return new `Graph` obtained by contracting the specified edge.

        Examples::

          >>> Graph([Vertex([2,2,0]), Vertex([0,1,1])]).contract(0)
          Graph([Vertex([1, 1, 0, 0])])
          >>> Graph([Vertex([2,1,0]), Vertex([2,0,1])]).contract(1)
          Graph([Vertex([0, 1, 1, 0])])

        The M_{1,1} trivalent graph yield the same result no matter
        what edge is contracted::

          >>> Graph([Vertex([2,1,0]), Vertex([2,1,0])]).contract(0)
          Graph([Vertex([1, 0, 1, 0])])
          >>> Graph([Vertex([2,1,0]), Vertex([2,1,0])]).contract(1)
          Graph([Vertex([0, 1, 0, 1])])
          >>> Graph([Vertex([2,1,0]), Vertex([2,1,0])]).contract(2)
          Graph([Vertex([1, 0, 1, 0])])

        If boundary components have already been computed, they are
        adapted and set in the contracted graph too::

          >>> g1 = Graph([Vertex([2,1,1]), Vertex([2,0,0])])
          >>> g1.boundary_components() # compute b.c.'s
          [CyclicTuple((2, 0, 2, 1)), CyclicTuple((1,)), CyclicTuple((0,))]
          >>> g2 = g1.contract(2)
          >>> g2.boundary_components()
          [CyclicTuple((0, 1)), CyclicTuple((1,)), CyclicTuple((0,))]

          >>> g1 = Graph([Vertex([2,1,0]), Vertex([2,0,1])])
          >>> g1.boundary_components() # compute b.c.'s
          [CyclicTuple((2, 0)), CyclicTuple((1, 2)), CyclicTuple((0, 1))]
          >>> g2 = g1.contract(2)
          >>> g2.boundary_components()
          [CyclicTuple((0,)), CyclicTuple((1,)), CyclicTuple((0, 1))]

        In the above examples, notice that any reference to edge `2`
        has been removed from the boundary cycles after contraction.

        """
        # check that we are not contracting a loop or an external edge
        assert self.endpoints_v[edgeno][0] != self.endpoints_v[edgeno][1], \
               "Graph.contract: cannot contract a loop."
        assert (self.endpoints_v[edgeno][0] is not None) \
               and (self.endpoints_v[edgeno][1] is not None), \
               "Graph.contract: cannot contract an external edge."
        assert (edgeno >= 0) and (edgeno < self.num_edges), \
               "Graph.contract: invalid edge number (%d):"\
               " must be in range 0..%d" \
               % (edgeno, self.num_edges)

        # store endpoints and arrow on the edge-to-be-contracted, and
        # possibly swap endpoints so that `v1 < v2`
        (v1, v2) = self.endpoints_v[edgeno]
        (pos1, pos2) = self.endpoints_i[edgeno]
        if v1 > v2:
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
        # 1. Rotate endpoints `v1`, `v2` so that the given edge
        #    appears *last* in `v1` and *first* in `v2` (*Note:*
        #    the contracted edge has already been deleted from
        #    `v1` and `v2`, so index positions need to be adjusted);
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

        # Properly orient contracted graph.

        # vertex order stays the same, but vertices numbered higher
        # than `self.orient_v[v2]` need to be shifted down one place.
        cut = self.orient_v[v2]
        renumber_orient_v = { cut:None }
        renumber_orient_v.update(dict((x,x) for x in xrange(cut)))
        renumber_orient_v.update(dict((x+1,x)
                                      for x in xrange(cut,self.num_vertices)))
        new_orient_v = [ renumber_orient_v[self.orient_v[x]]
                         for x in xrange(self.num_vertices)
                         if x != v2 ]
        
        # edges keep their arrows
        new_orient_a = self.orient_a[:edgeno] + self.orient_a[edgeno+1:]
        
        # if the contracted edge had a negative sign, then swap
        # orientation on the first edge of contracted graph, to
        # compensate
        new_orient_a[0] *= self.orient_a[edgeno]

        numbering = None
        bc = None
        #   - edge `edgeno` is removed (subst with `None`)
        renumber_edges[edgeno] = None  
        if self.numbering is not None:
            numbering = [ (n, CyclicTuple(itranslate(renumber_edges, bcy)))
                          for (bcy, n) in self.numbering.iteritems() ]
            bc = [ bcy for (n,bcy) in numbering ]
        elif self._boundary_components is not None:
            bc = [ CyclicTuple(itranslate(renumber_edges, bcy))
                   for bcy in self._boundary_components ]

        # consistency check
        if __debug__:
            assert len(new_endpoints_v) == self.num_edges - 1
            assert len(new_endpoints_i) == len(new_endpoints_v)
            assert len(new_orient_v) == self.num_vertices - 1
            assert len(new_orient_a) == self.num_edges - 1
            for x in xrange(self.num_edges - 1):
                assert new_orient_a[x] in [+1, -1]
            for x in xrange(self.num_vertices - 1):
                assert 0 <= new_orient_v[x] < self.num_vertices - 1
            g = Graph(new_vertices,
                 vertextype = self._vertextype,
                 endpoints = (new_endpoints_v, new_endpoints_i),
                 num_edges = self.num_edges - 1,
                 num_external_edges = self.num_external_edges,
                 )
            assert g.num_boundary_components() == self.num_boundary_components()
            if numbering is not None:
                assert set(bcy for (n,bcy) in numbering) \
                       == set(g.boundary_components())

        # build new graph 
        return Graph(new_vertices,
                     vertextype = self._vertextype,
                     endpoints = (new_endpoints_v, new_endpoints_i),
                     num_edges = self.num_edges - 1,
                     num_external_edges = self.num_external_edges,
                     numbering = numbering,
                     orientation = (new_orient_v, new_orient_a),
                     _boundary_components = bc,
                     )
            

    def genus(self):
        """Return the genus g of this `Graph` object."""
        # compute value if not already done
        if (self._genus is None):
            n = self.num_boundary_components()
            K = self.num_vertices
            L = self.num_edges
            # by Euler, K-L+n=2-2*g
            self._genus = (L - K - n + 2) / 2
        return self._genus


    def graft(self, G, v):
        """Return new `Graph` formed by grafting graph `G` into vertex
        with index `v`.  The number of"external" edges in `G` must match the
        valence of `v`.
        """
        assert G.num_external_edges == len(self.vertices[v]), \
               "Graph.graft:" \
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

        return Graph(new_vertices, vertextype=vertextype,
                     num_edges = self.num_edges + G.num_edges,
                     num_external_edges = self.num_external_edges)


    def is_canonical(self):
        """Return `True` if this `Graph` object is canonical.

        A graph is canonical iff:
        1) Each vertex is represented by the maximal sequence, among all
           sequences representing the same cyclic order.
        2) Vertices are sorted in lexicographic order.

        Examples::
          >>> Graph([Vertex([2,1,0]), Vertex([2,1,0])]).is_canonical()
          True             
          >>> Graph([Vertex([2,1,0]), Vertex([2,0,1])]).is_canonical()
          True             
          >>> Graph([Vertex([2,0,1]), Vertex([2,1,0])]).is_canonical()
          False
          >>> Graph([Vertex([0,1,2]), Vertex([2,1,0])]).is_canonical()
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
          >>> Graph([Vertex([3, 3, 0, 0]), Vertex([2, 2, 1, 1])]).is_connected()
          False
          >>> Graph([Vertex([3, 1, 2, 0]), Vertex([3, 0, 2, 1])]).is_connected()
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
        this `Graph` instance."""
        return (-1 == Graph._cmp_orient(self, self, automorphism))


    def is_oriented(self):
        """Return `True` if `Graph` is orientable.

        A ribbon graph is orientable iff it has no
        orientation-reversing automorphism.

        Enumerate all automorphisms of `graph`, end exits with `False`
        result as soon as one orientation-reversing one is found.

        Examples::

          >>> Graph([Vertex([1,0,1,0])]).is_oriented()
          True

          >>> Graph([Vertex([2, 0, 1]), Vertex([2, 0, 1])]).is_oriented()
          True
          
          >>> Graph([Vertex([2, 1, 0]), Vertex([2, 0, 1])]).is_oriented()
          True
          
          >>> Graph([Vertex([2, 1, 1]), Vertex([2, 0, 0])]).is_oriented()
          True

          >>> Graph([Vertex([3, 2, 2, 0, 1]), Vertex([3, 1, 0])], \
                    numbering=[(0, CyclicTuple((2,))),  \
                               (1, CyclicTuple((0, 1))),  \
                               (2, CyclicTuple((3, 1))),  \
                               (3, CyclicTuple((0, 3, 2))) ]) \
                               .is_oriented()
          True
          >>> Graph([Vertex([2, 3, 1]), Vertex([2, 1, 3, 0, 0])], \
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

    def isomorphisms(self, other):
        """Iterate over isomorphisms from `self` to `other`.

        An isomorphism is represented by a tuple `(pv, rot, pe)` where:

          - `pv` is a permutation of ther vertices: the `i`-th vertex
            of `g1` is sent to the `pv[i]`-th vertex of `g2`, rotated
            by `rot[i]` places leftwards;

          - `pe` is a permutation of the edge colors: edge `i` in `g1`
            is mapped to edge `pe[i]` in `g2`.

        This method can iterate over the automorphism group of a
        graph::

          >>> g1 = Graph([Vertex([2, 1, 1]), Vertex([2, 0, 0])])
          >>> for f in g1.isomorphisms(g1): print f
          ({0: 0, 1: 1}, [0, 0], {0: 0, 1: 1, 2: 2})
          ({0: 1, 1: 0}, [0, 0], {0: 1, 1: 0, 2: 2})

        Or it can find the isomorphisms between two given graphs::

          >>> g2 = Graph([Vertex([2, 2, 0]), Vertex([1, 1, 0])])
          >>> for f in g1.isomorphisms(g2): print f
          ({0: 0, 1: 1}, [2, 2], {0: 1, 1: 2, 2: 0})
          ({0: 1, 1: 0}, [2, 2], {0: 2, 1: 1, 2: 0})

        If there are no isomorphisms connecting the two graphs, then no
        item is returned by the iterator::

          >>> g3 = Graph([Vertex([2, 1, 0]), Vertex([2, 0, 1])])
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

        assert set(vsk) == set(vs2.keys()), \
               "Graph.isomorphisms: "\
               " graphs `%s` and `%s` differ in vertex valences: `%s` vs `%s`" \
               % (self, other, vsk, vs2.keys())
        assert dict((val, len(vs1[val])) for val in vsk) \
               == dict((val, len(vs2[val])) for val in vsk), \
               "Graph.isomorphisms: graphs `%s` and `%s`" \
               " have unequal vertex distribution by valence: `%s` vs `%s`" \
               % (self, other, vs1, vs2)

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
                    if self.numbering is not None:
                        assert other.numbering is not None, \
                               "Graph.isomorphisms: " \
                               "Numbered and un-numbered graphs mixed in arguments."
                        pe_does_not_preserve_bc = False
                        for bc1 in self.boundary_components():
                            bc2 = CyclicTuple(pe.itranslate(bc1))
                            # there are cases (see examples in the
                            # `Graph.__eq__` docstring, in which the
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
                            continue # to next `rot`
                    yield (pv, rot, pe)

    def numbering_get(self):
        """Return the numbering previously set on this instance via
        `.set_numbering()`.
        """
        return self._numbering
    def numbering_set(self, tuples):
        """Set the `.numbering` attribute from a sequence of tuples
        `(n, bcy)`.  Each `n` is a non-negative integer, and each
        `bcy` is an edge cycle representing a boundary component.

        The numbering attribute is set to a dictionary mapping the
        boundary cycle `bcy` to the integer `n`::

          >>> g0 = Graph([Vertex([1,2,0]), Vertex([1,0,2])])
          >>> bc = g0.boundary_components()  # three b.c.'s
          >>> g0.numbering_set(enumerate(bc))
          >>> g0.numbering_get()             \
              == { CyclicTuple((0, 2)): 2, \
                   CyclicTuple((1, 0)): 0, \
                   CyclicTuple((2, 1)): 1  }
          True

        When two boundary components are represented by the same edge
        cycle, that edge cycle is mapped to a `frozenset` instance
        containing the (distinct) indices assigned to it::
        
          >>> g1 = Graph([Vertex([1,2,0,1,2,0])])
          >>> bc = g1.boundary_components()  # two b.c.'s
          >>> g1.numbering_set([(0, bc[1]), (1, bc[0])])
          >>> g1.numbering_get() \
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
                numbering[bcy] = frozenset((n, numbering[bcy]))
            else:
                numbering[bcy] = n
        self._numbering = numbering
        if __debug__:
            if self._boundary_components is not None:
                assert set(self._boundary_components) \
                       == set(self.numbering), \
                       "Graph.numbering_set:"\
                       " Not all boundary components were numbered"
                indices = set()
                for n in self._numbering.itervalues():
                    if isinstance(n, frozenset):
                        indices = set.union(indices, n)
                    else: # `n` is integer
                        indices.add(n)
                assert indices == set(range(len(self._boundary_components))), \
                       "Graph.numbering_set:"\
                       " Numbering indices `%s` are not a permutation"\
                       " of the range 0..%d (no. of boundary components)"\
                       % (self._numbering.values(), len(self._numbering),)
    numbering = property(numbering_get, numbering_set)
    
    def num_boundary_components(self):
        """Return the number of boundary components of this `Graph` object.

        Each boundary component is represented by the list of (colored)
        edges.

        Examples::
          >>> Graph([Vertex([2,1,0]), Vertex([2,1,0])]).num_boundary_components()
          1
          >>> Graph([Vertex([2,1,0]), Vertex([2,0,1])]).num_boundary_components()
          3
        """
        # compute boundary components and cache result
        if self._num_boundary_components is None:
            self._num_boundary_components = len(self.boundary_components())

        return self._num_boundary_components


    def projection(self, other):
        """Return the component of the projection of `self` on the
        basis vector `other`.  This can be either 0 (if `self` and
        `other` are not isomorphic), or +1/-1 depending on comparison
        of the orientation of `self` with the pull-back orientation on
        `other`.

        If the two graphs are not isomorphic, then the result is 0::

          >>> g1 = Graph([Vertex([0,1,2]), Vertex([0,2,1])])
          >>> g2 = Graph([Vertex([0,1,2]), Vertex([0,1,2])])
          >>> Graph.projection(g1, g2)
          0

        Any graph obviously projects onto itself with coefficient `1`::

          >>> Graph.projection(g1, g1)
          1

        Flipping the orientation on an edge reverses the coefficient
        sign::

          >>> g2 = g1.clone()
          >>> Graph.projection(g1, g2)
          1
          >>> g2.endpoints_v = [tuple(reversed(g2.endpoints_v[0]))] \
                                + g1.endpoints_v[1:]
          >>> g2.endpoints_i = [tuple(reversed(g2.endpoints_i[0]))] \
                                + g1.endpoints_i[1:]
          >>> Graph.projection(g1, g2)
          -1
          
        The same happens if vertex order is changed by an odd
        permutation::
        
          >>> g2 = g1.clone()
          >>> g2.orient_v = list(reversed(g2.orient_v))
          >>> Graph.projection(g1, g2)
          -1

        """
        assert isinstance(other, Graph), \
               "Graph.__eq__:" \
               " called with non-Graph argument `other`: %s" % other
        try:
            iso = self.isomorphisms(other).next()
            return Graph._cmp_orient(self, other, iso)
        except StopIteration:
            # list of morphisms is empty, graphs are not equal.
            return 0
    
    def valence_spectrum(self):
        """Return a dictionary mapping valences into vertex indices.

        Examples::

           >>> Graph([Vertex([1,1,0,0])]).valence_spectrum()
           {4: [0]}

           >>> Graph([Vertex([1,1,0]), Vertex([2,2,0])]).valence_spectrum()
           {3: [0, 1]}

           >>> Graph([Vertex([3, 1, 0, 1]), \
                      Vertex([4, 4, 0]), Vertex([3, 2, 2])]).valence_spectrum()
           {3: [1, 2], 4: [0]}
        """
        # compute spectrum on first invocation
        if self._valence_spectrum is None:
            self._valence_spectrum = {}
            for (index, vertex) in enumerate(self.vertices):
                l = len(vertex)
                if l in self._valence_spectrum:
                    self._valence_spectrum[l].append(index)
                else:
                    self._valence_spectrum[l] = [index]
            # consistency checks
            assert set(self._valence_spectrum.keys()) == set(self._vertex_valences), \
                   "Graph.valence_spectrum:" \
                   "Computed valence spectrum `%s` does not exhaust all " \
                   " vertex valences %s" \
                   % (self._valence_spectrum, self._vertex_valences)
            assert set(concat(self._valence_spectrum.values())) \
                   == set(range(self.num_vertices)), \
                   "Graph.valence_spectrum:" \
                   "Computed valence spectrum `%s` does not exhaust all " \
                   " %d vertex indices" % (self._valence_spectrum, self.num_vertices)
        return self._valence_spectrum


def MakeNumberedGraphs(graph):
    """Return all distinct (up to isomorphism) decorations of `graph`
    with a numbering of the boundary cycles.

    Examples::

      >>> g1 = Graph([Vertex([2,0,0]), Vertex([2,1,1])])
      >>> for g in MakeNumberedGraphs(g1): print g
      Graph([Vertex([2, 0, 0]), Vertex([2, 1, 1])],    
            numbering={CyclicTuple((2, 1, 2, 0)): 0,   
                       CyclicTuple((0,)): 2,           
                       CyclicTuple((1,)): 1})
      Graph([Vertex([2, 0, 0]), Vertex([2, 1, 1])],    
            numbering={CyclicTuple((2, 1, 2, 0)): 1,   
                       CyclicTuple((0,)): 0,           
                       CyclicTuple((1,)): 2})        
      Graph([Vertex([2, 0, 0]), Vertex([2, 1, 1])],    
             numbering={CyclicTuple((2, 1, 2, 0)): 2,  
                        CyclicTuple((0,)): 0,          
                        CyclicTuple((1,)): 1})
       
    Note that, when only one numbering out of many possible ones is
    returned because of isomorphism, the returned numbering may not be
    the trivial one (it is actually the first permutation of 0..n
    returned by `InplacePermutationIterator`)::
      
      >>> g2 = Graph([Vertex([2,1,0]), Vertex([2,0,1])])
      >>> MakeNumberedGraphs(g2)
      [Graph([Vertex([2, 1, 0]), Vertex([2, 0, 1])], 
              numbering={CyclicTuple((0, 1)): 1,     
                         CyclicTuple((1, 2)): 2,     
                         CyclicTuple((2, 0)): 0})]

    When the graph has only one boundary component, there is only one
    possible numbering, which is actually returned::
    
      >>> g3 = Graph([Vertex([1,0,1,0])])
      >>> MakeNumberedGraphs(g3)
      [Graph([Vertex([1, 0, 1, 0])], 
              numbering={CyclicTuple((1, 0, 1, 0)): 0})]
      
    """
    graphs = []
    bc = graph.boundary_components()
    n = len(bc) # == graph.num_boundary_components()

    for numbering in InplacePermutationIterator(range(n)):
        # make a copy of `graph` and add the given numbering
        g = graph.clone()
        g.numbering_set((numbering[x], bc[x]) for x in xrange(n))

        # only add `g` to list if it is *not* isomorphic to a graph
        # already in the list
        if g not in graphs:
            graphs.append(g)

    return graphs


class ConnectedGraphsIterator(BufferingIterator):
    """Iterate over all connected numbered graphs having vertices of
    the prescribed valences.
    
    Examples::

      >>> for g in ConnectedGraphsIterator([4]): print g
      Graph([Vertex([1, 0, 1, 0])],                       
             numbering={CyclicTuple((1, 0, 1, 0)): 0})
      Graph([Vertex([1, 1, 0, 0])],                       
            numbering={CyclicTuple((0,)): 1,              
                       CyclicTuple((1, 0)): 0,
                       CyclicTuple((1,)): 2})    
      Graph([Vertex([1, 1, 0, 0])],                       
            numbering={CyclicTuple((0,)): 2,              
                       CyclicTuple((1, 0)): 1,            
                       CyclicTuple((1,)): 0})
      Graph([Vertex([1, 1, 0, 0])],                       
            numbering={CyclicTuple((0,)): 1,              
                       CyclicTuple((1, 0)): 2,            
                       CyclicTuple((1,)): 0})

      >>> for g in ConnectedGraphsIterator([3,3]): print g
      Graph([Vertex([2, 0, 1]), Vertex([2, 0, 1])],          
            numbering={CyclicTuple((2, 0, 1, 2, 0, 1)): 0}) 
      Graph([Vertex([2, 1, 0]), Vertex([2, 0, 1])],          
            numbering={CyclicTuple((0, 1)): 1,               
                       CyclicTuple((1, 2)): 2,
                       CyclicTuple((2, 0)): 0})  
      Graph([Vertex([2, 1, 1]), Vertex([2, 0, 0])],          
            numbering={CyclicTuple((2, 0, 2, 1)): 0,         
                       CyclicTuple((0,)): 1,                 
                       CyclicTuple((1,)): 2})              
      Graph([Vertex([2, 1, 1]), Vertex([2, 0, 0])],          
            numbering={CyclicTuple((2, 0, 2, 1)): 1,         
                       CyclicTuple((0,)): 2,                 
                       CyclicTuple((1,)): 0})               
      Graph([Vertex([2, 1, 1]), Vertex([2, 0, 0])],          
            numbering={CyclicTuple((2, 0, 2, 1)): 2,         
                       CyclicTuple((0,)): 1,                 
                       CyclicTuple((1,)): 0})

    Generation of all graphs with prescribed vertex valences `(v_1,
    v_2, ..., v_n)` goes this way:
    
      1) Generate all lists `L` of length `2*n` comprising the symbols
         `{0,...,n-1}`, each of which is repeated exactly twice;

      2) Pick such a list `L` and break it into smaller pieces of
         length `v_1`, ..., `v_n`, each one corresponding to a vertex
         (this is actually done in the `Graph` class constructor),
         effectively building a graph `G`.

      3) Test the graph `G` for connectedness: if it's not connected,
         then go back to step 2).

      4) Compare `G` with all graphs previously found: if there is a
         permutation of the edge labels that transforms `G` into an
         already-found graph, then go back to step 2).

    """

    __slots__ = [
        '_graphs',
        ]

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


class GivenValenceGraphsIterator(object):
    """Iterate over all connected (un-numbered) ribbon graphs having
    vertices of the prescribed valences.
    
    Examples::

      >>> for g in GivenValenceGraphsIterator([4]): print g
      Graph([Vertex([1, 0, 1, 0])])
      Graph([Vertex([1, 1, 0, 0])])

      >>> for g in GivenValenceGraphsIterator([3,3]): print g
      Graph([Vertex([2, 0, 1]), Vertex([2, 0, 1])])
      Graph([Vertex([2, 1, 0]), Vertex([2, 0, 1])])
      Graph([Vertex([2, 1, 1]), Vertex([2, 0, 0])])

    Generation of all graphs with prescribed vertex valences `(v_1,
    v_2, ..., v_n)` proceeds this way:
    
      1) Generate all lists `L` of length `2*n` comprising the symbols
         `{0,...,n-1}`, each of which is repeated exactly twice;

      2) Pick such a list `L` and break it into smaller pieces of
         length `v_1`, ..., `v_n`, each one corresponding to a vertex
         (this is actually done in the `Graph` class constructor),
         effectively building a graph `G`.

      3) Test the graph `G` for connectedness: if it's not connected,
         then go back to step 2).

      4) Compare `G` with all graphs previously found: if there is a
         permutation of the edge labels that transforms `G` into an
         already-found graph, then go back to step 2).

    """

    __slots__ = [
        'graphs',
        'vertextype',
        '_edge_seq_iterator',
        '_morphism_factory',
        '_vertex_valences',
        ]

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

            current = Graph(vertices,
                            vertextype=self.vertextype,)
            if not (current.is_canonical() and current.is_connected()):
                continue
            
            if not current in self.graphs:
                self.graphs.append(current)
                return current
            # otherwise, continue with next `current` graph

        # no more graphs to generate
        raise StopIteration


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


def Tree(nodeseq=[], vertextype=Vertex):
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
      Graph([Vertex([-1, 0, 1]), Vertex([0, -2, -3]),
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

    return Graph(vertices,
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


class MgnGraphsIterator(BufferingIterator):
    """Iterate over all connected numbered graphs having the
    prescribed genus `g` and number of boundary cycles `n`.
    
    Examples::

      >>> for g in MgnGraphsIterator(0,3): print g
      Graph([Vertex([1, 2, 1]), Vertex([2, 0, 0])],    
            numbering={CyclicTuple((2, 0, 2, 1)): 2,   
                       CyclicTuple((0,)): 1,   
                       CyclicTuple((1,)): 0})
      Graph([Vertex([1, 2, 1]), Vertex([2, 0, 0])],    
            numbering={CyclicTuple((2, 0, 2, 1)): 0,   
                       CyclicTuple((0,)): 2,
                       CyclicTuple((1,)): 1})
      Graph([Vertex([1, 2, 1]), Vertex([2, 0, 0])],    
            numbering={CyclicTuple((2, 0, 2, 1)): 1, 
                       CyclicTuple((0,)): 0,   
                       CyclicTuple((1,)): 2}) 
      Graph([Vertex([1, 0, 2]), Vertex([2, 0, 1])],    
            numbering={CyclicTuple((2, 0)): 1,         
                       CyclicTuple((0, 1)): 2,         
                       CyclicTuple((1, 2)): 0})
      Graph([Vertex([1, 1, 0, 0])],                    
            numbering={CyclicTuple((0,)): 1,         
                       CyclicTuple((0, 1)): 2,         
                       CyclicTuple((1,)): 0})      
      Graph([Vertex([1, 1, 0, 0])],                    
            numbering={CyclicTuple((0,)): 2,         
                       CyclicTuple((0, 1)): 0,         
                       CyclicTuple((1,)): 1})       
      Graph([Vertex([1, 1, 0, 0])],                    
            numbering={CyclicTuple((0,)): 0,         
                       CyclicTuple((0, 1)): 1,         
                       CyclicTuple((1,)): 2})

      >>> for g in MgnGraphsIterator(1,1): print g
      Graph([Vertex([1, 0, 2]), Vertex([2, 1, 0])],          
            numbering={CyclicTuple((1, 0, 2, 1, 0, 2)): 0})
      Graph([Vertex([1, 0, 1, 0])],                          
            numbering={CyclicTuple((0, 1, 0, 1)): 0})

    """

    __slots__ = [
        '_batch',
        '_current_edge',
        '_num_vertices',
        '_vertextype',
        'g',
        'n',
        ]

    def __init__(self, g, n, vertextype=VertexCache()):
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
        
        #: Factory method to build `Vertex` instances from the
        #  incoming edges list.
        self._vertextype = vertextype

        #: Unique (up to isomorphism) graphs found so far
        graphs = []
        
        #: Minimum number of edges of a (g,n)-graph
        max_valence = 2 * (2*g + n - 1)

        ## pass 1: Gather all roses.
        roses = []
        for rose in GivenValenceGraphsIterator((max_valence,)):
            if (rose.genus() != self.g) \
                   or (rose.num_boundary_components() != self.n) \
                   or (rose in roses):
                continue
            roses.append(rose)
            # a rose is a valid fatgraph too
            #graphs.extend(MakeNumberedGraphs(rose))
            
        ## pass 2: Gather all 3-valent graphs.
        trivalent = []
        #: Full binary trees
        trees = [ Tree(zip(l,r))
                  for l,r in AlgorithmB(max_valence - 3) ]
        for rose in roses:
            # now substitute the unique vertex with any possible tree
            # and any possible rotation
            for places in xrange(max_valence):
                # need to make a deep copy, because `Vertex` objects are shared
                rotated_rose = Graph([copy(rose[0])])
                rotated_rose[0].rotate(places)
                for tree in trees:
                    graph = rotated_rose.graft(tree, 0)
                    if (graph.genus() != self.g) \
                           or (graph.num_boundary_components() != self.n) \
                           or (graph in trivalent):
                        continue
                    trivalent.append(graph)
                    # insert decorated graphs into iterator buffer
                    graphs.extend(MakeNumberedGraphs(graph))

        #: Graphs to be contracted at next `.refill()` invocation
        self._batch = trivalent

        #: Number of edges of graphs that will be returned by next
        #  `.refill()` call.  Starts with `6*g + 3*n - 7`, which is the
        #  highest-numbered edge in trivalent graphs.
        self._current_edge = 6*g + 3*n - 7

        self._num_vertices = 4*g + 2*n - 4
        
        # initialize superclass with list of roses + trivalent graphs
        BufferingIterator.__init__(self, graphs)

    def refill(self):
        if self._num_vertices == 0:
            raise StopIteration
        
        result = []
        next_batch = []
        for graph in self._batch:
            # contract all edges
            for edge in xrange(graph.num_edges):
                if not graph.is_loop(edge):
                    dg = graph.contract(edge)
                    if dg not in next_batch:
                        # put decorated version into `result`
                        result.extend(MakeNumberedGraphs(dg))
                        # put graph back into next batch for processing
                        next_batch.append(dg)
        self._batch = next_batch
        self._num_vertices -= 1
        return result
    


## main: run tests

if "__main__" == __name__:
    import doctest
    doctest.testmod(name='rg',
                    optionflags=doctest.NORMALIZE_WHITESPACE)
