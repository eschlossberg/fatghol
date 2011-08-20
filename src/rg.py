#! /usr/bin/env python
#
"""Classes and functions to deal with ribbon graphs.
"""
__docformat__ = 'reStructuredText'


import debug

from combinatorics import (
    InplacePermutationIterator,
    SetProductIterator,
    Permutation
    )
from cyclicseq import CyclicList
from utils import (
    concat,
    deep_cmp,
    itranslate
    )

from copy import copy
from itertools import chain,count,izip
import operator


class VertexCache(object):
    """A caching factory of `Vertex` objects.
    """
    __slots__ = ('cache',)
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

    def __repr__(self):
        return "Vertex(%s)" % CyclicList.__repr__(self)
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


class Graph(object):
    """A fully-decorated ribbon graph.

    Exports a (read-only) sequence interface, through which vertices
    can be accessed.
    """
    # the only reason to use `__slots__` here is to keep a record of
    # all instance attribute names.
    __slots__ = (
        '_boundary_components',
        '_genus',
        '_num_boundary_components',
        '_valence_spectrum',
        '_vertex_factory',
        '_vertex_valences',
        'edge_seq',
        'edge_seq_aliases',
        'endpoints',
        'num_edges',
        'num_vertices',
        'vertices',
        )

    ## *Note:* the constructor for this class is overloaded.  We could
    ## move the code making a `Graph` object out of vertex valences
    ## and edge colorings to a separate factory function
    ## `make_graph_from_edge_seq`, but then the constructor would need
    ## to recompute `self.edge_seq_aliases` and `self._vertex_valences`.
    ## Thus, we favor performance over code elegance and choose to
    ## overload the ctor.
    def __init__(self, vertices, edge_seq=None, vertex_factory=Vertex):
        """Construct a `Graph` instance, taking either list of
        vertices or the linear list of edges plus list of vertex
        valences.

        If argument `edge_seq` is `None`, then `vertices` must be a
        sequence of `Vertex` class instances.  Note that the list of
        vertices is assigned, *not copied* into the instance variable.

        Otherwise, `vertices` must be a sorted list of vertex valences
        and `edge_seq` a sequence of edge colorings, which are grouped
        according to the vertex valences to form the Graph vertices
        (by calling `vertex_factory` on each group of edges).

        This constructor is overloaded.  You can build a `Graph`
        instance either by specifying an explicit list of vertices::

          >>> G1 = Graph([Vertex([2,0,1]), Vertex([2,1,0])])

        Or by passing vertex valences and the flat list representation
        of the graph (i.e., concatenate all vertex representations in
        a single list of numbers)::

          >>> G2 = Graph((3,3), (2,0,1,2,1,0))

        Both forms yield the same result::
        
          >>> G2 == G1
          True
        
        """
        # build graph for explicit vertex list
        if edge_seq is None:
            assert debug.is_sequence_of_type(Vertex, vertices), \
                   "Graph.__init__: parameter `vertices` must be" \
                   " sequence of `Vertex` instances."
            #: list of vertices 
            self.vertices = vertices # FIXME: should be tuple?
            #: list of vertex valences 
            self._vertex_valences = tuple(len(v) for v in vertices) # FIXME: should be sorted?
            #: edge sequence from which this is/could be built
            self.edge_seq = tuple(chain(*[iter(v) for v in vertices]))

        # build graph from linear list and vertex valences
        else:
            assert debug.is_sequence_of_integers(vertices), \
                   "Graph.__init__: parameter `vertices` must be" \
                   " sequence of integers, but got '%s' instead" \
                   % vertices
            #: list of vertex valences 
            self._vertex_valences = tuple(vertices)
            assert (sum(vertices) % 2 ) == 0, \
                   "Graph.__init__: invalid parameter `vertices`:"\
                   "sum of vertex valences must be even."

            assert debug.is_sequence_of_integers(edge_seq), \
                   "Graph.__init__: parameter `edge_seq` must be sequence of integers, "\
                   "but got '%s' instead" % edge_seq

            #: edge sequence from which this graph is built
            self.edge_seq = tuple(edge_seq)
            
            #: edge sequence(s) identifying this graph; item at
            #  position 0 is `self.edge_seq`.
            self.edge_seq_aliases = set(self.edge_seq)

            #: list of vertices 
            self.vertices = []
            # Break up `edge_seq` into smaller sequences corresponding
            # to vertices.
            base = 0
            for current_vertex_index in xrange(len(vertices)):
                VLEN = vertices[current_vertex_index]
                self.vertices.append(vertex_factory(edge_seq[base:base+VLEN]))
                base += VLEN

        # init code common to both ctor variants:

        self.num_edges = sum(self._vertex_valences) / 2
        self.num_vertices = len(self.vertices)
        # these values will be computed on-demand
        self._boundary_components = None
        self._num_boundary_components = None
        self._genus = None
        self._valence_spectrum = None
        self._vertex_factory = vertex_factory
        
        #: edge sequence(s) identifying this graph; `self.edge_seq`
        #  always counts as an alias.
        self.edge_seq_aliases = set(self.edge_seq)
        
        # `self.endpoints` is the adjacency list of this graph.  For
        # each edge, store a pair `(v1, v2)` where `v1` and `v2` are
        # indices of endpoints.
        self.endpoints = [ [] for dummy in xrange(self.num_edges) ]
        for current_vertex_index in xrange(len(self.vertices)):
            for edge in self.vertices[current_vertex_index]:
                self.endpoints[edge].append(current_vertex_index)

    def __eq__(self, other):
        """Return `True` if Graphs `self` and `other` are isomorphic.

        Examples::

          >>> Graph([Vertex([1,0,0,1])]) == Graph([Vertex([1,1,0,0])])
          True
        """
        assert isinstance(other, Graph), \
               "Graph.__eq__:" \
               " called with non-Graph argument `other`: %s" % other
        # shortcuts
        if self._vertex_valences != other._vertex_valences:
            return False
        # if there is any representative edge sequence in common,
        # then these must be different presentations of the same graph
        if 0 != len(self.edge_seq_aliases.intersection(other.edge_seq_aliases)):
            return True

        # else, go the long way: try to find an explicit isomorphims
        # between graphs `self` and `other`
        morphisms = MorphismIteratorFactory(self.valence_spectrum())
        try:
            # if there is any morphism, then return `True`
            morphisms(self, other).next()
            return True
        except StopIteration:
            # list of morphisms is empty, graphs are not equal.
            return False

    def __getitem__(self, index):
        return self.vertices[index]

    def __hash__(self):
        return hash(self.edge_seq)

    def __iter__(self):
        """Return iterator over vertices."""
        return iter(self.vertices)

    def __repr__(self):
        return "Graph(%s)" % (repr(self.vertices))
    
    def __str__(self):
        return str(self.vertices)

    def add_alias(self, alias):
        """Add an edge sequence alias.
        """
        assert isinstance(alias, tuple), \
               "Graph.add_alias:" \
               " first argument `alias` must be tuple," \
               " but got `%s` instead." % alias
        assert sum(self._vertex_valences) == len(alias), \
               "Graph.add_alias:" \
               " `alias` length (%d) does not match the number of" \
               " edges in this graph (%d)." \
               % (sum(self._vertex_valences), len(alias))
        assert debug.is_sequence_of_integers(alias), \
               "Graph.add_alias:" \
               " first argument `alias` must be a sequence of integers," \
               " but got `%s` instead." % alias
        self.edge_seq_aliases.add(alias)

    def automorphisms(self):
        """Enumerate automorphisms of this `Graph` object.

        See `MorphismIteratorFactory` for details of how a `Graph`
        isomorphism is represented.
        """
        return MorphismIteratorFactory(self.valence_spectrum())(self, self)

    
    def boundary_components(self):
        """Return the number of boundary components of this `Graph` object.

        Each boundary component is represented by the list of (colored)
        edges::

          >>> Graph([Vertex([2,1,0]),Vertex([2,0,1])]).boundary_components()
          [[2, 0], [1, 2], [0, 1]]

        If both sides of an edge belong to the same boundary
        component, that edge appears twice in the list::

          >>> Graph([Vertex([2,1,1]),Vertex([2,0,0])]).boundary_components()
          [[2, 0, 2, 1], [1], [0]]
          
          >>> Graph([Vertex([2,1,0]),Vertex([2,1,0])]).boundary_components()
          [[2, 1, 0, 2, 1, 0]]
          
        """
        # try to return the cached value
        if self._boundary_components is not None:
            return self._boundary_components

        # otherwise, compute it now...
        
        # micro-optimizations
        L = self.num_edges
        ends = self.endpoints

        # pass1: build a "copy" of `graph`, replacing each edge
        # coloring with a triplet `(other, index, edge)` pointing to
        # the other endpoint of that same edge: the element at
        # position `index` in vertex `other`.
        pass1 = []
        for (vertex_index, vertex) in enumerate(self.vertices):
            replacement = []
            for (current_index, edge) in enumerate(vertex):
                (v1, v2) = ends[edge]
                if v1 != v2:
                    if ends[edge][0] == vertex_index:
                        other_end = ends[edge][1]
                    else:
                        other_end = ends[edge][0]
                    other_index = self.vertices[other_end].index(edge)
                else:
                    other_end = v1 # == v2, that is *this* vertex
                    # presume `current_index` is *not* the first
                    # occurrence of edge
                    other_index = vertex.index(edge)
                    if other_index == current_index:
                        # indeed it is, take next occurrence
                        other_index = vertex.index(edge, current_index+1)
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

        # save result for later reference
        self._boundary_components = result
        
        # that's all, folks!
        return result

    def contract(self, edgeno):
        """Return new `Graph` obtained by contracting the specified edge.

        The returned graph is presented in canonical form, that is,
        vertices are ordered lexicographically, longest ones first.
        
        Examples::

          >>> Graph([Vertex([2,2,0]), Vertex([0,1,1])]).contract(0)
          Graph([Vertex([1, 1, 0, 0])])
          >>> Graph([Vertex([2,1,0]), Vertex([2,0,1])]).contract(1)
          Graph([Vertex([1, 0, 0, 1])])

        The M_{1,1} trivalent graph yield the same result no matter
        what edge is contracted::

          >>> Graph([Vertex([2,1,0]), Vertex([2,1,0])]).contract(0)
          Graph([Vertex([1, 0, 1, 0])])
          >>> Graph([Vertex([2,1,0]), Vertex([2,1,0])]).contract(1)
          Graph([Vertex([1, 0, 1, 0])])
          >>> Graph([Vertex([2,1,0]), Vertex([2,1,0])]).contract(2)
          Graph([Vertex([1, 0, 1, 0])])
        """
        # check that we are not contracting a loop
        assert self.endpoints[edgeno][0] != self.endpoints[edgeno][1], \
               "Graph.contract: cannot contract a loop."
        # store position of the edge to be contracted at the endpoints
        i1 = self.endpoints[edgeno][0]
        i2 = self.endpoints[edgeno][1]
        pos1 = self.vertices[i1].index(edgeno)
        pos2 = self.vertices[i2].index(edgeno)

        # build new list of vertices, removing the contracted edge and
        # shifting all indices above
        subst = { edgeno:None } # delete specified edge
        for i in xrange(0, edgeno):
            subst[i] = i        # edges with lower color index are unchanged
        for i in xrange(edgeno+1, self.num_edges+1):
            subst[i] = i-1       # edges with higher color index are shifted down
        new_vertices = [ self._vertex_factory(itranslate(subst, v))
                         for v in self.vertices ]

        # Mate endpoints of contracted edge:
        # 0. make copies of vertices `v1`, `v2` so that subsequent
        #    operations do not alter the (possibly) shared `Vertex`
        #    object.
        v1 = copy(new_vertices[i1])
        v2 = copy(new_vertices[i2])
        # 1. Rotate endpoints `v1`, `v2` so that the given edge
        #    appears *last* in `v1` and *first* in `v2` (*Note:*
        #    the contracted edge has already been deleted from
        #    `v1` and `v2`, so index positions need to be adjusted):
        if (0 < pos1) and (pos1 < len(v1)):
            v1.rotate(pos1)
        if (0 < pos2) and (pos2 < len(v2)):
            v2.rotate(pos2)
        # 2. Join vertices by concatenating the list of incident
        # edges:
        v1.extend(v2)

        # set new `v1` vertex in place of old first endpoint, 
        new_vertices[i1]= v1
        # and remove second endpoint from list of new vertices
        del new_vertices[i2]

        # build new graph in canonical form
        return Graph(sorted(v.make_canonical() for v in new_vertices),
                     vertex_factory=self._vertex_factory)
        
    def edges(self):
        """Iterate over edge colorings."""
        return xrange(0, self.num_edges)
    
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
        """
        for a in self.automorphisms():
            if self.is_orientation_reversing(a):
                return False
        # no orientation reversing automorphism found
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
          >>> Graph([4, 4], [3, 3, 0, 0, 2, 2, 1, 1]).is_connected()
          False
          >>> Graph([4, 4], [3, 1, 2, 0, 3, 0, 2, 1]).is_connected()
          True
        """
        endpoints = self.endpoints
        visited_edges = set()
        visited_vertices = set()
        vertices_to_visit = [0]
        for vi in vertices_to_visit:
            # enqueue neighboring vertices that are not connected by
            # an already-visited edge
            for l in self.vertices[vi]:
                if l not in visited_edges:
                    # add other endpoint of this edge to the to-visit list
                    if endpoints[l][0] == vi:
                        other = endpoints[l][1]
                    else:
                        other = endpoints[l][0]
                    if other not in visited_vertices:
                        vertices_to_visit.append(other)
                    visited_edges.add(l)
                visited_vertices.add(vi)
        return (len(visited_vertices) == len(self.vertices))

    def is_loop(self, edge):
        """Return `True` if `edge` is a loop (i.e., the two endpoint coincide).
        """
        return self.endpoints[edge][0] == self.endpoints[edge][1]
        

    def is_orientation_reversing(self, automorphism):
        """Return `True` if `automorphism` reverses orientation of
        this `Graph` instance."""
        pv = automorphism[0]
        result = pv.sign()
        for (e0, e1) in [ (pv[e[0]], pv[e[1]])
                          for e in self.endpoints ]:
            if e0 > e1:
               result = -result
##         def is_increasing(a,b):
##             if a <= b:
##                 return +1
##             else:
##                 return -1
##         for x in xrange(self.num_edges):
##             result *= is_increasing(* self.endpoints[x]) \
##                       * is_increasing(pv[self.endpoints[x][0]],
##                                       pv[self.endpoints[x][1]])
        return (-1 == result)

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
                if self._valence_spectrum.has_key(l):
                    self._valence_spectrum[l].append(index)
                else:
                    self._valence_spectrum[l] = [index]
        return self._valence_spectrum
        

class MorphismIteratorFactory(object):
    """Make iterators over isomorphisms of graphs with a given valence
    spectrum.

    An isomorphism is represented by a tuple `(pv, rot, pe)` where:

      - `pv` is a permutation of ther vertices: the `i`-th vertex
        of `g1` is sent to the `pv[i]`-th vertex of `g2`, rotated
        by `rot[i]` places leftwards;

      - `pe` is a permutation of the edge colors: edge `i` in `g1`
        is mapped to edge `pe[i]` in `g2`.

    Examples::

      >>> g = Graph([Vertex([2, 1, 1]), Vertex([2, 0, 0])])
      >>> morphisms=MorphismIteratorFactory(g.valence_spectrum())
      >>> for f in morphisms(g, Graph([Vertex([2, 2, 0]), Vertex([1, 1, 0])])):
      ...   print f
      ({0: 1, 1: 0}, (1, 1), {0: 2, 1: 1, 2: 0})
      ({0: 0, 1: 1}, (1, 1), {0: 1, 1: 2, 2: 0})


      >>> for f in morphisms(g, g):
      ...   print f
      ({0: 1, 1: 0}, (0, 0), {0: 1, 1: 0, 2: 2})
      ({0: 0, 1: 1}, (0, 0), {0: 0, 1: 1, 2: 2})
      
"""

    __slots__ = (
        '_candidate_pvs',
        )

    def __init__(self, valence_spectrum):
        """Constructor, taking graph valence spectrum.
        """
        ## Compute all permutations of vertices that preserve valence.
        ## (A permutation `p` preserves vertex valence if vertices `v`
        ## and `p[v]` have the same valence.)

        # save valences as we have no guarantees that keys() method
        # will always return them in the same order
        valences = valence_spectrum.keys()

        # it's easier to compute vertex-preserving permutations
        # starting from the order vertices are given in the valence
        # spectrum; will rearrange them later.
        domain = Permutation(concat([ valence_spectrum[v] for v in valences ]))
        permutations_of_vertices_of_same_valence = [
            [ copy(p) for p in InplacePermutationIterator(valence_spectrum[v]) ]
            for v in valences
            ]

        #: Permutations of the vertex order that preserve valence.
        self._candidate_pvs = [
            domain.rearrange(concat(ps))
            for ps
            in SetProductIterator(permutations_of_vertices_of_same_valence)
            ]

    def __call__(self, g1, g2):
        """Return iterator over all isomorphisms from `g1` to `g2`."""
        # use "repetition patterns" to avoid mapping loop-free
        # vertices into vertices with loops, and viceversa.
        # FIXME: as loop-free vertices are much more common than
        # looped ones, this might turn out to be slower than just
        # checking all possible rotations.  Maybe just check
        # the number of loops instead of building the full repetition pattern?
        rps1 = [ v.repetition_pattern() for v in g1.vertices ]
        rps2 = [ v.repetition_pattern() for v in g2.vertices ]
        for vertex_index_map in self._candidate_pvs:
            pvrots = [ [] for x in xrange(len(g1.vertices)) ]
            for (j, (b1, rp1), (b2, rp2)) \
                    in izip(count(),
                            rps1,
                            (rps2[i] for i in vertex_index_map)):
                # Items in pvrots are lists (of length
                # `g1.num_vertices`), composed of tuples
                # `(i1,b1+s,i2,b2,s)`, meaning that vertex at index
                # `i1` in `g1` should be mapped to vertex at index
                # `i2` in `g2` with shift `s` and bases `b1` and `b2`
                # (that is, `g1.vertices[i1][b1+s:b1+s+len]` should be
                # mapped linearly onto `g2.vertices[i2][b2:b2+len]`).
                # `rp1` is a kind of "derivative" of `v1`; we gather
                # the displacement `b1+s` for `v1` by summing elements
                # of rp1 up to -but not including- `rp_shift`.
                pvrots[j].extend([ (j,b1+sum(rp1[:s]),vertex_index_map[j],b2)
                                   for s
                                   in rp1.all_shifts_for_linear_eq(rp2) ])
            for pvrot in SetProductIterator(pvrots):
                pe = Permutation()
                pe_is_ok = True  # optimistic default
                for (i1,b1,i2,b2) in pvrot:
                    v1 = g1.vertices[i1]
                    v2 = g2.vertices[i2]
                    if not pe.extend(v1[b1:b1+len(v1)],
                                     v2[b2:b2+len(v2)]):
                        # cannot extend, proceed to next `pvrot`
                        pe_is_ok = False
                        break
                if pe_is_ok and (len(pe) > 0):
                    pv = Permutation(t[2] for t in pvrot)
                    # Check that the combined action of `pv` and `pe`
                    # preserves the adjacency relation.  Note:
                    #   - we make list comprehensions of both adjacency lists
                    #     to avoid inverting `pe`: that is, we compare the
                    #     the adjacency lists in the order they have in `g1`,
                    #     but with the vertex numbering from `g2`;
                    #   - elements of the adjacency lists are made into
                    #     `set`s for unordered comparison;
                    #   - *copy* the objects to pass through `pv.translate`,
                    #     as `pv.translate` does in-place modify of its argument.
                    if 0 != cmp([ set(g2.endpoints[pe[x]])
                                  for x in xrange(g2.num_edges) ],
                                [ set(pv.translate(copy(g1.endpoints[x])))
                                  for x in xrange(g1.num_edges) ]):
                        # continue with next `pvrot`
                        continue
                    rots = tuple(t[1]-t[3] for t in pvrot)
                    yield (pv, rots, pe)


class ConnectedGraphsIterator(object):
    """Iterate over all connected graphs having vertices of the
    prescribed valences.
    
    Examples::

      >>> list(ConnectedGraphsIterator([4]))
      [Graph([Vertex([1, 0, 1, 0])]),
       Graph([Vertex([1, 1, 0, 0])])]

      >>> list(ConnectedGraphsIterator([3,3]))
      [Graph([Vertex([2, 0, 1]), Vertex([2, 0, 1])]),
       Graph([Vertex([2, 1, 0]), Vertex([2, 0, 1])]),
       Graph([Vertex([2, 1, 1]), Vertex([2, 0, 0])])]
    

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

    __slots__ = (
        'graphs',
        '_edge_seq_iterator',
        '_morphism_factory',
        '_vertex_factory',
        '_vertex_valences',
        )

    def __init__(self, vertex_valences, vertex_factory=VertexCache()):
        assert debug.is_sequence_of_integers(vertex_valences), \
               "ConnectedGraphsIterator: parameter `vertex_valences` must be a sequence of integers, "\
               "but got %s" % vertex_valences
        assert 0 == sum(vertex_valences) % 2, \
               "ConnectedGraphsIterator: sum of vertex valences must be divisible by 2"

        self._morphism_factory = None
        self._vertex_factory = vertex_factory
        self._vertex_valences = vertex_valences
        self.graphs = []

        # build list [0,0,1,1,...,n-1,n-1]
        starting_edge_seq=[]
        for l in xrange(0, sum(vertex_valences)/2):
            starting_edge_seq += [l,l]
        self._edge_seq_iterator = InplacePermutationIterator(starting_edge_seq)

    def __iter__(self):
        return self
    
    def next(self):
        for edge_seq in self._edge_seq_iterator:
            current = Graph(self._vertex_valences,
                            edge_seq,
                            self._vertex_factory)
            if not (current.is_canonical() and current.is_connected()):
                continue

            # the valence spectrum is the same for all graphs in the list,
            # so only compute it once
            if self._morphism_factory is None:
                self._morphism_factory = MorphismIteratorFactory(current.valence_spectrum())

            # now walk down the list and remove isomorphs
            current_is_not_isomorphic_to_already_found = True
            for candidate in self.graphs:
                # if there is any isomorphism, then reject current
                try:
                    self._morphism_factory(candidate, current).next()
                    # if we get here, an isomorphism has been found,
                    # so try again with a new `current` graph
                    current_is_not_isomorphic_to_already_found = False
                    break
                except StopIteration:
                    # no isomorphism has been found, try with next
                    # `candidate` graph
                    pass
            if current_is_not_isomorphic_to_already_found:
                # add current to graph list
                self.graphs.append(current)
                return current
            else:
                # record `current` as alias of `candidate` and proceed to next graph
                candidate.add_alias(current.edge_seq)

        # no more graphs to generate
        raise StopIteration



## main: run tests

if "__main__" == __name__:
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
