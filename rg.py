#! /usr/bin/env python
#
"""
"""
__docformat__ = 'reStructuredText'


from utils import *
from cyclic import CyclicList

from itertools import *
import operator


def vertex_valences_for_given_g_and_n(g,n):
    """
    Examples::
      >>> vertex_valences_for_given_g_and_n(0,3)
      set([(3, 3), (4,)])
    """
    # with 1 vertex only, there are this many edges:
    L = 2*g + n - 1
    K = 1
    result = set()
    while True:
        ps = partitions(2*L, K)
        if len(ps) == 0:
            break
        for p in ps:
            p.sort()
            result.add(tuple(p))
        K += 1
        L = 2*g + n + K - 2
    return result


class Vertex(CyclicList):
    """A (representative of) a vertex of a ribbon graph.

    A vertex is represented by the cyclically ordered list of its
    (decorated) edges.  The edge colorings may be accessed through a
    (read-only) sequence interface.
    """
    __slots__ = ('_repetition_pattern',)
    def __init__(self, edge_seq, start=0, end=None):
        """Create `Vertex` instance by excerpting the slice `[start:end]` in `edge_seq`.
        """
        if end is None:
            end = len(edge_seq)
        CyclicList.__init__(self, edge_seq[start:end])

        # the following values will be computed when they are first requested
        self._repetition_pattern = None
##     def __new__(cls, edge_seq, start=0, end=None):
##         """Create `Vertex` instance by excerpting the slice `[start:end]` in `edge_seq`.
##         """
##         if end is None:
##             end = len(edge_seq)
##         return CyclicArray.__new__(cls, edge_seq[start:end])
##     def __init__(self, *args, **kwargs):
##         # the following values will be computed when they are first requested
##         self._repetition_pattern = None
       
    def __iter__(self):
        """Return iterator over edges."""
        return list.__iter__(self)

    #def __repr__(self):
    #    return str(self._edges)
    #def __str__(self):
    #    return str(self._edges)
    
    def is_maximal_representative(self):
        """Return `True` if this `Vertex` object is maximal among
        representatives of same cyclic sequence.
        
        Examples::
          >>> Vertex([3,2,1]).is_maximal_representative()
          True
          >>> Vertex([2,1,3]).is_maximal_representative()
          False
          >>> Vertex([1,1]).is_maximal_representative()
          True
          >>> Vertex([1]).is_maximal_representative()
          True
        """
        L = len(self)
        for i in xrange(1,L):
            for j in xrange(0,L):
                # k := (i+j) mod L
                k = i+j
                if k >= L:
                    k %= L
                if self[k] < self[j]:
                    # continue with next i
                    break
                elif self[k] > self[j]:
                    return False
                # else, continue comparing
        return True

    def repetition_pattern(self):
        """Return the repetition pattern of this `Vertex` object.
        Same as calling `repetition_pattern(this._edges)` but caches
        results.
        """
        if self._repetition_pattern is None:
            self._repetition_pattern = CyclicList.repetition_pattern(self) 
        return self._repetition_pattern


class Graph(object):
    """A fully-decorated ribbon graph.

    Exports a (read-only) sequence interface, through which vertices
    can be accessed.
    """
    def __init__(self, vertex_valences, edge_seq):
        assert is_sequence_of_integers(vertex_valences), \
               "Graph.__init__: parameter `vertex_valences` must be sequence of integers, "\
               "but got '%s' instead" % vertex_valences
        self._vertex_valences = vertex_valences
        assert (sum(vertex_valences) % 2 ) == 0, \
               "Graph.__init__: invalid parameter `vertex_valences`:"\
               "sum of vertex valences must be even."

        self._num_edges = sum(self._vertex_valences) / 2
        self._num_vertices = len(self._vertex_valences)
        # these values will be computed on-demand
        self._num_boundary_components = None
        self._genus = None
        self._valence_spectrum = None
        
        assert is_sequence_of_integers(edge_seq), \
               "Graph.__init__: parameter `edge_seq` must be sequence of integers, "\
               "but got '%s' instead" % edge_seq
        assert max(edge_seq) == self._num_edges - 1, \
               "Graph.__init__: invalid parameter `edge_seq`:"\
               "Sequence of edges %s doesn't match number of edges %d" \
               % (edge_seq, self._num_edges)
        # Break up `edge_seq` into smaller sequences corresponding to vertices.
        self.vertices = []
        base = 0
        for current_vertex_index in xrange(len(vertex_valences)):
            VLEN = vertex_valences[current_vertex_index]
            # FIXME: this results in `edge_seq` being copied into smaller
            # subsequences; can we avoid this by defining a list-like object
            # "vertex" as a "view" on a portion of an existing list?
            self.vertices.append(Vertex(edge_seq, base, base+VLEN))
            base += VLEN

    def __getitem__(self, index):
        return self.vertices[index]

    def __iter__(self):
        """Return iterator over vertices."""
        return iter(self.vertices)

    def __repr__(self):
        return "Graph(%s, %s)" \
               % (repr(self._vertex_valences), repr(self.vertices))
    
    def __str__(self):
        return str(self.vertices)

    def automorphisms(self):
        """Enumerate automorphisms of this `Graph` object.

        An automorhism is represented as a pair of ordered lists `(dests,
        rots)`: the i-th vertex of `graph` is to be mapped to the vertex
        `dests[i]`, and rotated by `rots[i]`.
        """
        # build enpoints vector for the final check that a constructed map
        # is an automorphism
        ev = [ self.endpoints(l) for l in self.edges() ]
        ev.sort()

        # gather valences and repetition pattern at
        # start for speedup
        valence = [ len(vertex) for vertex in self.vertices ]
        rp = [ vertex.repetition_pattern() for vertex in self.vertices ]

        ## pass 1: for each vertex, list all destinations it could be
        ## mapped to, in the form (dest. vertex, rotation).

        # pre-allocate list of right size; all elements must be empty
        # lists, that we fill with `.append()` later on
        candidates = [ [] ] * len(graph)
        
        # FIXME: if vertex `i` can be mapped into vertex `j`, with some
        # rotation delta, then vertex `j` can be mapped into vertex `i`
        # with rotation `-delta`, so rearrange this to only do
        # computations for `i>j` and use the values already available in
        # the other case...
        num_vertices = self.num_vertices()
        for i in xrange(num_vertices):
            for j in xrange(num_vertices):
                # if valences don't match, skip to next vertex in list
                if valence[i] != valence[j]:
                    continue
                # if repetition patterns don't match, skip to next vertex in list
                if not (rp[i] == rp[j]):
                    continue
                # append `(destination vertex, rotation shift)` to
                # candidate destinations list
                for delta in rp[i].all_shifts_for_list_eq(rp[j]):
                    candidates[i].append((j,delta))

        ## pass 2: for each vertex, pick a destination and return the resulting
        ##         automorphism. (FIXME: do we need to check that the adjacency 
        ##         matrix stays the same?)
        for a in enumerate_set_product(candidates):
            # check that map does not map two distinct vertices to the same one
            already_assigned = []
            a_is_no_real_map = False
            for dest in a:
                v = dest[0]
                if v in already_assigned:
                    a_is_no_real_map = True
                    break
                else:
                    already_assigned.append(v)
            if a_is_no_real_map:
                # try with next map `a`
                continue
            # check that the endpoints vector stays the same
            vertex_permutation = [ elt[0] for elt in a ]
            new_ev = [
                (vertex_permutation[e[0]], vertex_permutation[e[1]])
                for e in ev
                ]
            new_ev.sort()
            if 0 != deep_cmp(ev, new_ev):
                # this is no automorphism, skip to next one
                continue
            # return automorphism in (vertex_perm_list, rot_list) form
            yield ([ elt[0] for elt in a ],
                   [ elt[1] for elt in a ])

    def classify(self):
        """Return the pair (g,n) for this `Graph` object."""
        return (self.genus(), self.num_boundary_components())

    def edges(self):
        """Iterate over edge colorings."""
        return xrange(0, self.num_edges())
    
    def endpoints(self, n):
        """Return the endpoints of edge `n`.
    
        The endpoints are returned as a pair (v1,v2) where `v1` and `v2`
        are indices of endpoint vertices in this `Graph` object.
        """
        result = []
        for vi in xrange(len(self.vertices)):
            c = self.vertices[vi].count(n)
            if 2 == c:
                return (vi,vi)
            elif 1 == c:
                result.append(vi)
        if 0 == len(result):
            raise KeyError, "Edge %d not found in graph '%s'" % (n, repr(self))
        return result
    
    def genus(self):
        """Return the genus g of this `Graph` object."""
        # compute value if not already done
        if (self._genus is None):
            n = self.num_boundary_components()
            K = self.num_vertices()
            L = self.num_edges()
            # by Euler, K-L+n=2-2*g
            self._genus = (L - K - n + 2) / 2
        return self._genus

    def has_orientation_reversing_automorphism(self):
        """Return `True` if `Graph` has an orientation-reversing automorphism.

        Enumerate all automorphisms of `graph`, end exits with `True`
        result as soon as one orientation-reversing one is found.
        """
        for a in self.automorphisms():
            if is_orientation_reversing(a):
                return True
        return False

    def is_orientation_reversing(self, automorphism):
        """Return `True` if `automorphism` reverses orientation of this `Graph` instance."""
        def sign_of_rotation(l,r=1):
            """Return sign of a rotation of `l` elements, applied `r` times."""
            # evaluating a conditional is faster than computing (-1)**...
            if 0 == ((l-1)*r) % 2:
                return 1
            else:
                return -1
        def sign_of_permutation(p):
            """Return sign of permutation `p`.

            A permutation is represented as a linear list: `p` maps `i` to
            `p[i]`.  Items of `p` are required to be valid indices in `p`,
            that is, `p[i] =< max(p)` for all `i`.

            This is an adaptation of the `perm_sign` code by John Burkardt
            (see it among the collection at
            http://orion.math.iastate.edu/burkardt/f_src/subset/subset.f90
            or http://www.scs.fsu.edu/~burkardt/math2071/perm_sign.m ); it
            computes the sign by counting the number of interchanges
            required to change the given permutation into the identity
            one.

            Examples::
              >>> sign_of_permutation([1,2,3])
              1
              >>> sign_of_permutation([1,3,2])
              -1
              >>> sign_of_permutation([3,1,2])
              1
              >>> sign_of_permutation([1])
              1
              >>> sign_of_permutation([])
              1
            """
            n = len(p)
            s = +1
            # get elements back in their home positions
            for j in xrange(n):
                q = p[j]
                if q !=j :
                    p[j],p[q] = p[q],q # interchange p[j] and p[p[j]]
                    s = -s             # and account for the interchange
            # note that q is now in its home position
            # whether or not an interchange was required
            return s
        return reduce(operator.mul,
                      [ sign_of_rotation(l,r)
                        for l,r in
                        ([ len(v) for v in self.vertices],
                         automorphism[1]) ],
                      sign_of_permutation(automorphism[0]))

    def is_canonical(self):
        """Return `True` if this `Graph` object is canonical.

        A graph is canonical iff:
        1) Each vertex is represented by the maximal sequence, among all
           sequences representing the same cyclic order.
        2) Vertices are sorted in lexicographic order.

        Examples::
          >>> Graph([3,3],[2,1,0,2,1,0]).is_canonical()
          True             
          >>> Graph([3,3],[2,1,0,2,0,1]).is_canonical()
          True             
          >>> Graph([3,3],[2,0,1,2,1,0]).is_canonical()
          False
          >>> Graph([3,3],[0,1,2,2,1,0]).is_canonical()
          False 
        """
        previous_vertex = None
        for vertex in self.vertices:
            if not vertex.is_maximal_representative():
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
          >>> Graph([3,3], [2,1,0,2,1,0]).num_boundary_components()
          1
          >>> Graph([3,3], [2,1,0,2,0,1]).num_boundary_components()
          3
        """
        # try to return the cached value
        if self._num_boundary_components is not None:
            return self._num_boundary_components

        # otherwise, compute it now...
        
        L = self.num_edges()
        # for efficiency, gather all endpoints now
        ends = [ self.endpoints(l) for l in xrange(L) ]

        # pass1: build a "copy" of `graph`, replacing each edge coloring
        # with a pair `(other, index)` pointing to the other endpoint of
        # that same edge: the element at position `index` in vertex
        # `other`.
        pass1 = []
        def other(pair, one):
            """Return the member of `pair` not equal to `one`."""
            if pair[0] == one:
                return pair[1]
            else:
                return pair[0]
        for (vertex_index, vertex) in enumerate(self.vertices):
            replacement = []
            for (current_index, edge) in enumerate(vertex):
                (v1, v2) = ends[edge]
                if v1 != v2:
                    other_end = other(ends[edge], vertex_index)
                    other_index = self.vertices[other_end].index(edge)
                else:
                    other_end = v1 # == v2, that is *this* vertex
                    try:
                        # presume this is the first occurrence of edge...
                        other_index = vertex.index(edge, current_index+1)
                    except ValueError:
                        # it's not, take first
                        other_index = vertex.index(edge)
                # replace other_index with index of *next* edge
                # (in the vertex cyclic order)
                if other_index == len(self.vertices[other_end])-1:
                    other_index = 0
                else:
                    other_index += 1
                replacement.append((other_end, other_index))
            pass1.append(replacement)

        # pass2: now build a linear list, each element of the list
        # corresponding to an edge, of `(pos, seen)` where `pos` is the
        # index in this list where the other endpoint of that edge is
        # located, and `seen` is a flag indicating whether this side of
        # the edge has already been walked through.
        pass2 = []
        # build indices to the where each vertex begins in the linear list
        vi=[0]
        for vertex in self.vertices:
            vi.append(vi[-1]+len(vertex))
        # build list from collapsing the 2-level structure
        for vertex in pass1:
            for pair in vertex:
                pass2.append([vi[pair[0]]+pair[1],False])

        # pass3: pick up each element of the linear list, and follow it
        # until we come to an already marked one.
        result = 0
        pos = 0
        while pos < len(pass2):
            # fast forward to an element that we've not yet seen
            while (pos < len(pass2)) and (pass2[pos][1] == True):
                pos += 1
            if pos >= len(pass2):
                break
            # walk whole chain of edges
            i = pos
            while pass2[i][1] == False:
                pass2[i][1] = True
                i = pass2[i][0]
            result += 1
            pos += 1

        # save result for later reference
        self._num_boundary_components = result
        
        # that's all, folks!
        return result

    def num_edges(self):
        return self._num_edges

    def num_vertices(self):
        return self._num_vertices

    def valence_spectrum(self):
        """Return a dictionary mapping valences into vertex indices."""
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
        

def all_graphs(vertex_valences):
    """Return all graphs having vertices of the given valences.

    Examples::
      >>> all_graphs([4])
      [Graph([4], [[1, 0, 1, 0]]),
       Graph([4], [[1, 1, 0, 0]])]
      >>> all_graphs([3,3])
      [Graph([3, 3], [[2, 0, 1], [2, 0, 1]]),
       Graph([3, 3], [[2, 1, 0], [2, 0, 1]]),
       Graph([3, 3], [[2, 1, 1], [2, 0, 0]])]
    """
    assert is_sequence_of_integers(vertex_valences), \
           "all_graphs: parameter `vertex_valences` must be a sequence of integers, "\
           "but got %s" % vertex_valences

    total_edges = sum(vertex_valences) / 2

    ## pass 1: gather all canonical graphs built from edge sequences
    graphs = list(all_canonical_decorated_graphs(vertex_valences, total_edges))
    if len(graphs) == 0:
        return []

    ## pass 2: filter out sequences representing isomorphic graphs
    def vertex_permutations(valence_spectrum):
        """Return all permutations of vertices that preserve valence.
        (A permutation `p` preserves vertex valence if vertices `v`
        and `p[v]` have the same valence.)

        The passed parameter `valence_spectrum` is a dictionary,
        mapping vertex valence to the list of (indices of) vertices
        having that valence.

        Returns a pair `(list_of_mappings, domain)`.

        Examples::
          >>> vertex_permutation({ 3:[0,1] })
          ([[1, 0], [0, 1]],
           [0, 1])
          >>> vertex_permutation({ 3:[0], 5:[1] })
          ([[0, 1]],
           [0, 1])
        """
        # save valences as we have no guarantees that keys() method
        # will report them always in the same order
        valences = valence_spectrum.keys()
        vertex_to_vertex_mapping_domain = concat([ valence_spectrum[v] for v in valences ])
        permutations_of_vertices_of_same_valence = [
            [ p[:] for p in InplacePermutationIterator(valence_spectrum[v]) ]
            for v in valences
            ]
        vertex_to_vertex_mappings = [
            concat(ps)
            for ps in enumerate_set_product(permutations_of_vertices_of_same_valence)
            ]
        return (vertex_to_vertex_mappings, vertex_to_vertex_mapping_domain)
    # the valence spectrum is the same for all graphs in the list
    vertex_to_vertex_mappings, vertex_to_vertex_mapping_domain = \
                               vertex_permutations(graphs[0].valence_spectrum())
    
    pos = 0
    while pos < len(graphs):
        current = graphs[pos]
        pos2 = pos+1
        while pos2 < len(graphs):
            candidate = graphs[pos2]
            candidate_is_isomorphic_to_current = False
            if candidate == current:
                candidate_is_isomorphic_to_current = True
            else:
                for vertex_index_map in vertex_to_vertex_mappings:
                    morphism = Mapping(total_edges)
                    for i1,i2 in izip(vertex_to_vertex_mapping_domain, vertex_index_map):
                        v1 = current.vertices[i1]
                        (b1, rp1) = v1.repetition_pattern()
                        v2 = candidate.vertices[i2]
                        (b2, rp2) = v2.repetition_pattern()
                        rp_shift = rp1.shift_for_list_eq(rp2)
                        if rp_shift is None:
                            # cannot map vertices, quit looping on vertices
                            break
                        # rp1 is a kind of "derivative" of v1; we gather the 
                        #  displacement for v1 by summing elements of rp1 up to
                        # -but not including- `rp_shift`.
                        shift = sum(rp1[:rp_shift])
                        if not morphism.extend(v1[b1+shift:b1+shift+len(v1)],v2[b2:b2+len(v2)]):
                            # continue with next candidate
                            morphism = None
                            break
                    if morphism and morphism.is_complete():
                        # the two graphs are isomorphic
                        candidate_is_isomorphic_to_current = True
            if candidate_is_isomorphic_to_current:
                # delete candidate; do *not* advance `pos2`, as the
                # list would be shifted up because of the deletion in
                # the middle.
                del graphs[pos2]
            else:
                # advance to next candidate
                pos2 += 1
        pos += 1
    return graphs


class InplacePermutationIterator:
    """Iterate over all permutations of a given sequence.

    The given sequence `seq` is altered as new permutations are
    requested through the `next()` method.

    The code is a port of the C++ STL one, as described in:
      http://marknelson.us/2002/03/01/next-permutation

    Examples::
      >>> [ x[:] for x in InplacePermutationIterator([0]) ]
      [[0]]
      >>> [ x[:] for x in InplacePermutationIterator([0,1]) ]
      [[1, 0], [0, 1]]
      >>> [ x[:] for x in InplacePermutationIterator([0,1,2])]
      [[0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0], [0, 1, 2]]
    """
    def __init__(self, seq, start=0, end=None):
        self.seq = seq
        self.start = start
        if end is None:
            end = len(self.seq)
        self.end = end
        if (start == end):
            self.enumeration_finished = True
        else:
            self.enumeration_finished = False
    def __iter__(self):
        return self
    def next(self):
        """Return next permutation of initially given `sequence`."""
        if self.enumeration_finished:
            raise StopIteration
        i = self.end - 1
        while True:
            if (i == self.start):
                self.seq.reverse()
                self.enumeration_finished = True
                return self.seq
            j = i
            i -= 1
            if self.seq[i] < self.seq[j]:
                k = self.end-1
                while self.seq[i] >= self.seq[k]:
                    k -= 1
                # swap seq[i] and seq[k]
                self.seq[i],self.seq[k] = self.seq[k],self.seq[i]
                # reverse slice seq[j:] *in-place*
                if self.end-2 == j:
                    self.seq[j],self.seq[-1] = self.seq[-1],self.seq[j]
                else:
                    for l in xrange(0, (self.end - 1 - j) / 2):
                        self.seq[j+l],self.seq[-1-l] = self.seq[-1-l],self.seq[j+l]
                return self.seq


def all_edge_seq(n):
    """Iterate over lists representing edges of a ribbon graph.

    Each returned list has length `2*n` and comprises the symbols `{0,...,n-1}`,
    each of which is repeated exactly twice.
    """
    # build list [0,0,1,1,...,n-1,n-1]
    edges=[]
    for l in xrange(0,n):
        edges += [l,l]
    # iterate over all its permutations
    ps = InplacePermutationIterator(edges)
    for s in ps:
        yield s


def all_canonical_decorated_graphs(vertex_valences, total_edges):
    """Iterate over all canonical decorated graphs with `total_edges` edges."""
    for edge_seq in all_edge_seq(total_edges):
        g = Graph(vertex_valences, edge_seq)
        if g.is_canonical():
            yield g


class Mapping:
    """A mapping (in the mathematical sense) with a domain of fixed cardinality.

    Provides methods to incrementally construct the mapping by
    extending an existing map with new source->destination pairs.  The
    cardinality of the domain should be known in advance, and the
    `Mapping` object is *complete* when all items in the domain have been
    assigned a destination value.
    """
    def __init__(self, order):
        self.order = order
        self.map = {}
    def __call__(self, src):
        return self.map[src]
    def extend(self, srcs, dsts):
        """Return `True` if the mapping can be extended by mapping
        elements of `srcs` to corresponding elements of `dsts`.
        Return `False` if any of the new mappings conflicts with an
        already established one.
        """
        for i in range(0,len(srcs)):
            if self.map.has_key(srcs[i]):
                if self.map[srcs[i]] != dsts[i]:
                    return False
                else:
                    pass
            else:
                self.map[srcs[i]] = dsts[i]
        return True
    def extend_with_hash(self, mappings):
        """Return `True` if the mapping can be extended by mapping each
        key of `mappings` to the corresponding value.  Return `False`
        if any of the new mappings conflicts with an already
        established one.
        """
        for src in mappings.keys():
            if self.map.has_key(src):
                if self.map[src] != mappings[src]:
                    return False
                else:
                    pass
            else:
                self.map[src] = mappings[src]
        return True
    def is_complete(self):
        return self.order == len(self.map.keys())
    def is_assigned(self, src):
        return self.map.has_key(src)



## main: run tests

if "__main__" == __name__:
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
