#! /usr/bin/env python
#
"""
"""
__docformat__ = 'reStructuredText'


from utils import *
from cyclic import CyclicList,repetition_pattern

import operator


def vertex_valences_for_given_g_and_n(g,n):
    """
    Examples::
      >>> vertex_valences_for_given_g_and_n(0,3)
      [[4], [3, 3]]
    """
    # with 1 vertex only, there are this many edges:
    L = 2*g + n - 1
    K = 1
    result = []
    while True:
        ps = partitions(2*L, K)
        if len(ps) == 0:
            break
        for p in ps:
            result.append(p)
        K += 1
        L = 2*g + n + K - 2
    return result



class Graph:
    """A fully-decorated ribbon graph.
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

        assert is_sequence_of_integers(edge_seq), \
               "Graph.__init__: parameter `edge_seq` must be sequence of integers, "\
               "but got '%s' instead" % edge_seq
        assert max(edge_seq) == self._num_edges, \
               "Graph.__init__: invalid parameter `edge_seq`:"\
               "Sequence of edges %s doesn't match number of edges %d" \
               % (edge_seq, self._num_edges)
        self._edge_seq = edge_seq
        # Break up `edge_seq` into smaller sequences corresponding to vertices.
        self.vertices = []
        base = 0
        for current_vertex_index in xrange(len(vertex_valences)):
            VLEN = vertex_valences[current_vertex_index]
            # FIXME: this results in `edge_seq` being copied into smaller
            # subsequences; can we avoid this by defining a list-like object
            # "vertex" as a "view" on a portion of an existing list?
            self.vertices.append(edge_seq[base:base+VLEN])
            base += VLEN

    def __str__(self):
        return str(self.vertices)
            
    def num_edges(self):
        return self._num_edges

    def num_vertices(self):
        return len(self.vertices)

    def is_canonical(self):
        """Return `True` if this `Graph` object is canonical.

        A graph is canonical iff:
        1) Each vertex is represented by the maximal sequence, among all
           sequences representing the same cyclic order.
        2) Vertices are sorted in lexicographic order.

        Examples::
          >>> is_canonical([[3,2,1],[3,2,1]])
          True
          >>> is_canonical([[3,2,1],[3,1,2]])
          True
          >>> is_canonical([[3,1,2],[3,2,1]])
          False
          >>> is_canonical([[1,2,3],[3,2,1]])
          False
          >>> is_canonical([[1,2],[3,2,1]])
          False
        """
        previous_vertex = None
        for vertex in self.vertices:
            if not is_maximal_representative(vertex):
                return False
            if previous_vertex and (previous_vertex < vertex):
                return False
            previous_vertex = vertex
        return True
    

def all_graphs(vertex_valences):
    """Return all graphs having vertices of the given valences.

    Examples::
      >>> all_graphs([4])
      [[[2, 1, 2, 1]],
       [[2, 2, 1, 1]]]
      >>> all_graphs([3,3])
      [[[3, 1, 2], [3, 1, 2]],
       [[3, 2, 1], [3, 1, 2]],
       [[3, 2, 2], [3, 1, 1]]]
    """
    assert is_sequence_of_integers(vertex_valences), \
           "all_graphs: parameter `vertex_valences` must be a sequence of integers, "\
           "but got %s" % vertex_valences
    total_edges = sum(vertex_valences) / 2

    ## pass 1: gather all canonical graphs built from edge sequences
    graphs = list(all_canonical_decorated_graphs(vertex_valences, total_edges))

    ## pass 2: filter out sequences representing isomorphic graphs
    pos = 0
    while pos < len(graphs):
        current = graphs[pos]
        # slight optimization: since `current` is constant in the loop below,
        # pre-compute as much as we can...
        current_cy = [CyclicList(v) for v in current.vertices]
        current_b_rp = [repetition_pattern(v_cy) for v_cy in current_cy]

        pos2 = pos+1
        while pos2 < len(graphs):
            candidate = graphs[pos2]
            candidate_is_isomorhic_to_current = False
            if candidate == current:
                candidate_is_isomorhic_to_current = True
            else:
                perm = Map(total_edges)
                # FIXME: could save some processing time by caching
                # the cyclic list of vertices) for any `candidate`, at
                # the expense of memory usage...
                for v1,b1,rp1,v2 in [(w1,b1,rp1,CyclicList(w2)) for (w1,(b1,rp1),w2)
                              in map(None, current_cy, current_b_rp, candidate.vertices)]:
                    (b2, rp2) = repetition_pattern(v2)
                    rp_shift = rp1.shift_for_list_eq(rp2)
                    if rp_shift is None:
                        # cannot map vertices, quit looping on vertices
                        break
                    # rp1 is a kind of "derivative" of v1; we gather the
                    # displacement for v1 by summing elements of rp1 up to
                    # -but not including- `rp_shift`.
                    shift = sum(rp1[:rp_shift])
                    if not perm.extend(v1[b1+shift:b1+shift+len(v1)],v2[b2:b2+len(v2)]):
                        # continue with next candidate
                        perm = None
                        break
                if perm and perm.completed():
                    # the two graphs are isomorphic
                    candidate_is_isomorhic_to_current = True
            if candidate_is_isomorhic_to_current:
                # delete candidate; do *not* advance `pos2`, as the
                # list would be shifted up because of the deletion in
                # the middle.
                del graphs[pos2]
            else:
                # advance to next candidate
                pos2 += 1
        pos += 1
    return graphs


def all_edge_seq(n):
    """Iterate over lists representing edges of a ribbon graph.

    Each returned list has length `2*n` and comprises the symbols `{1,...,n}`,
    each of which is repeated exactly twice.
    """
    for s in permutations(2*n):
        tr_inplace(s, range(n+1,2*n+1), range(1,n+1))
        yield s


def all_canonical_decorated_graphs(vertex_valences, total_edges):
    """Iterate over all canonical decorated graphs with `total_edges` edges."""
    for edge_seq in all_edge_seq(total_edges):
        g = Graph(vertex_valences, edge_seq)
        if g.is_canonical():
            yield g
            

def is_maximal_representative(vertex):
    """Return `True` if `vertex` is maximal among representatives of same cyclic sequence.

    Examples:
    >>> is_maximal_representative([3,2,1])
    True
    >>> is_maximal_representative([2,1,3])
    False
    >>> is_maximal_representative([1,1])
    True
    >>> is_maximal_representative([1])
    True
    """
    def wrap_index(i,l):
        if i >= l:
            return i%l
        else:
            return i
    l = len(vertex)
    for i in range(1,l):
        for j in range(0,l):
            if vertex[wrap_index(i+j,l)] < vertex[j]:
                # continue with next i
                break
            elif vertex[wrap_index(i+j,l)] > vertex[j]:
                return False
            # else, continue comparing
    return True


def num_edges(g):
    """Return the total number of edges of graph `g`."""
    return max(map(max,g))


def endpoints(n,g):
    """Return the endpoints of edge `n` in graph `g`.

    The endpoints are returned as a pair (v1,v2) where `v1` and `v2`
    are indices of vertices in `g`.
    """
    result = []
    for v in range(0, len(g)):
        c = g[v].count(n)
        if 2 == c:
            return (v,v)
        elif 1 == c:
            result.append(v)
    if 0 == len(result):
        raise KeyError, "Edge %d not found in graph '%s'" % (n, repr(g))
    return result


def num_boundary_components(graph):
    """Return the number of boundary components of `graph`.

    Each boundary component is represented by the list of (colored)
    edges.

    Examples:
    >>> num_boundary_components([[3,2,1],[3,2,1]])
    1
    >>> num_boundary_components([[3,2,1],[3,1,2]])
    3
    """
    L = num_edges(graph)+1
    # for efficiency, gather all endpoints now
    ends = L*[None]
    for l in range(1,L):
        ends[l] = endpoints(l, graph)

    # pass1: build a "copy" of `graph`, replacing each edge coloring
    # with a pair `(other, index)` pointing to the other endpoint of
    # that same edge: the element at position `index` in vertex
    # `other`.
    pass1 = []
    for (vertex, vertex_index) in map(None, graph, range(0,len(graph))):
        replacement = []
        for (edge, current_index) in map(None, vertex, range(0,len(vertex))):
            (v1, v2) = ends[edge]
            if v1 != v2:
                other_end = other(ends[edge], vertex_index)
                other_index = graph[other_end].index(edge)
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
            if other_index == len(graph[other_end])-1:
                other_index = 0
            else:
                other_index += 1
            replacement.append((other_end, other_index))
        pass1.append(replacement)

    # pass2: now build a linear list, each element of the list
    # corresponding to an edge, of `(pos, seen)` where `pos` is the
    # index in this list where the other endpoint of that edge is
    # located, and `seen` is a flag indicating whether this side of
    # the edge has already been walked by.
    pass2 = []
    # build indices to the where each vertex begins in the linear list
    vi=[0]
    for vertex in graph:
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

    # that's all, folks!
    return result


def classify(graph):
    """Return the pair (g,n) for `graph`."""
    n = num_boundary_components(graph)
    K = len(graph)
    L = num_edges(graph)
    # by Euler, K-L+n=2-2*g
    g = (L - K - n + 2) / 2
    return (g,n)


def all_automorphisms(graph):
    """Enumerate all automorphisms of `graph`.

    An automorhism is represented as a pair of ordered lists `(dests,
    rots)`: the i-th vertex of `graph` is to be mapped to the vertex
    `dests[i]`, and rotated by `rots[i]`.
    """
    # build enpoints vector for the final check that a constructed map
    # is an automorphism
    ev = []
    for l in range(1,num_edges(graph)+1):
        ev.append(endpoints(l, graph))
    ev.sort()
    
    # gather valences and repetition pattern at
    # start for speedup
    valence = [len(vertex) for vertex in graph]
    rp = [repetition_pattern(vertex) for vertex in graph]

    ## pass 1: for each vertex, list all destinations it could be mapped to,
    ##         in the form (dest. vertex, rotation).
    candidates = [ [] ] * len(graph)
    # FIXME: if vertex i can be mapped into vertex j, with some rotation delta,
    # then vertex j can be mapped into vertex i with rotation -delta,
    # so rearrange this to only do computations for i>j and use the values already
    # available in the other case...
    for i in range(len(graph)):
        for j in range(len(graph)):
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
        def first(seq):
            return seq[0]
        for dest in a:
            v = first(dest)
            if v in already_assigned:
                a_is_no_real_map = True
            else:
                already_assigned.append(v)
        if a_is_no_real_map:
            # try with next map
            continue
        # check that the endpoints vector stays the same
        vertex_permutation = map(first, a)
        new_ev = [(vertex_permutation[e[0]], vertex_permutation[e[1]]) for e in ev]
        new_ev.sort()
        if 0 != deep_cmp(ev, new_ev):
            # this is no automorphism, skip to next one
            continue
        # return automorphism in (vertex_perm_list, rot_list) form
        def second(seq):
            return seq[1]
        yield (map(first, a), map(second, a))


def is_orientation_reversing(graph, automorphism):
    """Return `True` if `automorphism` reverses orientation of `graph`."""
    def sign_of_rotation(l,r=1):
        """Return the sign of a rotation of `l` elements, applied `r` times."""
        # evaluating a conditional is faster than computing (-1)**...
        if 0 == ((l-1)*r) % 2:
            return 1
        else:
            return -1
    def sign_of_permutation(p):
        """Recursively compute the sign of permutation `p`.

        A permutation is represented as a linear list: `p` maps
        `i` to `p[i]`.

        Examples:
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
        l = len(p)
        if 1 >= l:
            return 1
        # find highest-numbered element
        k = p.index(l)
        # remove highest-numbered element for recursion
        q = p[:]
        del q[k]
        # recursively compute
        if 0 == ((l+k) % 2):
            s = -1
        else:
            s = 1
        return s * sign_of_permutation(q)
    return reduce(operator.mul,
                  map(sign_of_rotation,
                      map(len, graph), automorphism[1]),
                  sign_of_permutation(automorphism[0]))
                                         

def has_orientation_reversing_automorphism(graph):
    """Return `True` if `graph` has an orientation-reversing automorphism.
    
    Enumerate all automorphisms of `graph`, end exits with `True`
    result as soon as one orientation-reversing one is found.
    """
    for a in all_automorphisms(graph):
        if is_orientation_reversing(a):
            return True
    return False


class Map:
    """A mapping (in the mathematical sense) with a domain of fixed cardinality.

    Provides methods to incrementally construct the mapping by
    extending an existing map with new source->destination pairs.  The
    cardinality of the domain should be known in advance, and the
    `Map` object is *complete* when all items in the domain have been
    assigned a destination value.
    """
    def __init__(self, order):
        self.order = order
        self.map = {}
    def __call__(self, src):
        return self.map[src]
    def extend(self, srcs, dsts):
        """Return `True` if the Map can be extended by mapping
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
        """Return `True` if the Map can be extended by mapping each
        key of `mappings` tothe corresponding value.  Return `False`
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
    def completed(self):
        return self.order == len(self.map.keys())
    def is_assigned(self, src):
        return self.map.has_key(src)



## main: run tests

if "__main__" == __name__:
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
