#! /usr/bin/env python

import operator
import types


def tr(s, t1, t2):
    """Change every occurrence (in sequence `s`) of an element of set `t1` with the corrisponding element of set `t2`.

    Examples:
    >>> tr([0,1,0,0],[0],[2])
    [2, 1, 2, 2]
    """
    def _tr(elt, t1, t2):
        """If `elt` equals some element in set `t1`, then return the corresponding element from set `t2`, otherwise return `elt` unchanged.
        """
        try:
            return t2[t1.index(elt)]
        except ValueError:
            return elt
    return [_tr(x, t1, t2) for x in s]


def dispo(l, n=None):
    """Return all possible dispositions of `n` symbols in a row of length `l`.

    If `n` is not specified, it is taken to be equal to `l`.

    Examples:
    >>> list(dispo(1,5))
    [[1], [2], [3], [4], [5]]
    >>> list(dispo(2))
    [[1, 2], [2, 1]]
    """
    if n is None:
        n = l
    for s0 in range(1,n+1):
        if 1 == l:
            yield [s0]
        else:
            for s1 in dispo(l-1,n-1):
                yield [s0]+tr(s1, [s0], [n])


class CyclicList(list):
    """List with indices wrapping around.

    Examples:
    >>> a=CyclicList([1,2,3])
    >>> len(a)
    3
    >>> a[3]
    1
    >>> a[4]
    2
    >>> a[4]=5
    >>> a
    [1, 5, 3]
    >>> a[1:2]
    [5]
    >>> a[1:2]=[7,8,9]
    >>> a
    [1, 7, 8, 9, 3]
    >>> len(a)
    5
    >>> a=CyclicList([1,2,3])
    >>> b=CyclicList([2,3,1])
    >>> b == a
    True
    >>> c=CyclicList([3,1,2])
    >>> c==a
    True
    >>> d=CyclicList([1,3,2])
    >>> d==a
    False
    """
    def __init__(self, sequence=None):
        if sequence is None:
            list.__init__(self)
        else:
            list.__init__(self, sequence)
    def __getitem__(self, i): return list.__getitem__(self, i % len(self))
    def __setitem__(self, i, item): list.__setitem__(self, i % len(self), item)
    def __delitem__(self, i): list.__delitem__(self, i%len(self))
    def __getslice__(self, i, j):
        """*Note:* returned slice is of `list` type!"""
        i = max(i, 0); j = max(j, 0); l = len(self)
        if abs(i-j) < l:
            return self.__class__(list.__getslice__(self, i%l, j%l))
        else:
            a = max(i,j)
            b = min(i,j)
            c = l*(a/l)     # nearest multiple of l below a
            d = l*(b/l + 1) # nearest multiple of l above b
            n = (a/l) - (b/l + 1) # how many times the whole seq is repeated in the middle
            return list.__getslice__(self, b%l,l) \
                   + (n*list.__getslice__(self, 0,l)) \
                    + list.__getslice__(self, 0,a%l)
    def __setslice__(self, i, j, other):
        i = max(i, 0); j = max(j, 0); l = len(self)
        if isinstance(other, list):
            list.__setslice__(self, i%len(self), j%len(self), other)
        else:
            list.__setslice__(self, i%len(self), j%len(self), list(other))
    def __delslice__(self, i, j):
        i = max(i, 0); j = max(j, 0)
        list.__delslice__(self, i%len(self), j%len(self))
    def shift_for_list_eq(self, other, start=0):
        """Return minimum shift index `b >= start` such that `self[b:b+len]==other` as Python lists.

        Examples:
        >>> a=CyclicList([1,2,3])
        >>> b=CyclicList([2,3,1])
        >>> a.shift_for_list_eq(b)
        1
        >>> a.shift_for_list_eq(b,2)
        None
        """
        l = len(self)
        def _eq_shifted(first, second, shift):
            l=len(first)
            i=0
            while i < l:
                if first[i+shift] != second[i]:
                    return False
                else:
                    i += 1
            return True
        shift = start
        while shift < l:
            if _eq_shifted(self, other, shift):
                return shift
            else:
                shift += 1
        return None
    def all_shifts_for_list_eq(self, other):
        start = 0
        l = len(self)
        while start < l:
            shift = self.shift_for_list_eq(other, start)
            if shift is None:
                break
            else:
                yield shift
            start = shift+1
    def __eq__(self, other):
        """Compare `self` with all possible translations of `other`."""
        if len(other) != len(self):
            return False
        elif None == self.shift_for_list_eq(other):
            return False
        else:
            return True

def repetition_pattern(c):
    """Return the repetition pattern of a cyclic list `c`.

    The repetition pattern is a *cyclic* list of integers: each number
    `n` in the repetition list corresponds to `n` equal-valued items
    in the list `c`.

    Examples:
    >>> c=CyclicList([1,2,3])
    >>> repetition_pattern(c)
    (1, [1, 1, 1])
    >>> c=CyclicList([4,4,4])
    >>> repetition_pattern(c)
    (0, [3])
    >>> c=CyclicList([4,4,4,1])
    >>> repetition_pattern(c)
    (3, [1, 3])
    >>> c=CyclicList([1,4,4,4,1])
    >>> repetition_pattern(c)
    (1, [3, 2])
    >>> c=CyclicList([4,4,4,3,3])
    >>> repetition_pattern(c)
    (3, [2, 3])
    """
    l=len(c)
    # find the first transition
    b=0
    while (b < l) and (c[b] == c[b+1]):
        b += 1
    # all items are equal
    if b == l:
        return (0, CyclicList([l]))
    # else, start building pattern from here
    result=CyclicList()
    b += 1
    i=0
    while i < l:
        j=0
        while (j < l) and (c[b+i+j] == c[b+i+j+1]):
            j += 1
        result.append(j+1)
        i += j+1
    return (b, result)

class RotatedList(CyclicList):
    """Access a given `CyclicList` with indices displaced by a fixed amount.

    Examples:
    >>> a = [1, 7, 8, 9, 3]
    >>> b=RotatedList(1,a)
    >>> len(b)
    5
    >>> b[0]
    7
    >>> b[5]
    7
    >>> b[1]
    8
    >>> b[0]=6
    >>> b[0]
    6
    """
    def __init__(self, displacement, initial=None):
        self.shift = displacement
        CyclicList.__init__(self, initial)
    def __getitem__(self, i): return CyclicList.__getitem__(self, i + self.shift)
    def __setitem__(self, i, item): CyclicList.__setitem__(self, i + self.shift, item)
    def __delitem__(self, i): CyclicList.__delitem__(self, i + self.shift)
    def __getslice__(self, i, j):
        return CyclicList.__getslice__(self, i+self.shift, j+self.shift)
    def __setslice__(self, i, j, other):
        CyclicList.__setslice__(self, i+self.shift, j+self.shift, other)
    def __delslice__(self, i, j):
        CyclicList.__delslice__(self, i+self.shift, j+self.shift)


class Map:
    def __init__(self, order):
        self.order = order
        self.map = {}
    def __call__(self, src):
        return self.map[src]
    def extend(self, srcs, dsts):
        """Return `True` if the Map can be extended by mapping
        elements of `srcs` to corresponding elements of `dsts`.
        Return `False` if any of the new mappings conflicts
        with an already established mapping.
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
    

def all_edge_seq(n):
    """Return all lists representing edges of a ribbon graph.

    Each returned list has length `2*n` and comprises the symbols `{1, ..., n}`,
    each of which is repeated exactly twice.
    """
    for s in dispo(2*n):
        s = tr(s, range(n+1,2*n+1), range(1,n+1))
        yield s


def vertices_from_edge_seq(edge_seq, vertex_list):
    """Break up `edge_seq` into smaller sequences corresponding to vertices.

    Second argument `vertex_list` must be a list of vertex valences.
    """
    result = []
    vl = vertex_list[:]
    while len(vl)>0:
        v = vl.pop(0)
        result.append(edge_seq[:v])
        del edge_seq[:v]
    return result


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


def is_canonical(graph):
    """Return `True` if `graph` is canonical.
    A graph is canonical iff:
    1) Each vertex is represented by the maximal sequence, among all
       sequences representing the same cyclic order.
    2) Vertices are sorted in lexicographic order.

    Examples:
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
    for vertex in graph:
        if not is_maximal_representative(vertex):
            return False
        if previous_vertex and (previous_vertex < vertex):
            return False
        previous_vertex = vertex
    return True


def all_graphs(vertex_list):
    """Return all graphs having vertices of the given valences.
    """
    total_edges = sum(vertex_list) / 2
    # gather all distinct edge sequences
    graphs = []
    for edge_seq in all_edge_seq(total_edges):
        g = vertices_from_edge_seq(edge_seq, vertex_list)
        if is_canonical(g):
            graphs.append(g)
    # filter out sequences representing isomorphic graphs
    pos = 0
    while pos < len(graphs):
        current = graphs[pos]
        # slight optimization: since `current` is constant in the loop below,
        # pre-compute as much as we can...
        current_cy = [CyclicList(v) for v in current]
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
                              in map(None, current_cy, current_b_rp, candidate)]:
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
                # list will shift up because of the deletion in the
                # middle.
                del graphs[pos2]
            else:
                # advance to next candidate
                pos2 += 1
        pos += 1
    return graphs


def other(pair, one):
    """Return the member of `pair` not equal to `one`."""
    if pair[0] == one:
        return pair[1]
    else:
        return pair[0]


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
    def enumerate_certesian_product(p):
        """Iterate over all elements in the cartesian products of elements of items in `p`.

        Examples:
        >>> list(enumerate_cartesian_product([[1],[1]])
        [[1,1]]
        >>> list(enumerate_cartesian_product([[1,2],[1]])
        [[1,1],[1,2)]
        >>> list(enumerate_cartesian_product([[1,2],[1,2]])
        [[1,1],[1,2],[2,1],[2,2]]
        """
        if len(p) == 0:
            yield []
        else:
            for i in p[-1]:
                for js in enumerate_cartesian_product(p[:-1]):
                    yield js+[i]

    def deep_cmp(s1,s2):
        """Compare items in `s1` and `s2`, recursing into subsequences.

        Examples:
        >>> deep_cmp(1,1)
        0
        >>> deep_cmp([1],[1])
        0
        >>> deep_cmp([1,1],[1,1])
        0
        >>> deep_cmp([1,[1]],[1,[1]])
        0
        >>> deep_cmp([1,[1]],[1,[2]])
        -1
        """
        if not (type(s1) == type(s2)):
            raise TypeError, \
                "Comparing arguments of different type: %s vs %s" \
                % (repr(type(s1)), repr(type(s2)))
        else:
            try:
                # assume s1,s2 are sequences and recursively apply this
                # function to pairs of corresponding elements...
                def _first_nonzero(x,y):
                    if 0 != x:
                        return x
                    else:
                        return y
                return reduce(_first_nonzero, map(deep_cmp, s1, s2), 0)
            except TypeError:
                # ...if s1,s2 are not sequences, then do a builtin comparison
                return cmp(s1,s2)

    for a in enumerate_cartesian_product(candidates):
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


def graph_to_xypic(graph):
    r"""Print XY-Pic code snippet to render graph `graph`.

    Examples:
    #>>> print graph_to_xypic([[3,2,1],[3,1,2]])
    \xy 0;<2cm,0cm>:%
    (1,1)="v1",%
    (-1,1)="v2",%
    "v1",{\xypolygon4"v1l"{~:{(1.20,0):(0,-1)::}~={90}~>{}}},%
    "v2",{\xypolygon4"v2l"{~:{(1.20,0):}~={270}~>{}}},%
    "v1"*\txt{[321]},%
    "v2"*\txt{[312]},%
    "v1";"v2"**\crv{"v1l3"&"v2l3"},%
    "v1";"v2"**\crv{"v1l2"&"v2l1"},%
    "v1";"v2"**\crv{"v1l1"&"v2l2"},%
    \endxy
    #>>> print graph_to_xypic([[1,2,3,1],[1,2,3,2],[1,2,3,3]])
    \xy 0;<2cm,0cm>:%
    {\xypolygon3"v"{~={0}~>{}}},% mark v1,v2,v3
    "v1",{\xypolygon8"v1l"{~:{(1.20,0):}~={90}~>{}}},%
    "v2",{\xypolygon8"v2l"{~:{(1.20,0):}~={210}~>{}}},%
    "v3",{\xypolygon8"v3l"{~:{(1.20,0):}~={330}~>{}}},%
    "v1"*\txt{[1231]},
    "v2"*\txt{[1232]},
    "v3"*\txt{[1233]},
    "v1";"v2"**\crv{"v1l2"&"v2l1"},%
    "v1";"v3"**\crv{"v1l3"&"v3l1"},%
    "v2";"v3"**\crv{"v2l3"&"v3l2"},%
    "v1";"v1"**\crv{"v1l1"&"v1l4"},%
    "v2";"v2"**\crv{"v2l2"&"v2l4"},%
    "v3";"v3"**\crv{"v3l3"&"v3l4"},%
    \endxy
    """
    def vertex_label(v):
        return '['+("".join(map(str,v)))+']'
    label = map(vertex_label, graph)
    K = len(graph) # number of vertices
    result = r'\xy 0;<2cm,0cm>:%'+'\n'
    # put graph vertices on a regular polygon
    if K < 3:
        result += '(1,1)="v1",%\n(-1,1)="v2",%\n'
    else:
        result += r'{\xypolygon%d"v"{~={0}~>{}}},%% mark vertices' \
                  % (max(K,3)) \
                  + '\n'
    # mark invisible "control points" for bezier curves connecting vertices
    def rotation_angle(K,k):
        return 90+(k-1)*360/K
    for k in range(1,K+1):
        # "vK",{\xypolygonL"vKl"{~{1.20,0):~={...}~>{}}},%
        result += r'"v%d",{\xypolygon%d"v%dl"{~:{(1.20,0):q}~={%d}~>{}}},%%' \
                  % (k, 2*len(graph[k-1])-2, k, rotation_angle(K,k)) \
                  + '\n'
    for k in range(1,K+1):
        # "vK"*\txt{[...]},%
        result += r'"v%d"*\txt{%s},%%' \
                  % (k, label[k-1]) \
                  + '\n'
    for l in range(1,num_edges(graph)+1):
        (v1,v2) = endpoints(l,graph)
        if v1 != v2:
            result += r'"v%d";"v%d"**\crv{"v%dl%d"&"v%dl%d"},%%' \
                      % (v1+1, v2+1, v1+1, 1+graph[v1].index(l),
                         v2+1, 1+graph[v2].index(l)) \
                      + '\n'
        else:
            h = graph[v1].index(l)
            result += r'"v%d";"v%d"**\crv{"v%dl%d"&"v%dl%d"},%%' \
                      % (v1+1, v2+1, v1+1, h+1, v2+1, 1+graph[v1].index(l,h+1)) \
                      + '\n'
    (g,n) = classify(graph)
    result += r'0*\txt{g=%d,n=%d}' % (g,n) + '\n'
    result += r'\endxy'
    return result


## main

if "__main__" == __name__:
    # parse command-line options
    from optparse import OptionParser
    parser = OptionParser(usage="""Usage: %prog [options] action [arg ...]

    Actions:

      mgn G N
        Print the vertex valences occurring in M_{g,n} graphs

      vertices V1,V2,...
        Print the graphs having only vertices of the specified valences.

      test
        Run internal code tests and report results.
        
    """)
    parser.add_option("-v", "--verbose",
                      action='store_true', dest='verbose', default=False,
                      help="Report verbosely on progress.")
    parser.add_option("-L", "--latex",
                      action='store_true', dest='latex', default=False,
                      help="Print Xy-Pic code to draw graphs.")
    parser.add_option("-o", "--output", dest="outfile", default=None,
                      help="Output file for `vertices` action.")
    (options, args) = parser.parse_args()

    if 0 == len(args):
        parser.print_help()
        
    elif 'test' == args[0]:
        import doctest
        doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    elif 'mgn' == args[0]:
        if len(args) < 3:
            parser.print_help()
        g = int(args[1])
        n = int(args[2])
        # with 1 vertex only, there are this many edges:
        L = 2*g + n - 1
        K = 1
        def part(N,K):
            if K == 1:
                if N >= 3:
                    return [[N]]
                else:
                    return []
            result = []
            for k in range(3,N-2):
                for p in part(N-k,K-1):
                    result.append([k]+p)
            return result
        while True:
            p = part(2*L, K)
            if len(p) == 0:
                break
            for pi in p:
                print pi
            K += 1
            L = 2*g + n + K - 2

    elif 'vertices' == args[0]:
        del args[0]
        import sys
        # open output file
        if options.outfile is None:
            outfile = sys.stdout
        else:
            outfile = open(options.outfile, 'w')
        # compute graphs matching given vertex sequences
        graphs = []
        for vertex_pattern in args:
            vertex_list = map(int, vertex_pattern.split(','))
            if options.verbose:
                sys.stderr.write("Computing graphs with vertex pattern %s\n" \
                                 % vertex_list)
            graphs += all_graphs(vertex_list)

        # print latex code
        if options.latex:
            outfile.write(r"""
\documentclass[a4paper,twocolumn]{article}
\usepackage[curve,poly,xdvi]{xy}
\begin{document}
""")
        for g in graphs:
            if options.latex:
                outfile.write(graph_to_xypic(g)+'\n')
            else:
                outfile.write("%s\n" % g)
        outfile.write("\n")
        outfile.write("Found %d graphs.\n" % len(graphs))
        outfile.write("\n")
        if options.latex:
            outfile.write(r"\end{document}")
            outfile.write("\n")
