#! /usr/bin/env python
#
"""Functions related to valences of ribbon graphs.
"""
__docformat__ = 'reStructuredText'


from combinatorics import PartitionIterator


def vertex_valences_for_given_g_and_n(g,n):
    """Return all possible valences of fatgraphs appearing in the
    `(g,n)` complex.

    The returned partitions are sorted in ascending order.

    Examples::
      >>> vertex_valences_for_given_g_and_n(0,3)
      [(3, 3), (4,)]
      >>> vertex_valences_for_given_g_and_n(1,1)
      [(3, 3), (4,)]
    """
    # with 1 vertex only, there are this many edges:
    L = 2*g + n - 1
    V = 1
    result = set(PartitionIterator(2*L, V, 3))
    while (2*L >= 3*V):
        V += 1
        L = 2*g + n + V - 2 # by Euler's formula
        ## Each edge has *two* vertex endpoints (possibly the same vertex),
        ## so the set of vertex valences is a partition of `2*L`, where `L`
        ## is the number of edges, with some further constraints:
        ##   - each vertex has valence at least 3;
        ##   - thus, no vertex can have valence greater than `2*L-3*(V-1)`,
        ##     that is, supposing all vertices except one are 3-valent,
        ##     the highest-valent vertex must have valence `2L - 3(V-1)`.
        result.update(set(PartitionIterator(2*L, V, 3, 2*L-3*(V-1))))
    return sorted(result)





## main: run tests

if "__main__" == __name__:
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
