#! /usr/bin/env python
#
"""Functions related to valences of ribbon graphs.
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


#import cython


from fatghol.combinatorics import PartitionIterator


#@cython.ccall
#@cython.locals(g=cython.int, n=cython.int,
#               L=cython.int, V=cython.int, result=set)
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
