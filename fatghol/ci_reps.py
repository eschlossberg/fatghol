#! /usr/bin/env python
#
"""Classes for computing the representations of S_n on M_g,n
"""
#
#   Copyright (C) 2019-2023 Eli Schlossberg <eschlossb@gmail.com>
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

from fatghol.graph_homology import (
    NumberedFatgraph,
)

from combinatorics import (
    Permutation,
    PartitionIterator
)


# Permute the markings of the fatgraph fg by the permutation perm
def permute_marked_fatgraph(perm, fg):
    for bc in fg.boundary_cycles:
        if fg.numbering[bc] in perm:
            fg.numbering[bc] = perm[fg.numbering[bc]]


# Create a permutation of the given cycle type
def cycle_type_to_perm(cycle_type):
    d = {}
    cur = 0
    for k in cycle_type:
        for i in range(0, k):
            d[i] = cur + (i + 1) % k
        cur += k
    return Permutation(d)


# Compute the character of a given degree for M_g,n
def ci_char(ci, n):
    partitions = list(PartitionIterator(n, n))
    char = {}
    for partition in partitions:
        # add an entry to the character table for each partition
        perm = cycle_type_to_perm(partition)
        char[partition] = 0

        # iterate over the basis fatgraphs
        for i in range(len(ci)):
            fg = ci[i]

            # create a copy of the fatgraph and permute its numbering
            g = NumberedFatgraph(fg.underlying, fg.numbering.copy())
            permute_marked_fatgraph(perm, g)

            # check if it is mapped to itself
            isoms = list(NumberedFatgraph.isomorphisms(g, fg))
            if len(isoms) > 0:
                print(len(isoms))
                char[partition] += isoms[0].compare_orientations()*perm.sign()
