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
    FatgraphComplex
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
            d[cur + i] = cur + (i + 1) % k
        cur += k
    return Permutation(d)


# Compute the sign of a permutation given its cycle type
def sign(perm):
    par = 1
    for i in perm:
        if (i % 2 == 0):
            par *= -1
    return par


# Compute the character of a given degree for M_g,n
class MgnChars:
    def __init__(self, g, n):
        self.g = g
        self.n = n
        self.mgn = FatgraphComplex(g, n)
        self.chars = []

        self.ci_permutations = self._compute_permutations(n)
        
        for i in range(self.mgn.length):
            self.chars.append(self._ci_char(i))


    # Computes the action of S_n on the basis of each C_i
    # Returns a list of the format i->perm_cycle_type->fg->\sigma(fg)
    # with i being the degree of C_i, perm_cycle_type the cycle type of the
    # desired permutation, fg the fatgraph being acted on, and \sigma(fg)
    # the result of the action
    def _compute_permutations(self, n):
        partitions = list(PartitionIterator(n, n))
        
        # Generate the permutations of all cycle types
        perms = []
        for partition in partitions:
            perms.append((partition, cycle_type_to_perm(partition)))

        ci_perms = []
        for i in range(self.mgn.length):
            ci_fg_perms = {}
            ci = self.mgn.module[i]
            
            for perm in perms:
                fg_perms = {}
                j = 0
                for fg in ci:
                    g = NumberedFatgraph(fg.underlying, fg.numbering.copy())
                    permute_marked_fatgraph(perm[1], g)
                    fg_perms[j] = g
                    j += 1

                ci_fg_perms[perm[0]] = fg_perms
            ci_perms.append(ci_fg_perms)
        return ci_perms


    def _ci_char(self, i):
        ci = self.mgn.module[i]
        perms = self.get_action(i)
        chars = {}
        for perm in perms:
            chars[perm] = 0
            sigma = perms[perm]
            j = 0
            for fg in ci:
                # find the image of fg under sigma
                g = sigma[j]

                # check if it is mapped to itself
                isoms = list(NumberedFatgraph.isomorphisms(g, fg))
                if len(isoms) > 0:
                    chars[perm] += isoms[0].compare_orientations()*sign(perm)
                j += 1
        return chars

    # Returns the action of S_n on degree i
    def get_action(self, i):
        return self.ci_permutations[i]


# Verify the computations of the characters of S_n on M_g,n for small g,n
def ci_char_test():
    for g in [0,1,2]:
        for n in [2,3,4,5]:
            if (g,n) in [(0,2),(1,4),(2,2),(1,5),(2,3),(2,4),(2,5)]:
                continue
            print("g:", g, "n:", n)
            mgn = MgnChars(g, n)
            chars = mgn.chars
            for char in chars:
                print(char)

if __name__=="__main__":
    ci_char_test()
