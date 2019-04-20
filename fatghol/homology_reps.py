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
import numpy as np
import scipy as sp
from scipy.io import mmread
from scipy.sparse.linalg import svds
from scipy.linalg import null_space
import os
from fatghol.graph_homology import (
        FatgraphComplex,
        MgnChainComplex,
        NumberedFatgraph
        )
from fatghol.combinatorics import (
        Permutation,
        PartitionIterator
)


# Convert SimpleMatrix to Numpy matrix. Returns a sparse scipy matrix
def simple_matrix_convert(mat):
    mat.save("./.temp.mat")
    M = mmread("./.temp.mat")
    os.remove(".temp.mat")
    return M


def permute_marked_fatgraph(fg, perm):
    for bc in fg.boundary_cycles:
        if fg.numbering[bc] in perm:
            fg.numbering[bc] = perm[fg.numbering[bc]]


def cycle_type_to_perm(cycle_type):
    d = {}
    cur = 0
    for k in cycle_type:
        for i in range(0, k):
            d[cur + i] = cur + (i + 1) % k
        cur += k
    return Permutation(d)


class NullSpaceComplex:
    def __init__(self, g, n):
        self.complex = FatgraphComplex(g, n)
        self.null_spaces = self._compute_null_spaces()
        self.n = n
        self.ci_characters = []

    def _compute_null_spaces(self):
        bnds = self.compute_boundary_operators()
        bases = []
        for d in range(len(bnds)):
            M = simple_matrix_convert(bnds[d][0])
            if(M.shape[0] == 0):
                continue
            bases.append((d, null_space(M.todense())))
        return bases

    def compute_ci_characters(self):
        characters = [{} for _ in xrange(len(self))]
        partitions = list(PartitionIterator(self.n, self.n))
        for partition in partitions:
            for i in xrange(len(self)):
                m = self.module[i]
                perm = cycle_type_to_perm(partition)

                characters[i][partition] = 0

                for fg in m:
                    g = NumberedFatgraph(fg.underlying, fg.numbering.copy())
                    permute_marked_fatgraph(g, perm)

                    isoms = list(NumberedFatgraph.isomorphisms(g, fg))
                    if (len(isoms) > 0):
                        characters[i][partition] += isoms[0].compare_orientations() * perm.sign()
        self.ci_characters = characters

    def _permute_vector(self, degree, vector, perm):
        assert len(vector) == len(self.module[degree]), \
            "Vector has smaller length than number of basis elements"
        m = self.module[degree]
        permuted_vector = [0 for _ in range(len(vector))]
        index = 0
        for pool in m.iterblocks():
            for k in xrange(len(pool)):
                (j, a) = pool._index(pool.numberings[k])
                permuted_vector[j] = vector[index] * a.compare_orientations() * perm.sign()
                index += 1
        return permuted_vector

    def compute_boundary_operators(self):
        return self.complex.compute_boundary_operators()


def kernel_comp_tests():
    for g in [0,1,2]:
        for n in [2,3,4,5]:
            if (g,n) in [(0,2),(1,4),(2,2),(1,5),(2,3),(2,4),(2,5),(0,5)]:
                continue
            print("g:", g, ", n:", n)
            C = FatgraphComplex(g, n)


if __name__=="__main__":
    kernel_comp_tests()
