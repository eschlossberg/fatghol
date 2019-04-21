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
        self.boundary_operators = []
        self.null_spaces = self._compute_null_spaces()
        self.n = n
        self.ci_perms = self._compute_permutations(self.n)
        self.ci_characters = self._compute_ci_characters()
        self.null_space_characters = []
        self.compute_null_space_characters()
        self.homology_characters = []
        self.compute_homology_characters()

    def __len__(self):
        return len(self.complex)

    def _compute_permutations(self, n):
        partitions = list(PartitionIterator(n, n))

        # Generate the permutations of all cycle types
        perms = []
        for partition in partitions:
            perms.append((partition, cycle_type_to_perm(partition)))

        ci_perms = []
        for i in range(self.complex.length):
            ci_fg_perms = {}
            ci = self.complex.module[i]

            for perm in perms:
                fg_perms = {}
                for j in range(len(ci)):
                    fg = ci[j]
                    g = NumberedFatgraph(fg.underlying, fg.numbering.copy())
                    permute_marked_fatgraph(g, perm[1])

                    isom = None
                    for k in range(len(ci)):
                        isoms = list(NumberedFatgraph.isomorphisms(g, ci[k]))
                        if len(isoms) > 0:
                            index = k
                            isom = isoms[0]
                            break
                    assert isom is not None
                    sign = isom.compare_orientations()

                    fg_perms[j] = (k, sign * perm[1].sign())

                ci_fg_perms[perm[0]] = fg_perms
            ci_perms.append(ci_fg_perms)
        return ci_perms

    def _compute_null_spaces(self):
        bnds = self.compute_boundary_operators()
        bases = []
        for d in range(len(bnds)):
            M = simple_matrix_convert(bnds[d][0])
            self.boundary_operators.append(M)
            if(M.shape[0] == 0):
                bases.append((d, []))
                continue
            bases.append((d, null_space(M.todense())))
        return bases

    def _compute_ci_characters(self):
        characters = [{} for _ in xrange(len(self))]
        for i in xrange(len(self)):
            for partition in self.ci_perms[i]:
                characters[i][partition] = 0
                for j in range(len(self.complex.module[i])):
                    if j == self.ci_perms[i][partition][j][0]:
                        characters[i][partition] += 1* self.ci_perms[i][partition][j][1]
        return characters

    def _permute_vector(self, degree, vector, perm):
        assert vector.shape[0] == len(self.complex.module[degree]), \
            "Vector has smaller length than number of basis elements"
        m = self.complex.module[degree]
        permuted_vector = [0 for _ in range(len(vector))]
        for i in range(len(vector)):
            (j, sign) = self._permute_basis_vector(degree, i, perm)
            permuted_vector[j] = vector[i] * sign
        return permuted_vector

    # Returns the tuple (index', sign) of \sigma(C_i[index]) where \sigma is the
    # representative of its cycle type
    def _permute_basis_vector(self, degree, index, cycle_type):
        return self.ci_perms[degree][cycle_type][index]

    def compute_boundary_operators(self):
        return self.complex.compute_boundary_operators()

    # Compute the characters of the null spaces of the boundary map of degree i
    def null_space_character(self, degree):
        basis = self.null_spaces[degree][1]
        char = {}
        partitions = list(PartitionIterator(self.n, self.n))

        for perm in partitions:
            char[perm] = 0
            if(type(basis) == list):
                continue
            for i in range(basis.shape[1]):
                x = basis[:,i]
                x = np.array(x)
                x_prime = self._permute_vector(degree, x, perm)
                char[perm] += np.dot(x, x_prime)
            char[perm] = int(round(char[perm]))
        return char

    # Compute all the null space characters
    def compute_null_space_characters(self):
        for i in range(len(self)):
            self.null_space_characters.append(self.null_space_character(i))

    def compute_homology_character(self, degree):
        assert degree >= 0
        chi = self.null_space_characters[degree]
        if degree == 0:
            return chi

        for cycle_type in chi:
            chi[cycle_type] -= self.ci_characters[degree-1][cycle_type]\
                    - self.null_space_characters[degree-1][cycle_type]
        return chi

    def compute_homology_characters(self):
        for i in range(len(self)):
            self.homology_characters.append(self.compute_homology_character(i))




def kernel_comp_tests():
    for g in [0,1,2]:
        for n in [2,3,4,5]:
            if (g,n) in [(0,2),(1,4),(2,2),(1,5),(2,3),(2,4),(2,5),(0,5)]:
                continue
            print("g:", g, ", n:", n)
            C = FatgraphComplex(g, n)


if __name__=="__main__":
    kernel_comp_tests()
