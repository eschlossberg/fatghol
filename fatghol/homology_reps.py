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
        )


# Convert SimpleMatrix to Numpy matrix. Returns a sparse scipy matrix
def simple_matrix_convert(mat):
    mat.save("./.temp.mat")
    M = mmread("./.temp.mat")
    os.remove(".temp.mat")
    return M


class NullSpaceComplex:
    def __init__(self, FgComplex):
        self.complex = FgComplex
        self.null_spaces = self._compute_null_spaces()

    def _compute_null_spaces(self):
        bnds = self.complex.compute_boundary_operators()
        bases = []
        for d in range(len(bnds)):
            M = simple_matrix_convert(bnds[d][0])
            if(M.shape[0] == 0):
                continue
            bases.append((d, null_space(M.todense())))
        return bases


def kernel_comp_tests():
    for g in [0,1,2]:
        for n in [2,3,4,5]:
            if (g,n) in [(0,2),(1,4),(2,2),(1,5),(2,3),(2,4),(2,5),(0,5)]:
                continue
            print("g:", g, ", n:", n)
            C = FatgraphComplex(g, n)
            print(compute_null_spaces(C))


if __name__=="__main__":
    kernel_comp_tests()
