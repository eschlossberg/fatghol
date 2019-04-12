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
import sympy
from scipy.io import mmread
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


# Compute the kernel of the boundary map
def boundary_map_kernel(boundary_map):
    return sp.linalg.null_space(M)


def kernel_comp_tests():
    for g in [0,1,2]:
        for n in [2,3,4,5]:
            if (g,n) in [(0,2),(1,4),(2,2),(1,5),(2,3),(2,4),(2,5)]:
                continue
            print("g:", g, ", n:", n)
            bases = []
            C = FatgraphComplex(g, n)
            bnds = C.compute_boundary_operators()
            for d in range(len(bnds)):
                M = simple_matrix_convert(bnds[d][0]).todense()
                if(M.shape[0] == 0):
                    continue
                M = sympy.Matrix(M)
                bases.append((d, np.matrix(M.nullspace())))
            print(bases)

if __name__=="__main__":
    kernel_comp_tests()
