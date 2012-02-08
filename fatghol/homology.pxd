
cdef class VectorSpace(object):
    cdef readonly list base
    cdef readonly int dimension
    cpdef dict coordinates(self, list combo)

cdef class DifferentialComplex(list):
    # FIXME: `bd` should be of `SimpleMatrix` type
    cpdef append(self, bd, int ddim, int cdim)
    cpdef compute_homology_ranks(self)

cdef class ChainComplex(object):
    cdef readonly int length
    cdef list differential
    cdef list module
    cpdef DifferentialComplex compute_boundary_operators(self)
    cpdef list compute_homology_ranks(self)
     
