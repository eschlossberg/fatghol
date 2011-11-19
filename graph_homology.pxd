from fractions import Fraction

from homology import ChainComplex, DifferentialComplex
from rg import Fatgraph


cdef class MgnChainComplex(ChainComplex):
    cdef readonly Fraction orbifold_euler_characteristics
    cpdef DifferentialComplex compute_boundary_operators(MgnChainComplex self)


cdef class NumberedFatgraph(Fatgraph):
    cpdef readonly Fatgraph underlying
    cpdef readonly dict numbering
    cdef NumberedFatgraph contract(self, int edgeno)
    cpdef isomorphisms(NumberedFatgraph self, NumberedFatgraph other)


cdef class NumberedFatgraphPool(object):
    cdef readonly Fatgraph graph
    cdef readonly bint is_orientable
    cdef readonly list P
    cdef readonly list P_
    cdef readonly list numberings
    cdef readonly int num_automorphisms
    cdef list facets(NumberedFatgraphPool self, int edge, NumberedFatgraphPool other)
    cdef _push_fwd(f, g1, g2)
    cdef tuple _compute_nb_map(pb_map, NumberedFatgraphPool src, NumberedFatgraphPool dst)
    cdef _index(self, list numbering)


cpdef MgnChainComplex FatgraphComplex(int g, int n)


