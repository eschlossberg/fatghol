from native cimport *

from cache cimport Cacheable
from cycliclist cimport CyclicList
from iterators cimport BufferingIterator

cdef class Vertex(CyclicList):
    cdef readonly int num_loops

cdef class EqualIfIsomorphic(Cacheable):
    cdef readonly invariants
    cdef dict _cache_eq

cdef class Fatgraph(EqualIfIsomorphic):
    cdef __weakref__
    cdef readonly list _cache_bcys
    cdef readonly dict _cache_edge_orbits
    cdef readonly dict _cache_edge_pair_orbits
    cdef readonly dict _cache_isomorphisms
    cdef readonly list _cache_starting_vertices
    cdef readonly dict _cache_valence_spectrum
    cdef readonly list boundary_cycles
    cdef readonly int genus
    cdef readonly int num_boundary_cycles
    cdef readonly int num_edges
    cdef readonly int num_external_edges
    cdef readonly int num_vertices
    cdef readonly list edge_numbering
    cdef readonly list endpoints_i
    cdef readonly list endpoints_v
    cdef readonly list vertices
    cpdef list automorphisms(self)
    cpdef list compute_boundary_cycles(self)
    cdef int _cmp_orient(self, Fatgraph other, iso)
    cdef bridge(self, int edge1, int side1, int edge2, int side2)
    cdef Fatgraph contract_fg(self, int edgeno)
    cdef dict edge_orbits(self)
    cdef dict edge_pair_orbits(self)
    cdef endpoints(self, int edgeno)
    cdef Fatgraph hangcircle(self, int edge, int side)
    cpdef bint is_loop(self, int edge)
    cpdef bint is_orientation_reversing(self, automorphism)
    cpdef bint is_oriented(self)
    cpdef list isomorphisms(g1, g2)
    cpdef int num_automorphisms(self)
    cpdef int num_boundary_cycles(self)
    cdef other_end(graph, int edge, vertex, int attachment)
    cpdef int projection(self, other)
    cpdef dict valence_spectrum(self)
    cpdef frozenset vertex_valences(self)
    cpdef dict vertex_valence_distribution(self)

cdef class BoundaryCycle(frozenset):
    cdef readonly graph  # weakref proxy to Fatgraph
    cpdef BoundaryCycle set_graph(self, Fatgraph graph)
    cpdef contract_bcy(self, vi1, vi2, Fatgraph graph)
    cpdef transform(self, iso, Fatgraph graph)

cpdef list MgnTrivalentGraphsRecursiveGenerator(int g, int n)

cdef class MgnGraphsIterator(BufferingIterator):
    cdef readonly int g
    cdef readonly int n
    cdef readonly int _num_vertices
    cdef list _batch
    cdef list refill(self)
