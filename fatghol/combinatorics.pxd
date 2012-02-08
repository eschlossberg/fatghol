cpdef object bernoulli(int n)
cpdef object choose(int n, int k)
cpdef long factorial(int n)
cpdef int minus_one_exp(int m)
cdef class OrderedSetPartitionsIterator:
    cdef readonly list sizes
    cdef items
    cdef state
    cpdef list next(self)
    cdef _next_state(self)
cdef class SetProductIterator_:
    cdef bint __closed
    cdef object __factors
    cdef int __L
    cdef int __i
    cdef list __M
    cdef list __m
cdef class Permutation(dict):
    cpdef Permutation inverse(self)
    cpdef bint is_identity(self)
    cpdef list rearranged(self, seq)
    cpdef int sign(self)
    cpdef translate(self, seq)
    cpdef list ltranslate(self, iterable)
    cpdef bint extend(self, srcs, dsts)
    cpdef bint update(self, dict mappings)    
cdef class PermutationList:
    cdef int __order
    cdef set __base_set
cdef class PermutationIterator:
    cdef list seq
    cdef long rank
    cdef int order
cdef class InplacePermutationIterator:
    cdef list seq
    cdef int start
    cdef int end
    cdef bint enumeration_finished
cdef class FixedLengthPartitionIterator:
    cdef long _N
    cdef int _K
    cdef int _k
    cdef long _min
    cdef long _max
    cdef list _p
    cdef bint done

cpdef PartitionIterator(N, K, min_=?, max_=?)

