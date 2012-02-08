cdef class AggregateList(object):
     cdef list __components
     cdef list __lengths
     cpdef aggregate1(self, seq)
     cdef append(self, item)
     cdef extend(self, seq)
     cdef int index(self, value, int start=?)
     cpdef iterblocks(self)
     cpdef itervalues(self)
