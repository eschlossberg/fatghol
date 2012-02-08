from collections import Iterator

cdef class BufferingIterator(Iterator):
    cdef __buffer
    cdef list refill(self)

cdef class chunks(Iterator):
    cdef int current_chunk
    cdef sizes
    cdef iterable

cdef class IndexedIterator(Iterator):
    cdef list __lst
    cdef int __cur
