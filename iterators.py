#! /usr/bin/env python
#
"""Iterator classes and utility functions.
"""
__docformat__ = 'reStructuredText'

import cython

from collections import Iterator


## main content

@cython.cclass
class BufferingIterator(Iterator):
    """Iterate over items stored in an internal buffer; when all items
    in the buffer have been handed out to caller, refill the buffer by
    calling `self.refill()` and start over again.

    This is intended as a base class for iterators that can generate
    more than one value per invocation; still, by the iterator
    protocol, they should return only one value to caller.  Subclasses
    of `BufferingIterator` should only need to define the `refill()`
    method, returning a list (or other iterable) with items that
    should be inserted in the buffer.

    The base class implementation just returns the items passed in the
    `initial` constructor argument and then raise `StopIteration`::

      >>> b = BufferingIterator([1,2,3])
      >>> for x in b: print x
      1
      2
      3

      >>> list(BufferingIterator())
      []
      
    """
    
    __slots__ = ('__buffer',)

    def __init__(self, initial=None):
        """Create a `BufferingIterator` instance and fill the internal
        buffer with items from `initial` (if supplied).
        """
        if initial is None:
            self.__buffer = []
        else:
            self.__buffer = list(initial)


    def next(self):
        """Return next item from queue, refilling queue if empty."""
        # try to refill buffer if empty
        if 0 == len(self.__buffer):
            self.__buffer.extend(self.refill())

        # if still empty after refill, then iteration has ended
        if 0 == len(self.__buffer):
            raise StopIteration

        return self.__buffer.pop(0)


    def refill(self):
        """Return new items to store in the buffer.

        At end of iteration, `refill` may either raise
        `StopIteration`, or just return an empty list.

        Sub-classes should override this method: the default
        implementation just signals `StopIteration`.
        """
        raise StopIteration


@cython.cclass
class chunks(Iterator):
    """Lump items from iterable into chunks of specified size.

    Instanciate the iterator passing a sequence of chunk sizes in
    argument 1 and an iterable to consume in argument 2::

      >>> for c in chunks([1,1,1], xrange(3)): print c
      [0]
      [1]
      [2]

    The list of chunk sizes may be any kind of sequence, for instance
    a tuple or even a (possibly infinite) iterable::
    
      >>> list(chunks((1,2,3), range(6)))
      [[0], [1, 2], [3, 4, 5]]

    The total size of the chunks may be less than the size of the
    iterator: remaining items in the iterator are not consumed::

      >>> for c in chunks([1,2], range(6)): print c
      [0]
      [1, 2]

    As a special case, if a chunk has size 0, then an empty list is
    returned in its place and no item from iterable is consumed::
    
      >>> for c in chunks([2,0,2], range(4)): print c
      [0, 1]
      []
      [2, 3]
      
    """
    def __init__(self, sizes, iterable):
        """Constructor, taking sequence of chunk sizes and iterable to
        consume."""
        self.current_chunk = -1
        self.sizes = sizes
        self.iterable = iter(iterable)

    @cython.locals(x=cython.int)
    def next(self):
        """Return next chunk."""
        self.current_chunk += 1
        if self.current_chunk >= len(self.sizes):
            raise StopIteration
        return [ self.iterable.next()
                   for x in xrange(self.sizes[self.current_chunk]) ]


@cython.cclass
class IndexedIterator(Iterator):
    """Return the items corresponding to indices `start`, `start+1`, etc.
    in the initialization sequence `lst`.

    Iteration stops as soon as an `IndexError` (indicating
    out-of-bounds) is returned.

    Examples::

      >>> lst = [0, 1, 2, 3]
      >>> for x in IndexedIterator(lst): print x,
      0 1 2 3
      >>> for y in IndexedIterator(lst, 2): print y,
      2 3
      
    """
    def __init__(self, lst, start=0):
        self.__lst = lst
        self.__cur = start-1

    def next(self):
        self.__cur += 1
        try:
            return self.__lst[self.__cur]
        except IndexError:
            raise StopIteration



## main: run tests

if "__main__" == __name__:
    import doctest
    doctest.testmod(name="iterators",
                    optionflags=doctest.NORMALIZE_WHITESPACE)
