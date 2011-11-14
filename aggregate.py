#! /usr/bin/env python
#
"""
The `AggregateList` class.
"""
__docformat__ = 'reStructuredText'

import cython

import itertools

@cython.cclass
class AggregateList(object):
    """Act on a set of lists as if they were concatenated together.

    The aggregate list acts as a collective: the length is the sum of
    all the components' length, iteration over the composite list
    reveals all the components members in the order the components
    were added.  References to elements of the aggregate list are
    routed to the component lists, and iteration over the aggregate
    list returns the elements of the component lists (in the order the
    components were added to the aggregate)::

      >>> a1 = [0,1,2]
      >>> a2 = [3, 4]
      >>> a = AggregateList(a1, a2)
      >>> len(a)
      5
      >>> a[1]
      1
      >>> a[4]
      4
      >>> a[4] = 77
      >>> 77 in a
      True
      >>> list(a)
      [0, 1, 2, 3, 77]

    Operations are routed to the appropriate component list: e.g.,
    `append` appends to the last list in the group.

      >>> a.append(5)
      >>> a2
      [3, 77, 5]

    *Note:* Currently not implemented:
      - remove()
      - pop()
      - count()
      - sort()
      - reverse()
      - All slicing operations.
    """
    def __init__(self, *seqs):
        """Constructor, taking initial list of
        """
        self.__components = list(seqs)
        self.__lengths = [ len(s) for s in seqs ]

    def __contains__(self, value):
        for component in self.__components:
            if component.__contains__(value):
                return True

    @cython.locals(i=cython.int, l=cython.int, n=cython.int)
    def __delitem__(self, i):
        for (n, component) in enumerate(self.__components):
            l = len(component)
            if i >= l:
                i -= l
            else:
                component.__delitem__(i)
                self.__lengths[n] -= 1
        raise IndexError("AggregateList.__delitem__(): list assignment out of range")

    @cython.locals(i=cython.int, l=cython.int)
    def __getitem__(self, i):
        for component in self.__components:
            l = len(component)
            if i >= l:
                i -= l
            else:
                return component.__getitem__(i)
        raise IndexError("AggregateList.__getitem__(): list assignment out of range")

    @cython.ccall
    def iterblocks(self):
        return iter(self.__components)

    @cython.ccall
    def itervalues(self):
        return itertools.chain(* self.__components)
    def __iter__(self):
        return self.itervalues()

    def __len__(self):
        return sum(self.__lengths)

    @cython.locals(i=cython.int, l=cython.int)
    def __setitem__(self, i, value):
        for component in self.__components:
            l = len(component)
            if i >= l:
                i -= l
            else:
                component.__setitem__(i, value)
                return
        raise IndexError("AggregateList.__setitem__(): list assignment out of range")

    def aggregate(self, *seqs):
        """Append each sequence in `seqs` to the aggregated list.

        Example::

          >>> a = AggregateList([0,1,2])
          >>> a.aggregate([3,4],[5,6])
          >>> list(a)
          [0, 1, 2, 3, 4, 5, 6]
        """
        for s in seqs:
            self.__components.append(s)
            self.__lengths.append(len(s))

    @cython.ccall
    def aggregate1(self, seq):
        """Append sequence `seq` to the aggregate list.

        Example::

          >>> a = AggregateList([0,1,2])
          >>> a.aggregate([3,4])
          >>> a.aggregate([5,6])
          >>> list(a)
          [0, 1, 2, 3, 4, 5, 6]
        """
        self.__components.append(seq)
        self.__lengths.append(len(seq))

    @cython.cfunc
    def append(self, item):
        """Append `item` to the last component."""
        self.__components[-1].append(item)
        self.__lengths[-1] += 1
        
    @cython.cfunc
    def extend(self, seq):
        """Extend the last component with `seq`."""
        self.__components[-1].extend(seq)
        self.__lengths[-1] = len(self.__components[-1])

    #@cython.cfunc(cython.int)
    @cython.locals(index=cython.int, l=cython.int)
    def index(self, value, start=0):
        """Return index of the first element equal to `value`."""
        index = 0
        for component in self.__components:
            l = len(component)
            if start > l:
                start -= l
                index += l
                continue
            try:
                return index + component.index(value, start)
            except ValueError:
                start = 0
                continue
        raise ValueError("AggregateList.index(): x not in list")
