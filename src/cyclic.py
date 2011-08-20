#! /usr/bin/env python
#
"""Cyclic data structures: lists, arrays, and related functions.
"""
__docformat__ = 'reStructuredText' 


from array import array, ArrayType


class CyclicList(list):
    """List with indices wrapping around.

    Examples::
      >>> a=CyclicList([1,2,3])
      >>> len(a)
      3
      >>> a[3]
      1
      >>> a[4]
      2
      >>> a[4]=5
      >>> a
      [1, 5, 3]
      >>> a[1:2]
      [5]
      >>> a[1:2]=[7,8,9]
      >>> a
      [1, 7, 8, 9, 3]
      >>> len(a)
      5
      >>> a=CyclicList([1,2,3])
      >>> b=CyclicList([2,3,1])
      >>> b == a
      True
      >>> c=CyclicList([3,1,2])
      >>> c==a
      True
      >>> d=CyclicList([1,3,2])
      >>> d==a
      False
    """
    def __init__(self, sequence=None):
        if sequence is None:
            list.__init__(self)
        else:
            list.__init__(self, sequence)

    def __getitem__(self, i): return list.__getitem__(self, i % len(self))
    def __setitem__(self, i, item): list.__setitem__(self, i % len(self), item)
    def __delitem__(self, i): list.__delitem__(self, i%len(self))

    def __getslice__(self, i, j):
        """*Note:* returned slice is of `list` type!"""
        i = max(i, 0); j = max(j, 0); l = len(self)
        if abs(i-j) < l:
            return self.__class__(list.__getslice__(self, i%l, j%l))
        else:
            a = max(i,j)
            b = min(i,j)
            c = l*(a/l)     # nearest multiple of l below a
            d = l*(b/l + 1) # nearest multiple of l above b
            n = (a/l) - (b/l + 1) # how many times the whole seq is repeated in the middle
            return list.__getslice__(self, b%l,l) \
                   + (n*list.__getslice__(self, 0,l)) \
                    + list.__getslice__(self, 0,a%l)
    def __setslice__(self, i, j, other):
        i = max(i, 0); j = max(j, 0); l = len(self)
        if isinstance(other, list):
            list.__setslice__(self, i%len(self), j%len(self), other)
        else:
            list.__setslice__(self, i%len(self), j%len(self), list(other))
    def __delslice__(self, i, j):
        i = max(i, 0); j = max(j, 0)
        list.__delslice__(self, i%len(self), j%len(self))

    def __eq__(self, other):
        """Compare `self` with all possible translations of `other`."""
        if len(other) != len(self):
            return False
        elif None == self.shift_for_list_eq(other):
            return False
        else:
            return True

    def shift_for_list_eq(self, other, start=0):
        """Return minimum shift index `b >= start` such that `self[b:b+len]==other` as Python lists.

        Examples::
          >>> a=CyclicList([1,2,3])
          >>> b=CyclicList([2,3,1])
          >>> a.shift_for_list_eq(b)
          1
          >>> a.shift_for_list_eq(b,2) is None
          True
        """
        l = len(self)
        def _eq_shifted(first, second, shift):
            l=len(first)
            i=0
            while i < l:
                if first[i+shift] != second[i]:
                    return False
                else:
                    i += 1
            return True
        shift = start 
        while shift < l:
            if _eq_shifted(self, other, shift):
                return shift
            else:
                shift += 1
        return None
    def all_shifts_for_list_eq(self, other):
        start = 0
        l = len(self)
        while start < l:
            shift = self.shift_for_list_eq(other, start)
            if shift is None:
                break
            else:
                yield shift
            start = shift+1


class CyclicArray(ArrayType):
    """Array with indices wrapping around.

    Examples::
      >>> a=CyclicArray('i', [1,2,3])
      >>> len(a)
      3
      >>> a[3]
      1
      >>> a[4]
      2
      >>> a[4]=5
      >>> a
      array('i', [1, 5, 3])
      >>> a[1:2]
      array('i', [5])
      >>> a[1:2]=array('i', [7,8,9])
      >>> a
      array('i', [1, 7, 8, 9, 3])
      >>> len(a)
      5
      >>> a=CyclicArray('i', [1,2,3])
      >>> b=CyclicArray('i', [2,3,1])
      >>> b == a
      True
      >>> c=CyclicArray('i', [3,1,2])
      >>> c==a
      True
      >>> d=CyclicArray('i', [1,3,2])
      >>> d==a
      False
    """
    def __new__(cls, sequence=None, typecode='i', *args, **kwargs):
        if sequence is None:
            return ArrayType.__new__(cls, typecode, *args, **kwargs)
        else:
            return ArrayType.__new__(cls, typecode, sequence, *args, **kwargs)

    def __getitem__(self, i): return array.__getitem__(self, i % len(self))
    def __setitem__(self, i, item): array.__setitem__(self, i % len(self), item)
    def __delitem__(self, i): array.__delitem__(self, i%len(self))

    def __getslice__(self, i, j):
        """*Note:* returned slice is of `array` type!"""
        i = max(i, 0); j = max(j, 0); l = len(self)
        if abs(i-j) < l:
            return self.__class__(array.__getslice__(self, i%l, j%l))
        else:
            a = max(i,j)
            b = min(i,j)
            c = l*(a/l)     # nearest multiple of l below a
            d = l*(b/l + 1) # nearest multiple of l above b
            n = (a/l) - (b/l + 1) # how many times the whole seq is repeated in the middle
            return array.__getslice__(self, b%l,l) \
                   + (n*array.__getslice__(self, 0,l)) \
                    + array.__getslice__(self, 0,a%l)
    def __setslice__(self, i, j, other):
        i = max(i, 0); j = max(j, 0); l = len(self)
        if isinstance(other, array):
            array.__setslice__(self, i%len(self), j%len(self), other)
        else:
            array.__setslice__(self, i%len(self), j%len(self),
                               array(self.typecode, other))
    def __delslice__(self, i, j):
        i = max(i, 0); j = max(j, 0)
        array.__delslice__(self, i%len(self), j%len(self))

    def __eq__(self, other):
        """Compare `self` with all possible translations of `other`."""
        if len(other) != len(self):
            return False
        elif None == self.shift_for_list_eq(other):
            return False
        else:
            return True

    def shift_for_list_eq(self, other, start=0):
        """Return minimum shift index `b >= start` such that `self[b:b+len]==other` as Python lists.

        Examples::
          >>> a=CyclicArray([1,2,3])
          >>> b=CyclicArray([2,3,1])
          >>> a.shift_for_list_eq(b)
          1
          >>> a.shift_for_list_eq(b,2) is None
          True
        """
        l = len(self)
        def _eq_shifted(first, second, shift):
            l=len(first)
            i=0
            while i < l:
                if first[i+shift] != second[i]:
                    return False
                else:
                    i += 1
            return True
        shift = start 
        while shift < l:
            if _eq_shifted(self, other, shift):
                return shift
            else:
                shift += 1
        return None
    def all_shifts_for_list_eq(self, other):
        start = 0
        l = len(self)
        while start < l:
            shift = self.shift_for_list_eq(other, start)
            if shift is None:
                break
            else:
                yield shift
            start = shift+1


class RotatedList(CyclicList):
    """Access a given `CyclicList` with indices displaced by a fixed amount.

    Examples::
      >>> a = [1, 7, 8, 9, 3]
      >>> b=RotatedList(1,a)
      >>> len(b)
      5
      >>> b[0]
      7
      >>> b[5]
      7
      >>> b[1]
      8
      >>> b[0]=6
      >>> b[0]
      6
    """
    def __init__(self, displacement, initial=None):
        self.shift = displacement
        CyclicList.__init__(self, initial)
    def __getitem__(self, i): return CyclicList.__getitem__(self, i + self.shift)
    def __setitem__(self, i, item): CyclicList.__setitem__(self, i + self.shift, item)
    def __delitem__(self, i): CyclicList.__delitem__(self, i + self.shift)
    def __getslice__(self, i, j):
        return CyclicList.__getslice__(self, i+self.shift, j+self.shift)
    def __setslice__(self, i, j, other):
        CyclicList.__setslice__(self, i+self.shift, j+self.shift, other)
    def __delslice__(self, i, j):
        CyclicList.__delslice__(self, i+self.shift, j+self.shift)


def repetition_pattern(c):
    """Return the repetition pattern of a cyclic list `c`.

    The repetition pattern is a *cyclic* list of integers: each number
    `n` in the repetition list corresponds to `n` equal-valued items
    in the list `c`.

    Examples::
      >>> c=CyclicList([1,2,3])
      >>> repetition_pattern(c)
      (1, [1, 1, 1])
      >>> c=CyclicList([4,4,4])
      >>> repetition_pattern(c)
      (0, [3])
      >>> c=CyclicList([4,4,4,1])
      >>> repetition_pattern(c)
      (3, [1, 3])
      >>> c=CyclicList([1,4,4,4,1])
      >>> repetition_pattern(c)
      (1, [3, 2])
      >>> c=CyclicList([4,4,4,3,3])
      >>> repetition_pattern(c)
      (3, [2, 3])
    """
    l=len(c)
    # find the first transition
    b=0
    while (b < l) and (c[b] == c[b+1]):
        b += 1
    # all items are equal
    if b == l:
        return (0, CyclicList([l]))
    # else, start building pattern from here
    result=CyclicList()
    b += 1
    i=0
    while i < l:
        j=0
        while (j < l) and (c[b+i+j] == c[b+i+j+1]):
            j += 1
        result.append(j+1)
        i += j+1
    return (b, result)



## main: run tests

if "__main__" == __name__:
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
