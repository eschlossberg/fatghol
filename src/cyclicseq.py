#! /usr/bin/env python
#
"""Cyclic sequence factories.

Defines a type factory `CyclicSequence`, which can be applied to any
sequence-like type to produce a"cyclic" variant, that is, indices wrap
around: if `a` is such a cyclic sequence, then `a[i] == a[i + len(a)]`
for any index `i`.  For example::

  >>> clist = CyclicSequence(list)
  >>> a = clist([0,1,2,3])
  >>> len(a)
  4
  >>> a[1] == a[5]
  True
  >>> a[1] == a[9]
  True

More generally, `a[i] == a[j]` if `i = j mod len(a)`::

  >>> a[5] # == a[2]
  1
  >>> a[-7] # == a[2]
  1

By default, this module exports cyclic variants of the standard
`tuple` and `list` sequences::

  CyclicTuple = CyclicSequence(tuple)
  CyclicList = CyclicSequence(list)

`CyclicSequence` assumes that an instance of the base (non-cyclic)
type can be instanciated by passing to it only one argument, namely,
the initialization sequence.  (This is how the standard sequence types
`list` and `tuple` behave.)  For types needing a more complex
construction, a factory function must be passed as second argument:
the factory function accepts the initialization sequence and returns
an instance of the base type.  For example, this is how one can define
cyclic variants of the sequences provided by the `array` module::

  # cyclic versions of sequences in the `array` module
  from array import array
  CyclicIntArray = CyclicSequence(array, lambda seq: array('i', seq))
  CyclicUIntArray = CyclicSequence(array, lambda seq: array('I', seq))
  CyclicLongArray = CyclicSequence(array, lambda seq: array('l', seq))
  CyclicULongArray = CyclicSequence(array, lambda seq: array('L', seq))

"""
__docformat__ = 'reStructuredText'


import sys


def CyclicSequence(base, factory=None):
    if factory is None:
        factory=base
        
    class _CyclicSequence(base):

        def __init__(self, sequence=None):
            self._repetition_pattern = None
            if sequence is None:
                base.__init__(self)
            else:
                base.__init__(self, sequence)

        def __repr__(self):
            return "%s(%s)" % (self.__class__.__name__,
                               base.__repr__(self))

        ## item accessors
        ##        
        def __delitem__(self, i):
            base.__delitem__(self, i % len(self))
        def __getitem__(self, i):
            return base.__getitem__(self, i % len(self))
        def __setitem__(self, i, item):
            base.__setitem__(self, i % len(self), item)

        ## slice accessors
        ##    
        def __delslice__(self, i, j):
            if sys.maxint == j:
                # __delslice__ is passed `sys.maxint` as second value
                # when taking slices of the form `[n:]`.  In this
                # case, just fall back to standard Python operation.
                return base.__setslice__(self, i, j)
            l = len(self)
            #: sup of the range
            a = max(i,j)
            #: inf of the range
            b = min(i,j)
            #: highest multiple of `l` below `a`
            c = l*(a//l)
            a -= c
            b -= c
            return base.__delslice__(self, b, a)

        def __getslice__(self, i, j):
            """Return [i:j] slice as a base class instance."""
            if j == sys.maxint:
                # __getslice__ is passed `sys.maxint` as second value
                # when taking slices of the form `[n:]`.  In this
                # case, just fall back to standard Python operation.
                return base.__getslice__(self, i, j)
            l = len(self)
            #: sup of the range
            a = max(i,j)
            #: inf of the range
            b = min(i,j)
            if (b >= 0) and (a-b < l):
                #: highest multiple of `l` below `a`
                c = l*(a//l)
                a -= c
                b -= c
                return base.__getslice__(self, b, a)
            else:
                #: how many times the whole seq is repeated in the middle
                n = (a//l) - (b//l + 1) 
                return base.__getslice__(self, b%l,l) \
                       + (n * base.__getslice__(self, 0,l)) \
                       + base.__getslice__(self, 0,a%l)

        def __setslice__(self, i, j, other):
            if sys.maxint == j:
                # __setslice__ is passed `sys.maxint` as second value
                # when taking slices of the form `[n:]`.  In this
                # case, fall back to standard Python operation.
                return base.__setslice__(self, i, j, other)
            l = len(self)
            #: sup of the range
            a = max(i,j)
            #: inf of the range
            b = min(i,j)
            #: highest multiple of `l` below `a`
            c = l*(a//l)
            a -= c
            b -= c
            return base.__setslice__(self, b, a, other)

        ## equality predicates
        ##
        def __hash__(self):
            """Return sum of items in the sequence.

            This is invariant for all permutations, not just the
            cyclic ones, but it is fast and easy to compute.
            """
            try:
                # this works for numeric sequences
                return int(sum(self))
            except:
                # slower but more general version works for all kinds of sequences
                return reduce(operator.xor, [hash(x) for x in self])
        
        def __eq__(self, other):
            """Return `True` if `self` is sequentially equal to `other`
            with all indices shifted by a fixed amount."""
            if len(other) != len(self):
                return False
            elif None == self._shift_for_sequential_eq(other):
                return False
            else:
                return True

        # both `__eq__` and `__ne__` are needed for testing equality of objects;
        # see `<http://www.voidspace.org.uk/python/articles/comparison.shtml>`
        def __ne__(self, other):
            """The opposite of `__eq__` (which see)."""
            return not self.__eq__(other)

        @staticmethod
        def _eq_shifted(first, second, shift):
            """Return `True` if `first` is sequentially equal to `second`
            when shifting all indices by the fixed `shift` amount.
            """
            l=len(first)
            i=0
            while i < l:
                if first[i+shift] != second[i]:
                    return False
                else:
                    i += 1
            return True

        def _shift_for_sequential_eq(self, other, start=0):
            """Return minimum shift index `b >= start` such that
            `self[b:b+len]==other` as sequential sequences.

            Examples::
              >>> a=CyclicList([1,2,3])
              >>> b=CyclicList([2,3,1])
              >>> a._shift_for_sequential_eq(b)
              1
              >>> a._shift_for_sequential_eq(b,2) is None
              True
            """
            l = len(self)
            shift = start 
            while shift < l:
                if self._eq_shifted(self, other, shift):
                    return shift
                else:
                    shift += 1
            return None

        def all_shifts_for_sequential_eq(self, other):
            """Iterate over all shift amounts such that `self` is
            sequentially equal to `other` when shifted by that amount.
            """
            start = 0
            l = len(self)
            while start < l:
                shift = self._shift_for_sequential_eq(other, start)
                if shift is None:
                    break
                else:
                    yield shift
                start = shift+1

        ## additions to the sequence type
        ##
        def downcast(self):
            """Return instance of base type, with same contents of this object."""
            return base(self)

        def rotate(self, n):
            """Rotate sequence leftwards by `n` positions, *in-place*.

            Examples::
              >>> a=CyclicList([3,2,1])
              >>> a.rotate(1); a
              CyclicList([2, 1, 3])

            If `n` is negative, rotate rightwards by `abs(n)`
            positions::
              >>> a.rotate(-1); a
              CyclicList([3, 2, 1])

            Rotating by 0 or by a multiple of the list length, will
            leave the sequence unchanged::
              >>> a=CyclicList([3,2,1])
              >>> a.rotate(0); a
              CyclicList([3, 2, 1])
              >>> a.rotate(3); a
              CyclicList([3, 2, 1])
              >>> a.rotate(6); a
              CyclicList([3, 2, 1])
            
            """
            n %= len(self)
            self[:] = self[n:] + self[:n]
          
    return type("Cyclic" + str(base.__name__).capitalize(),
                (_CyclicSequence, base) + base.__bases__,
                dict(_CyclicSequence.__dict__))


## cyclic versions of standard sequences
CyclicTuple = CyclicSequence(tuple)
CyclicList = CyclicSequence(list)

## cyclic versions of sequences in the `array` module
#from array import array
#CyclicIntArray = CyclicSequence(array, factory=lambda seq: array('i', seq))
#CyclicUIntArray = CyclicSequence(array, factory=lambda seq: array('I', seq))
#CyclicLongArray = CyclicSequence(array, factory=lambda seq: array('l', seq))
#CyclicULongArray = CyclicSequence(array, factory=lambda seq: array('L', seq))


## main: run tests

if "__main__" == __name__:
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
    
