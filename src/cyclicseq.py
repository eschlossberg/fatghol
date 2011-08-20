#! /usr/bin/env python
#
"""Cyclic sequences.
"""
__docformat__ = 'reStructuredText'


from sys import maxint as SYS_MAXINT


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
            if SYS_MAXINT == j:
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
            if j == SYS_MAXINT:
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
            if SYS_MAXINT == j:
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
            return sum(hash(self[x:x+len(self)]) for x in xrange(len(self))) \
                   % SYS_MAXINT
        
        def __eq__(self, other):
            """Return `True` if `self` is linearly equal to `other`
            with all indices shifted by a fixed amount."""
            if len(other) != len(self):
                return False
            elif None == self._shift_for_linear_eq(other):
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
            """Return `True` if `first` is linearly equal to `second`
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

        def _shift_for_linear_eq(self, other, start=0):
            """Return minimum shift index `b >= start` such that
            `self[b:b+len]==other` as linear sequences.

            Examples::
              >>> a=CyclicList([1,2,3])
              >>> b=CyclicList([2,3,1])
              >>> a._shift_for_linear_eq(b)
              1
              >>> a._shift_for_linear_eq(b,2) is None
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

        def all_shifts_for_linear_eq(self, other):
            """Iterate over all shift amounts such that `self` is
            linearly equal to `other` when shifted by that amount.
            """
            start = 0
            l = len(self)
            while start < l:
                shift = self._shift_for_linear_eq(other, start)
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

        def repetition_pattern(self, cls=None):
          """Return the repetition pattern of a cyclic sequence.

          The repetition pattern of a cyclic sequence `c` is again a
          *cyclic* sequence of integers: each number `n` in the
          repetition list corresponds to `n` equal-valued items in
          `c`.

          The repetition pattern needs an "anchor point" to be meaningful:
          sequences `[0,1,0]` and `[0,0,1]` are equal as cyclic
          sequences and have the same repetition pattern `[2,1]`
          (which is == `[1,2]`), but one needs to know if 2 is the
          number of repetitions of `1`s or `0`s.

          So, the returned value is actually a pair `(anchor, rep)`
          where `rep` is the repetition pattern (a `CyclicList` instance)
          of the linearized sequences starting at index `anchor`.
          
          Examples::
            >>> CyclicList([1,2,3]).repetition_pattern()
            (1, [1, 1, 1])
            >>> CyclicList([4,4,4]).repetition_pattern()
            (0, [3])
            >>> CyclicList([4,4,4,1]).repetition_pattern()
            (3, [1, 3])
            >>> CyclicList([1,4,4,4,1]).repetition_pattern()
            (1, [3, 2])
            >>> CyclicList([4,4,4,3,3]).repetition_pattern()
            (3, [2, 3])
          """
          if self._repetition_pattern is None:
              if cls is None:
                  cls = CyclicList
              l=len(self)
              # find the first transition
              b=0
              while (b < l) and (self[b] == self[b+1]):
                  b += 1
              # all items are equal
              if b == l:
                  return (0, cls([l]))
              # else, start building pattern from here
              result = cls()
              b += 1
              i=0
              while i < l:
                  j=0
                  while (j < l) and (self[b+i+j] == self[b+i+j+1]):
                      j += 1
                  result.append(j+1)
                  i += j+1
              self._repetition_pattern = (b, result)
              
          return self._repetition_pattern
      
        def rotate(self, n):
            """Rotate sequence leftwards by `n` positions, *in-place*.

            Examples::
              >>> a=CyclicList([3,2,1])
              >>> a.rotate(1); a
              [2, 1, 3]

            If `n` is negative, rotate rightwards by `abs(n)`
            positions::
              >>> a.rotate(-1); a
              [3, 2, 1]

            Rotating by 0 or by a multiple of the list length, will
            leave the sequence unchanged::
              >>> a=CyclicList([3,2,1])
              >>> a.rotate(0); a
              [3, 2, 1]
              >>> a.rotate(3); a
              [3, 2, 1]
              >>> a.rotate(6); a
              [3, 2, 1]
            
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
    
