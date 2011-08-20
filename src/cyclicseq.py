#! /usr/bin/env python
#
"""Cyclic sequences.
"""
__docformat__ = 'reStructuredText'



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
            i = max(i, 0); j = max(j, 0)
            base.__delslice__(self, i%len(self), j%len(self))
        def __getslice__(self, i, j):
            """Return [i:j] slice, as a `"""+base.__class__.__name__+"""` instance."""
            i = max(i, 0); j = max(j, 0); l = len(self)
            if abs(i-j) < l:
                return base.__getslice__(self, i%l, j%l)
            else:
                a = max(i,j)
                b = min(i,j)
                c = l*(a/l)     # nearest multiple of l below a
                d = l*(b/l + 1) # nearest multiple of l above b
                n = (a/l) - (b/l + 1) # how many times the whole seq is repeated in the middle
                return base.__getslice__(self, b%l,l) \
                       + (n * base.__getslice__(self, 0,l)) \
                       + base.__getslice__(self, 0,a%l)
        def __setslice__(self, i, j, other):
            i = max(i, 0); j = max(j, 0); l = len(self)
            if isinstance(other, base):
                base.__setslice__(self, i%len(self), j%len(self), other)
            else:
                base.__setslice__(self, i%len(self), j%len(self), list(other))

        ## equality predicates
        ##
        def __eq__(self, other):
            """Return `True` if `self` is linearly equal to `other`
            with all indices shifted by a fixed amount."""
            if len(other) != len(self):
                return False
            elif None == self._shift_for_linear_eq(other):
                return False
            else:
                return True

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
            >>> c=CyclicList([1,2,3])
            >>> _repetition_pattern(c)
            (1, [1, 1, 1])
            >>> c=CyclicList([4,4,4])
            >>> _repetition_pattern(c)
            (0, [3])
            >>> c=CyclicArray([4,4,4,1])
            >>> _repetition_pattern(c, CyclicArray)
            (3, array('i', [1, 3]))
            >>> c=CyclicList([1,4,4,4,1])
            >>> _repetition_pattern(c, CyclicArray)
            (1, array('i', [3, 2]))
            >>> c=CyclicArray([4,4,4,3,3])
            >>> _repetition_pattern(c)
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
