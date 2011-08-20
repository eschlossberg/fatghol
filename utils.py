#! /usr/bin/env python
#
"""A collection of small utility functions.

These were mostly ripped out of `rg.py` for readability.
"""
__docformat__ = 'reStructuredText' 


import operator
import types



def concat(seqs):
    """Return concatenation of all sequences in `seqs`.

    Examples::
      >>> concat([[0]])
      [0]
      >>> concat([[0],[1]])
      [0, 1]
      >>> concat(['a','b'])
      'ab'
    """
    return reduce(operator.add, seqs)
        

def deep_cmp(s1,s2):
    """Compare items in `s1` and `s2`, recursing into subsequences.

    Examples::
      >>> deep_cmp(1,1)
      0
      >>> deep_cmp([1],[1])
      0
      >>> deep_cmp([1,1],[1,1])
      0
      >>> deep_cmp([1,[1]],[1,[1]])
      0
      >>> deep_cmp([1,[1]],[1,[2]])
      -1
    """
    if not (type(s1) == type(s2)):
        raise TypeError, \
            "Comparing arguments of different type: %s vs %s" \
            % (repr(type(s1)), repr(type(s2)))
    else:
        try:
            # assume s1,s2 are sequences and recursively apply this
            # function to pairs of corresponding elements...
            def _first_nonzero(x,y):
                if 0 != x:
                    return x
                else:
                    return y
            return reduce(_first_nonzero, map(deep_cmp, s1, s2), 0)
        except TypeError:
            # ...if s1,s2 are not sequences, then do a builtin comparison
            return cmp(s1,s2)


def enumerate_set_product(p):
    """Iterate over all elements in the cartesian products of elements of items in `p`.

    Examples::
      >>> list(enumerate_set_product([]))
      [[]]
      >>> list(enumerate_set_product([[1]]))
      [[1]]
      >>> list(enumerate_set_product([[1],[1]]))
      [[1, 1]]
      >>> list(enumerate_set_product([[1,2],[]]))
      [[]]
      >>> list(enumerate_set_product([[1,2],[1]]))
      [[1, 1], [2, 1]]
      >>> list(enumerate_set_product([[1,2],[1,2]]))
      [[1, 1], [2, 1], [1, 2], [2, 2]]
    """
    L = len(p)
    M = [ len(s)-1 for s in p ]
    if (0 == L) or (-1 in M):
        # there are no factors, or one of them has no elements
        yield []
    else:
        m = [0] * L
        i = 0
        while i < L:
            # return element corresponding to current multi-index
            yield [ s[m[i]] for (i,s) in enumerate(p) ]
            # advance multi-index
            i = 0
            while (i < L):
                if m[i] == M[i]:
                    m[i] = 0
                    i += 1
                else:
                    m[i] += 1
                    break


def positive_int(arg):
    """Convert a string or number to a positive integer, if possible.
    Behaves just like the built-in `int` (which see), and additionally
    raises `ValueError` if the converted integer is less-then or equal
    to 0.
    """
    result = int(arg)
    if result <= 0:
        raise ValueError("non-positive integer literal: %d" % result)
    return result


class itranslate:
    """Return items from a sequence, substituting them as specified.

    First c'tor argument `subst` is a dictionary, specifying
    substitutions to be applied.  If an item matches a key of the
    `subst` dictionary, the associated dictionary value is returned
    instead; unless the value is `None`, in which case the item is
    skipped altogether.

    *Note:* you should use an appropriate `dict`-subclass if you want
     to translate items which are not immutable.
    
    Examples::
      >>> list(itranslate({0:None, 3:2}, [2,1,0,0,1,3]))
      [2, 1, 1, 2]
    """
    def __init__(self, subst, iterable):
        self.mappings = subst
        self.iterable = iter(iterable)
    def __iter__(self):
        return self
    def next(self):
        while True:
            next = self.iterable.next()
            if next not in self.mappings:
                return next
            translated = self.mappings[next]
            if translated is None:
                # skip this item
                continue
            return translated



## main: run tests

if "__main__" == __name__:
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
