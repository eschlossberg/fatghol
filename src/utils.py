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
    """Perform substitutions on a iterable, just like `str.translate`.
    """
    def __init__(self, subst, iterable):
        self.mappings = subst
        self.iterable = iterable
    def __iter__(self):
        return self
    def next(self):
        while True:
            next = self.iterable.next()
            if not self.mappings.has_key(next):
                return next
            translated = self.mappings[next]
            if translated is None:
                # skip this item
                continue
            return translated


def _tr(elt, t1, t2):
    """If `elt` equals some element in set `t1`, then return the corresponding element from set `t2`, otherwise return `elt` unchanged.
    """
    try:
        return t2[t1.index(elt)]
    except ValueError:
        return elt

def tr(s, t1, t2):
    """Return a copy of sequence `s` where each occurence of an
    element of set `t1` has been changed with the corrisponding
    element of set `t2`.

    Examples::
      >>> tr([0,1,0,0],[0],[2])
      [2, 1, 2, 2]
    """
    return [_tr(x, t1, t2) for x in s]


def itr(s, t1, t2):
    """Return an iterator over `s` where each occurrence of an element
    of set `t1` is changed with the corrisponding element of set `t2`.

    Examples::
      >>> list(itr([0,1,0,0],[0],[2]))
      [2, 1, 2, 2]
    """
    for x in s:
        yield _tr(x, t1, t2)

def tr_inplace(s, t1, t2):
    """Change every occurrence of an element of set `t1` in sequence
    `s` with the corrisponding element of set `t2`.

    Return the altered `s`.
    
    Examples::
      >>> x=[0,1,0,0]
      >>> tr_inplace(x,[0],[2])
      [2, 1, 2, 2]
      >>> x
      [2, 1, 2, 2]
    """
    for i in xrange(len(s)):
        s[i] = _tr(s[i], t1, t2)
    return s



## main: run tests

if "__main__" == __name__:
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
