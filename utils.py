#! /usr/bin/env python
#
"""A collection of small utility functions.

These were mostly ripped out of `rg.py` for readability.
"""
__docformat__ = 'reStructuredText' 


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
      >>> list(enumerate_set_product([[1],[1]]))
      [[1, 1]]
      >>> list(enumerate_set_product([[1,2],[1]]))
      [[1, 1], [2, 1]]
      >>> list(enumerate_set_product([[1,2],[1,2]]))
      [[1, 1], [2, 1], [1, 2], [2, 2]]
    """
    if len(p) == 0:
        yield []
    else:
        for i in p[-1]:
            for js in enumerate_set_product(p[:-1]):
                yield js+[i]


def other(pair, one):
    """Return the member of `pair` not equal to `one`."""
    if pair[0] == one:
        return pair[1]
    else:
        return pair[0]


def _tr(elt, t1, t2):
    """If `elt` equals some element in set `t1`, then return the corresponding element from set `t2`, otherwise return `elt` unchanged.
    """
    try:
        return t2[t1.index(elt)]
    except ValueError:
        return elt

def tr(s, t1, t2):
    """Change every occurrence (in sequence `s`) of an element of set `t1` with the corrisponding element of set `t2`.

    Examples::
      >>> tr([0,1,0,0],[0],[2])
      [2, 1, 2, 2]
    """
    return [_tr(x, t1, t2) for x in s]


def itr(s, t1, t2):
    """Change every occurrence (in iterable `s`) of an element of set `t1` with the corrisponding element of set `t2`.

    Examples::
      >>> list(tr([0,1,0,0],[0],[2]))
      [2, 1, 2, 2]
    """
    for x in s:
        yield _tr(x, t1, t2)



## main: run tests

if "__main__" == __name__:
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
