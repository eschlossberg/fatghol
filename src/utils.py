#! /usr/bin/env python
#
"""A collection of small utility functions.

These were mostly ripped out of `rg.py` for readability.
"""
__docformat__ = 'reStructuredText' 


## stdlib imports

import operator
import types


## main content

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


## conditional application of decorators
def maybe(deco, default=False, cond=None):
    if cond is None:
        try:
            cond = deco.enabled
        except AttributeError:
            cond = default
    if cond:
        return deco
    else:
        return (lambda fn: fn)


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


def sign(x):
    """Return the sign of `x`: that is, 0 if `x` is zero, +1 iff it is
    positive, -1 iff it is negative.

    Examples::

      >>> sign(0)
      0
      >>> sign(3)
      1
      >>> sign(-5)
      -1
      
    """
    if x == 0:
        return 0
    elif x == abs(x):
        return +1
    else:
        return -1


## main: run tests

if "__main__" == __name__:
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
