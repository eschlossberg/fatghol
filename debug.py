#! /usr/bin/env python
#
"""Miscellaneous support functions for `assert` statements and
run-time debugging.
"""
__docformat__ = 'reStructuredText'


import operator
import sys
import types


def is_sequence_of_type(t, seq):
    """Return `True` if all items of sequence `s` are of type `t`.

    Examples::
      >>> is_sequence_of_type(types.IntType, [1,2,3])
      True
      >>> is_sequence_of_type(types.IntType, [1,2,"xxx"])
      False
      >>> is_sequence_of_type(types.StringType, ["xxx","yyy"])
      True
    """
    def is_type_t(item):
        return (type(item) is t)
    return reduce(operator.and_, map(is_type_t, seq), True)


def is_sequence_of_integers(seq):
    """Return `True` if all items of sequence `s` are of type `IntType`.

    Examples::
      >>> is_sequence_of_integers([1,2,3])
      True
      >>> is_sequence_of_integers([1,2,"xxx"])
      False
    """
    return is_sequence_of_type(types.IntType, seq)


def trace(func):
    def __tracing_wrapper(*args, **kwargs):
        print "DEBUG: Entering %s(%s, %s)" % (func.__name__, args, kwargs)
        result = func(*args, **kwargs)
        print "DEBUG: Got result %s from %s(%s, %s)" % (result, func.__name__, args, kwargs)
        return result
    return __tracing_wrapper




## main: run tests

if "__main__" == __name__:
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
