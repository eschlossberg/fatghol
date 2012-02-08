#! /usr/bin/env python
#
"""A collection of small utility functions.

These were mostly ripped out of `rg.py` for readability.
"""
#
#   Copyright (C) 2008-2012 Riccardo Murri <riccardo.murri@gmail.com>
#   All rights reserved.
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
__docformat__ = 'reStructuredText' 


#import cython

## stdlib imports

from collections import Iterator
import itertools
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

      >>> concat(['ab','c'])
      'abc'
    """
    return reduce(operator.add, seqs)


#@cython.cclass
class itranslate(Iterator):
    """Return items from a sequence, substituting them as specified.

    First argument `subst` is a dictionary, specifying substitutions
    to be applied.  If an item matches a key of the `subst`
    dictionary, the associated dictionary value is returned instead;
    unless the value is `None`, in which case the item is skipped
    altogether.

    *Note:* you should use an appropriate `dict`-subclass if you want
     to translate items which are not immutable.
    
    Examples::
      >>> list(itranslate({0:None, 3:2}, [2,1,0,0,1,3]))
      [2, 1, 1, 2]
    """

    __slots__ = ('mappings', 'iterable')
    
    def __init__(self, subst, iterable):
        self.mappings = subst
        self.iterable = iter(iterable)
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


#@cython.ccall(list)
def lconcat(seqs):
    """Return list concatenation of all sequences in `seqs`.

    Examples::
    
      >>> lconcat([[0]])
      [0]

      >>> lconcat([[0],[1]])
      [0, 1]

      >>> lconcat(['ab','c'])
      ['a', 'b', 'c']
    """
    return list(itertools.chain(*seqs))


#@cython.ccall(list)
def ltranslate(subst, iterable):
    """Return list of items from a sequence, substituting them as specified.

    First argument `subst` is a dictionary, specifying substitutions
    to be applied.  If an item matches a key of the `subst`
    dictionary, the associated dictionary value is returned instead;
    unless the value is `None`, in which case the item is skipped
    altogether.

    *Note:* you should use an appropriate `dict`-subclass if you want
     to translate items which are not immutable.
    
    Examples::
      >>> ltranslate({0:None, 3:2}, [2,1,0,0,1,3])
      [2, 1, 1, 2]
    """
    return [ (subst[item] if (item in subst) else item)
             for item in iterable
             if (item not in subst) or (subst[item] is not None) ]


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


#@cython.ccall(int)
#@cython.locals(arg=str, result=cython.int)
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


#@cython.locals(L=list, p=cython.int)
#@cython.ccall(list)
def rotated(L, p):
    """
    Return a copy of the given list `L`, with items rotated `p`
    positions leftwards.

    Examples::

      >>> rotated([1,2,3,4,5], 2)
      [3, 4, 5, 1, 2]
      
    """
    return L[p:] + L[:p]


#@cython.ccall(int)
#@cython.locals(x=cython.int)
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
