#! /usr/bin/env python
#
"""Utilities for caching computation results.
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


## stdlib imports
from collections import defaultdict
import functools
from time import time
import weakref


## additional imports

#import cython

## local imports

from fatghol.iterators import Iterator


## auxiliary classes

#@cython.cclass
class _IteratorRecorder(object):
    """Cache results from a given iterator.  Client classes provided
    by the `replay()` method will then replay the iterator history;
    multiple players can replay from one recorded source
    independently.

    Example::

      >>> L = [1,2,3]
      >>> R = _IteratorRecorder(iter(L))
      >>> for x in R: print x
      1
      2
      3
      >>> for x in R.replay(): print x
      1
      2
      3
      
    **WARNING:** This is not thread-safe!
    """

    __slots__ = ['done', 'iterable', 'history', '__weakref__']
    
    def __init__(self, iterable):
        self.iterable = iterable
        self.history = []
        self.done = False

    def __iter__(self):
        return self.replay()

    def advance(self):
        """Record next item from the source iterator."""
        if self.done:
            raise StopIteration
        else:
            try:
                self.history.append(self.iterable.next())
            except StopIteration:
                self.done = True
                del self.iterable
                raise

    def replay(self):
        """Return a new player."""
        return _IteratorReplayer(self)
    

#@cython.cclass
class _IteratorReplayer(Iterator):
    """Replay values recorded into a given `_IteratorRecorder` class.
    Multiple players can replay from one recorded source
    independently.

    Instances of `_IteratorReplayer`:class: are only produced as a
    result of the `_IteratorRecorder.replay` method.  See
    `_IteratorRecorder`:class: for examples.

    **WARNING:** This is not thread-safe!
    """

    __slots__ = ['done', 'master', 'pos']
    
    def __init__(self, master):
        self.master = master
        self.pos = -1
        self.done = False

    def next(self):
        if self.done:
            raise StopIteration
        try:
            if self.pos == len(self.master.history) - 1:
                self.master.advance()
            self.pos += 1
            return self.master.history[self.pos]
        except StopIteration:
            self.done = True
            del self.master
            raise


#@cython.cclass
class Caching(object):
    """Instances of this class provide an interface for use by caches
    based on a `dict` subclass.
    """

    __slots__ = [
        '_cache0',             # cache for nullary methods
        '_cache_contract',     # cache `Fatgraph.contract` results
        '_cache_eq',           # cache `Fatgraph.__eq__` results
        '_cache_isomorphisms', # cache `Fatgraph.isomorphisms` results
        ]



## caching functions

# store 
_func_cache = defaultdict(dict)

#@cython.ccall
def fcache(func):
    """Cache result of a generic function.

    This decorator can cache results of calls `func(*args)`.

    CAVEATS:

    1. The whole argument tuple is cached, so any object referenced
      there will *not* be garbage collected!
    2. No keyword arguments are allowed in the cached function.
    """
    @functools.wraps(func)
    def wrapper(*args):
        cache = _func_cache[id(func)]
        key = args
        try:
            return cache[key]
        except KeyError:
            result = func(*args)
            cache[key] = result
            return result
    return wrapper


#@cython.ccall
def ocache0(func):
    """Cache result of a nullary object method.
    
    This decorator can cache results of calls `obj.method()`. The
    result cache is held in the object itself; therefore, to cache
    result from methods of objects using a '__slots__' declaration, a
    '_cache' slot must be present and writable.
    """
    @functools.wraps(func)
    def wrapper(obj):
        try:
            cache = obj._cache0
        except AttributeError:
            obj._cache0 = cache = dict()
        try:
            return cache[func.func_name]
        except KeyError:
            result = func(obj)
            cache[func.func_name] = result
            return result
    return wrapper


#@cython.ccall
def ocache_contract(func):
    """Cache result of the `Fatgraph.contract` method.
    
    Only a weak reference to the method's result is held, so the cache
    entry is automatically freed when the result Fatgraph goes out of
    scope.

    Results of calls are cached in the `Fatgraph` instance on which
    the method is called, so they are automatically dropped when that
    object is collected.
    """
    @functools.wraps(func)
    def wrapper(obj, edgeno):
        try:
            cache = obj._cache_contract
        except AttributeError:
            obj._cache_contract = cache = weakref.WeakValueDictionary()
        try:
            return cache[edgeno]
        except KeyError:
            result = func(obj, edgeno)
            cache[edgeno] = result
            return result
    return wrapper


#@cython.ccall
def ocache_eq(func):
    """Cache result of the 2-ary symmetric `Fatgraph.__eq__` method.

    Results of calls are cached in the `Fatgraph` instance on which
    the method is called, *and* (simmetrically) in the `Fatgraph` to
    which that is compared to.
    
    Only a weak reference to the compared-to graph is held, so the
    cache entry is automatically freed when the compared-to Fatgraph
    goes out of scope.
    """
    @functools.wraps(func)
    def wrapper(o1, o2):
        try:
            cache1 = o1._cache_eq
        except AttributeError:
            o1._cache_eq = cache1 = weakref.WeakKeyDictionary()
        try:
            cache2 = o2._cache
        except AttributeError:
            o2._cache_eq = cache2 = weakref.WeakKeyDictionary()
        try:
            return cache1[o2]
        except KeyError:
            result = func(o1, o2)
            cache1[o2] = result
            cache2[o1] = result
            return result
    return wrapper


#@cython.ccall
def ocache_isomorphisms(func):
    """Cache results of `isomorphisms(g1,g2)` methods, which return an
    iterator/generator.

    Iterator results cannot be cached like any other object, because
    they need to return the same set of values each time the
    generating function is invoked.

    Results of calls are cached in the `Fatgraph` instance on which
    the method is called, so they are automatically dropped when that
    object is collected.

    Only a weak reference to the target Fatgraph is held, so the cache
    entry is automatically freed when it goes out of scope.
    """
    @functools.wraps(func)
    def wrapper(o1, o2):
        try:
            cache = o1._cache_isomorphisms
        except AttributeError:
            cache = o1._cache_isomorphisms = weakref.WeakKeyDictionary()
        try:
            return cache[o2].replay()
        except KeyError:
            result = _IteratorRecorder(func(o1, o2))
            cache[o2] = result
            return result.replay()
    return wrapper



## main: run tests

if "__main__" == __name__:
    import doctest
    doctest.testmod(name="cache",
                    optionflags=doctest.NORMALIZE_WHITESPACE)
