#! /usr/bin/env python
#
"""Utilities for caching computation results.
"""
__docformat__ = 'reStructuredText'


## stdlib imports
from collections import defaultdict
import functools
from time import time
import weakref


## additional imports

import cython

## local imports

from iterators import Iterator


## auxiliary classes

@cython.cclass
class _TimeBasedUnique(Iterator):
    """Return a new unique ID at each iteration step.

    The returned ID is the number of nanoseconds elapsed since the
    system epoch.  The returned ID is guaranteed to be monotonically
    increasing::

      >>> u = _TimeBasedUnique()
      >>> u1 = u.next()
      >>> u2 = u.next()
      >>> u3 = u.next()
      >>> u1 < u2 < u3
      True

    """
    def __init__(self):
        pass

    def next(self):
        return int(time() * 1000000)

_unique = _TimeBasedUnique()


@cython.cclass
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
    

@cython.cclass
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


@cython.cclass
class Caching(object):
    """Instances of this class provide an interface for use by caches
    based on a `dict` subclass.
    """

    __slots__ = [
        '__id',
        '_cache0', # cache for nullary methods
        '_cache',  # strong result cache
        '_wcache', # WeakValueRef cache
        ]

    def __init__(self):
        self.__id = _unique.next()

    @cython.ccall
    def cache_id(self):
        """Return a integer value which is unique and guaranteed not
        to be re-used.
        
        (In contrast, Python's standard `id()` function returns a
        value based on the instance memory address, which *can* be
        re-used during the lifetime of a program.)
        """
        return self.__id


@cython.ccall
def cache_id(o):
    try:
        return o.cache_id()
    except AttributeError:
        return id(o)



## caching functions

# store 
_func_cache = defaultdict(dict)

@cython.ccall
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


@cython.ccall
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


@cython.ccall
def ocache_weakref(func):
    """Cache result of a generic object method.
    
    Only a weak reference to the method's result is held, so the cache
    entry is automatically freed when the result object goes out of
    scope.

    On the contrary, function arguments are stored as strong
    references, so any object referenced in the arguments will be kept
    alive as long as it is in the cache.

    This decorator can cache results of calls `obj.method()`. The
    result cache is held in the object itself; therefore, to cache
    result from methods of objects using a '__slots__' declaration, a
    '_cache' slot must be present and writable.
    """
    @functools.wraps(func)
    def wrapper(obj, *args):
        try:
            cache = obj._wcache
        except AttributeError:
            obj._wcache = cache = defaultdict(weakref.WeakValueDictionary)
        try:
            return cache[func.func_name][args]
        except KeyError:
            result = func(obj, *args)
            cache[func.func_name][args] = result
            return result
    return wrapper


@cython.ccall
def ocache_symmetric(func):
    """Cache result of 2-ary symmetric method.

    This decorator can cache results of calls `obj1.method(obj2)`,
    subject to the constraint that swapping `obj1` and `obj2` gives
    the *same result* (e.g., equality tests).  When
    `obj1.method(obj2)` is first called, the result is also stored in
    the location corresponding to `obj2.method(obj1)`.    
    """
    @functools.wraps(func)
    def wrapper(o1, o2):
        try:
            cache1 = o1._cache
        except AttributeError:
            o1._cache = cache1 = defaultdict(dict)
        try:
            cache2 = o2._cache
        except AttributeError:
            o2._cache = cache2 = defaultdict(dict)
        try:
            return cache1[func.func_name][cache_id(o2)]
        except KeyError:
            result = func(o1, o2)
            cache1[func.func_name][cache_id(o2)] = result
            cache2[func.func_name][cache_id(o1)] = result
            return result
    return wrapper


@cython.ccall
def ocache_iterator(func):
    """Cache results of `isomorphism(g1,g2)` methods, which return an
    iterator/generator.

    Iterator results cannot be cached like any other object, because
    they need to return the same set of values each time the
    generating function is invoked.
    """
    @functools.wraps(func)
    def wrapper(o1, o2):
        try:
            cache = o1._wcache
        except AttributeError:
            cache = o1._wcache = defaultdict(weakref.WeakValueDictionary)
        try:
            return cache[func.func_name][cache_id(o2)].replay()
        except KeyError:
            result = _IteratorRecorder(func(o1, o2))
            cache[func.func_name][cache_id(o2)] = result
            return result.replay()
    return wrapper



## main: run tests

if "__main__" == __name__:
    import doctest
    doctest.testmod(name="cache",
                    optionflags=doctest.NORMALIZE_WHITESPACE)
