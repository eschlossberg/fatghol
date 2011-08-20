#! /usr/bin/env python
#
"""Utilities for caching computation results.
"""
__docformat__ = 'reStructuredText'


## additional imports

from decorator import decorator
from time import time
import weakref

## local imports

from iterators import Iterator
from persist import PersistedIterator


## auxiliary classes

class _TimeBasedUnique(Iterator):
    """Return a new unique ID at each iteration step.

    The returned ID is the number of nanoseconds elapsed since the
    system epoch.  The returned ID is guaranteed to be monotonically
    increasing::

      >>> u = TimeBasedUnique()
      >>> u.next() < u.next() < u.next()
      True

    """
    def __init__(self):
        pass

    def next(self):
        return int(time() * 1000000)

_unique = _TimeBasedUnique()


class _IteratorRecorder(object):
    """Cache results from a given iterator.  Client classes provided
    by the `replay()` method will then replay the iterator history;
    multiple players can replay from one recorded source
    independently.

    **WARNING:** This is not thread-safe!
    """

    __slots__ = ['iterable', 'history']
    
    def __init__(self, iterable):
        self.iterable = iterable
        self.history = []
        #Cacheable.__init__(self)

    def __iter__(self):
        return self.replay()

    def advance(self):
        """Record next item from the source iterator."""
        self.history.append(self.iterable.next())

    def replay(self):
        """Return a new player."""
        return _IteratorReplayer(self)
    

class _IteratorReplayer(Iterator):
    """Replay values recorded into a given `_IteratorRecorder` class;
    multiple players can replay from one recorded source
    independently.

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
            raise


class Cacheable(object):
    """Instances of this class provide an interface for use by caches
    based on a `dict` subclass.
    """

    __slots__ = ['__id', '_cache']

    def __init__(self):
        self.__id = _unique.next()

    def cache_id(self):
        """Return a integer value which is unique and guaranteed not
        to be re-used.
        
        (In contrast, Python's standard `id()` function returns a
        value based on the intance memory address, which *can* be
        re-used during the lifetime of a program.)
        """
        return self.__id


def cache_id(o):
    try:
        return o.cache_id()
    except AttributeError:
        return id(o)



## caching functions

@decorator
def fcache(func, *args):
    """Cache result of a generic function.

    This decorator can cache results of calls `func(*args)`.

    CAVEAT: The result cache is held in the function object itself, so
    any object present in the cache will *not* be garbage collected!
    """
    try:
        rcache = func._cache
    except AttributeError:
        func._cache = rcache = {}
    key = (args,)
    try:
        return rcache[key]
    except KeyError:
        result = func(*args)
        rcache[key] = result
        return result


@decorator
def fcache_iterator(func, *args):
    """Cache results of a function returning an iterator/generator.

    Iterator results cannot be cached like any other object, because
    they need to return the same set of values each time the
    generating function is invoked.
    """
    try:
        rcache = func._cache
    except AttributeError:
        rcache = func._cache = {}
    key = (func.func_name, args)
    try:
        return rcache[key].replay()
    except KeyError:
        result = _IteratorRecorder(func(*args))
        rcache[key] = result
        return result.replay()



@decorator
def ocache0(func, obj):
    """Cache result of a nullary object method.

    This decorator can cache results of calls `obj.method()`. The
    result cache is held in the object itself; therefore, to cache
    result from methods of objects using a '__slots__' declaration, a
    '_cache' slot must be present and writable.
    """
    try:
        rcache = obj._cache
    except AttributeError:
        obj._cache = rcache = {}
    key = func.func_name
    try:
        return rcache[key]
    except KeyError:
        result = func(obj)
        rcache[key] = result
        return result


@decorator
def ocache_weakref(func, obj, *args):
    """Cache result of a generic object method.

    This decorator can cache results of calls `obj.method()`. The
    result cache is held in the object itself; therefore, to cache
    result from methods of objects using a '__slots__' declaration, a
    '_cache' slot must be present and writable.
    """
    try:
        rcache = obj._cache
    except AttributeError:
        obj._cache = rcache = {}
    key = (func.func_name, args)
    try:
        result = rcache[key]()
    except KeyError:
        result = func(obj, *args)
        def cleanup(ref): del rcache[key]
        rcache[key] = weakref.ref(result, cleanup)
    return result


@decorator
def ocache_symmetric(func, o1, o2):
    """Cache result of 2-ary symmetric method.

    This decorator can cache results of calls `obj1.method(obj2)`,
    subject to the constraint that swapping `obj1` and `obj2` gives
    the *same result* (e.g., equality tests).  When
    `obj1.method(obj2)` is first called, the result is also stored in
    the location corresponding to `obj2.method(obj1)`.

    The result cache is held in the object itself; therefore, to cache
    result from methods of objects using a '__slots__' declaration, a
    '_cache' slot must be present and writable.
    """
    try:
        rcache1 = o1._cache
    except AttributeError:
        o1._cache = rcache1 = {}
    try:
        rcache2 = o2._cache
    except AttributeError:
        o2._cache = rcache2 = {}
    key = (func.func_name, cache_id(o2))
    try:
        result = rcache1[key]
    except KeyError:
        result = func(o1, o2)
        rcache1[key] = result
        rcache2[(func.func_name, cache_id(o1))] = result
    return result


@decorator
def ocache_iterator(func, o1, o2):
    """Cache results of `isomorphism(g1,g2)` methods, which return an
    iterator/generator.

    Iterator results cannot be cached like any other object, because
    they need to return the same set of values each time the
    generating function is invoked.
    """
    try:
        rcache = o1._cache
    except AttributeError:
        rcache = o1._cache = {}
    key = (func.func_name, o2.cache_id())
    try:
        return rcache[key].replay()
    except KeyError:
        result = _IteratorRecorder(func(o1, o2))
        rcache[key] = result
        return result.replay()



## main: run tests

if "__main__" == __name__:
    import doctest
    doctest.testmod(name="cache",
                    optionflags=doctest.NORMALIZE_WHITESPACE)
