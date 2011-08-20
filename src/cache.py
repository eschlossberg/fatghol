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


## main content

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


class Cacheable(object):
    """Instances of this class provide an interface for use by caches
    based on a `dict` subclass.

    When an object that has been cached goes out of scope, its
    `.invalidate()` method should be called, and it will notify all
    caches holding a reference to it to release the reference, so that
    Python GC can do its job.

    Conforming caches should implement a `release(key)` method, by
    which they release any reference to object identified by `key`.
    """

    def __init__(self):
        self.__id = _unique.next()
        self.__cachers = []

    def cache_id(self):
        """Return a integer value which is unique and guaranteed not
        to be re-used.
        
        (In contrast, Python's standard `id()` function returns a
        value based on the intance memory address, which *can* be
        re-used during the lifetime of a program.)
        """
        return self.__id

    def referenced(self, cache, key):
        """Add `cache` to the set of stores that hold a reference to
        this object.
        """
        self.__cachers.append((cache, key))

    def release(self):
        """Notify all caches that this object should be removed from
        their store, so that Python's GC can collect it.
        """
        for (c, k) in self.__cachers:
            try:
                del c[k]
            except KeyError:
                pass



def cache_id(o):
    try:
        return o.cache_id()
    except AttributeError:
        return id(o)



__cache1 = {}
@decorator
def cache1(func, obj):
    """Cache result of a 1-ary object method.

    This decorator can cache results of calls `obj.method()` or
    `func(obj)`; it will yield possibly incorrect results in any other
    case.

    Unlike the `memoize` decorator, it can also work with objects that
    use a `__slots__` declaration.
    """
    rcache = __cache1.setdefault(func.func_name, {})
    key = cache_id(obj)
    if key in rcache:
        return rcache[key]
    else:
        result = func(obj)
        rcache[key] = result
        if isinstance(obj, Cacheable): obj.referenced(rcache, key)
        return result


__cache2 = {}
@decorator
def cache(func, obj, *args):
    """Cache result of a generic object method.

    This decorator can cache results of calls `obj.method(*args)` or
    `func(obj, *args)`.

    Unlike the `memoize` decorator, it can also work with objects that
    use a `__slots__` declaration.
    """
    rcache = __cache2.setdefault(func.func_name, {})
    oid = cache_id(obj)
    key = (oid, args)
    if key in rcache:
        return rcache[key]
    else:
        result = func(obj, *args)
        rcache[key] = result
        if isinstance(obj, Cacheable): obj.referenced(rcache, key)
        return result


__cache_symmetric = {}
@decorator
def cache_symmetric(func, o1, o2):
    """Cache result of 2-ary symmetric method.

    This decorator can cache results of `obj1.method(obj2)` or
    `func(obj1, obj2)` invocations, subject to the constraint that
    swapping `obj1` and `obj2` gives the *same result* (e.g., equality
    tests).

    Unlike the `memoize` decorator, it can also work with objects that
    use a `__slots__` declaration.
    """
    rcache = __cache_symmetric.setdefault(func.func_name, {})
    key = frozenset([cache_id(o1),
                     cache_id(o2)])
    if key in rcache:
        return rcache[key]
    else:
        result = func(o1, o2)
        rcache[key] = result
        if isinstance(o1, Cacheable): o1.referenced(rcache, key)
        if isinstance(o2, Cacheable): o2.referenced(rcache, key)
        return result


class _IteratorRecorder(Cacheable):
    """Cache results from a given iterator.  Client classes provided
    by the `replay()` method will then replay the iterator history;
    multiple players can replay from one recorded source
    independently.

    **WARNING:** This is not thread-safe!
    """
    
    def __init__(self, iterable):
        self.iterable = iterable
        self.history = []
        Cacheable.__init__(self)

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


__cache_iterator = {}
@decorator
def cache_iterator(func, obj, *args):
    """Cache results of a method or function returning an iterator/generator.

    Iterator results cannot be cached like any other object, because
    they need to return the same set of values each time the
    generating function is invoked.
    """
    rcache = __cache_iterator.setdefault(func.func_name, {})
    oid = cache_id(obj)
    key = (oid, args)
    if key in rcache:
        return rcache[key].replay()
    else:
        result = _IteratorRecorder(func(obj, *args))
        rcache[key] = result
        try:
            if isinstance(obj, Cacheable): obj.referenced(rcache, key)
        except AttributeError:
            pass
        return result.replay()



## main: run tests

if "__main__" == __name__:
    import doctest
    doctest.testmod(name="cache",
                    optionflags=doctest.NORMALIZE_WHITESPACE)
