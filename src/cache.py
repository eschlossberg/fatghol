#! /usr/bin/env python
#
"""Utilities for caching computation results.
"""
__docformat__ = 'reStructuredText'


## additional imports

from decorator import decorator


## local imports

from iterators import Iterator
from persist import PersistedIterator


## main content

class PersistentCount(Iterator):
    """Like `itertools.count`, except reached count is saved
    to disk and resumes where last session left.
    """

    def __init__(self, filename=None, firstval=0):
        """Initialize counter.

        Parameters:
          `filename`
            The filename to store the reched count in;
            if not provided, then the class name is used.
          `firstval`
            The value count should start at, if no saved
            value is found in the provided file.
        """
        if filename is None:
            self.__filename = self.__class__.__name__ + '.txt'
        try:
            stored = open(self.__filename, 'r')
            self.counted = int(stored.read(), 16)
            stored.close()
        except (IOError, ValueError):
            self.counted = firstval - 1
        self.__file = open(self.__filename, 'w')

    def next(self):
        self.counted += 1

        self.__file.seek(0)
        self.__file.write("0x%08x\n" % self.counted)
        
        return self.counted

_unique = PersistentCount()


def persistent_id(o):
    """Return the persistent ID of object `o`.
    If `o` has no `._persistent_id` attribute, then set it and return a
    new persistent id value.

    Persistent IDs are unique and are guaranteed not to be re-used
    during the lifetime of a program.
    """
    try:
        return o._persistent_id
    except AttributeError:
        id = _unique.next()
        o._persistent_id = id
        return id


@decorator
def cache1(func, obj, fcache={}):
    """Cache result of 1-ary object method.

    This decorator can cache results of calls `obj.method()` or
    `func(obj)`; it will yield possibly incorrect results in any other
    case.

    Unlike the `memoize` decorator, it can also work with objects that
    use a `__slots__` declaration.
    """
    rcache = fcache.setdefault(func.func_name, {})
    key = persistent_id(obj)
    if key in rcache:
        return rcache[key]
    else:
        result = func(obj)
        rcache[key] = result
        return result



@decorator
def cache2(func, o1, o2, fcache={}):
    """Cache result of 2-ary symmetric method.

    This decorator can cache results of `obj1.method(obj2)` or
    `func(obj1, obj2)` invocations, subject to the constraint that
    swapping `obj1` and `obj2` gives the *same result* (e.g., equality
    tests).

    Unlike the `memoize` decorator, it can also work with objects that
    use a `__slots__` declaration.
    """
    rcache = fcache.setdefault(func.func_name, {})
    key = frozenset([persistent_id(o1),
                     persistent_id(o2)])
    if key in rcache:
        return rcache[key]
    else:
        result = func(o1, o2)
        rcache[key] = result
        return result
 


## main: run tests

if "__main__" == __name__:
    import doctest
    doctest.testmod(name="cache",
                    optionflags=doctest.NORMALIZE_WHITESPACE)
