#! /usr/bin/env python
#
"""
Global shared runtime state.

Implemented using the `Singleton` and `Monostate` patterns, see
http://code.activestate.com/recipes/66531/ for original code and a
useful discussion.
"""
__docformat__ = 'reStructuredText'


class Singleton(object):
    def __new__(cls, *p, **k):
        if not '_the_instance' in cls.__dict__:
            cls._the_instance = object.__new__(cls)
        return cls._the_instance


class Monostate(object):
    __state = {}
    def __new__(cls, *p, **k):
        self = object.__new__(cls, *p, **k)
        self.__dict__ = cls.__state
        return self


runtime = Monostate()


## main: run tests

if "__main__" == __name__:
    import doctest
    doctest.testmod(name="monostate",
                    optionflags=doctest.NORMALIZE_WHITESPACE)
