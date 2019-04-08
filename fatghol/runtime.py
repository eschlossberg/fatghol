#! /usr/bin/env python
#
"""
Global shared runtime state.

Implemented using the `Singleton` and `Monostate` patterns, see
http://code.activestate.com/recipes/66531/ for original code and a
useful discussion.
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
