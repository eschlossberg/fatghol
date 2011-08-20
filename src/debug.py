#! /usr/bin/env python
#
"""Miscellaneous support functions for `assert` statements and
run-time debugging.
"""
__docformat__ = 'reStructuredText'


import operator
import sys
import types


from decorator import decorator


def is_sequence_of_type(t, seq):
    """Return `True` if all items of sequence `s` are of type `t`.

    Examples::
      >>> is_sequence_of_type(types.IntType, [1,2,3])
      True
      >>> is_sequence_of_type(types.IntType, [1,2,"xxx"])
      False
      >>> is_sequence_of_type(types.StringType, ["xxx","yyy"])
      True
      >>> is_sequence_of_type(types.IntType, [0])
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


@decorator
def trace(func, *args, **kwargs):
    incantation = "%s(%s, %s)" % (func.func_name,
                                  str.join(", ",
                                           ["<0x%x>:%s" % (id(a), repr(a))
                                            for a in args]),
                                  str.join(", ",
                                           ["%s=%s," % (k,v)
                                            for k,v in kwargs.iteritems()]))
    print "DEBUG: Entering %s" % incantation
    result = func(*args)
    print "DEBUG: Got result %s from %s" % (result, incantation)
    return result


# from: http://www.friday.com/bbum/2005/10/27/tracing-python-execution/
import sys
import linecache
import inspect

def _traceit(frame, event, arg):
    if event == 'line':
        lineno = frame.f_lineno
        if '__file__' in frame.f_globals:
            filename = frame.f_globals['__file__']
            if (filename.endswith('.pyc') or
                filename.endswith('.pyo')):
                filename = filename[:-1]
            name = frame.f_globals['__name__']
            line = linecache.getline(filename, lineno)
        else:
            name = '[unknown]'
            try:
                src = inspect.getsourcelines(frame)
                line = src[lineno]
            except IOError:
                line = 'Unknown code named [%s].  VM instruction #%d' % \
                    (frame.f_code.co_name, frame.f_lasti)
        print '%s:%s: %s' % (name, lineno, line.rstrip())
    return _traceit

def start_tracing():
    sys.settracing(_traceit)

def stop_tracing():
    sys.settracing(lambda x,y,z: None)
    


## main: run tests

if "__main__" == __name__:
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
