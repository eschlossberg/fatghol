#! /usr/bin/env python
#
"""
Utilities for timing sections of code.
"""
__docformat__ = 'reStructuredText'


import time


_times = { }
_started = { }


def start(label):
    """Start a clock.

    :param str label: Arbitrary string identifying this clock for
    later reference.
    """
    _started[label] = time.clock()


def stop(label):
    """
    Stop a clock, given its reference label.

    :param str label: Arbitrary string identifying this clock;
    must match one of the arguments given to `timing.start`.

    :raise KeyError: If `label` was never passed to `timing.start`.
    
    """
    _times[label] = time.clock() - _started[label]
    del _started[label]


def get(label):
    """Return the time elapsed between the `start` and `stop` calls to
    the clock named `label`, as single floating-point number
    expressing the duration in seconds.

    :raise KeyError: If `label` was never passed to `timing.stop`.

    """
    return _times[label]
    #return (delta.day*86400 + delta.hour*3600 + delta.minute*60 
    #        + delta.second + (delta.microsecond / 1000000.0))


## main: run tests

if "__main__" == __name__:
    import doctest
    doctest.testmod(name="timing",
                    optionflags=doctest.NORMALIZE_WHITESPACE)
