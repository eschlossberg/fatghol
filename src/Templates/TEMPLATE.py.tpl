#! /usr/bin/env python
#
"""(>>>POINT<<<)
"""
__docformat__ = 'reStructuredText'




## main: run tests

if "__main__" == __name__:
    import doctest
    doctest.testmod(name="(>>>FILE_SANS<<<)",
                    optionflags=doctest.NORMALIZE_WHITESPACE)
