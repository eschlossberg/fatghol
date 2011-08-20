#! /usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension

setup (
    ext_modules = [
        Extension("_simplematrix",
                  ["simplematrix.i"],
                  libraries=['givaro', 'gmp', 'linbox'],
                  swig_opts=['-c++', '-modern'],
                  )
        ]
)
