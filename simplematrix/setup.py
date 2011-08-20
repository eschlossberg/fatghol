#!/usr/bin/env python

"""
setup.py file for SWIG `simplematrix` interface
"""

from distutils.core import (
    setup,
    Extension
    )


simplematrix_module = Extension('_simplematrix',
                           sources=['simplematrix_wrap.cxx'],
                           )

setup (name = 'simplematrix',
       version = '0.1',
       author      = "Riccardo Murri <riccardo.murri@gmail.com>",
       description = """Simple swig example from docs""",
       ext_modules = [simplematrix_module],
       py_modules = ["simplematrix"],
       )
