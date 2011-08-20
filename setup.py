#! /usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension

import os
swdir = os.getcwd()+'/sw'

setup (
    ext_modules = [
        Extension("_simplematrix", ["simplematrix.i"],
                  swig_opts=['-c++', '-modern'],
                  include_dirs=[swdir+'/include'],
                  libraries=['linbox', 'gmp', 'givaro', 'lapack', 'cblas'],
                  library_dirs=[swdir+'/lib'],
                  define_macros=[
                      # other options include:
                      #  FATGHOL_USE_LINBOX_ELIMINATION_PIVOT_NONE
                      #  FATGHOL_USE_LINBOX_DEFAULT (black-box, as of LinBox 1.1.7)
                      #  FATGHOL_USE_RHEINFALL
                      ('FATGHOL_USE_LINBOX_ELIMINATION_PIVOT_LINEAR', 1),
                      ],
                )
        ]
)
