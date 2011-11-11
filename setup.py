#! /usr/bin/env python

import os
swdir = os.getcwd()+'/sw'

from distutils.core import setup
from distutils.extension import Extension


ext_commands = { }
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
# if Cython is available, use it to compile modules
try:
    import Cython.Distutils
    ext_commands['build_ext'] = Cython.Distutils.build_ext
    NO_ASSERTS = [
        ('PYREX_WITHOUT_ASSERTIONS', 1),
        ('NDEBUG', 1),
        ]
    ext_modules.extend([
        Extension("aggregate",      ["aggregate.py"],      define_macros=NO_ASSERTS),
        Extension("cache",          ["cache.py"],          define_macros=NO_ASSERTS),
        Extension("combinatorics",  ["combinatorics.py"],  define_macros=NO_ASSERTS),
        Extension("cycliclist",     ["cycliclist.py"],     define_macros=NO_ASSERTS),
        Extension("graph_homology", ["graph_homology.py"], define_macros=NO_ASSERTS),
        Extension("homology",       ["homology.py"],       define_macros=NO_ASSERTS),
        Extension("iterators",      ["iterators.py"],      define_macros=NO_ASSERTS),
        Extension("loadsave",       ["loadsave.py"],       define_macros=NO_ASSERTS),
        Extension("rg",             ["rg.py"],             define_macros=NO_ASSERTS),
        Extension("utils",          ["utils.py"],          define_macros=NO_ASSERTS),
        Extension("valences",       ["valences.py"],       define_macros=NO_ASSERTS),
        ])
except ImportError:
    pass
    

setup (
    cmdclass = ext_commands,
    ext_modules = ext_modules,
)
