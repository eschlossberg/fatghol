#! /usr/bin/env python
#
"""Poor man's substitute for "pickle", supporting Cython classes.
"""
__docformat__ = 'reStructuredText'


#import cython

import logging
import os
import os.path
from zlib import adler32


#@cython.locals(filename=str,
#               checksum=cython.long, result=list,
#               checkpoint_file=object, checksum_file=object,
#               error=Exception)
#@cython.cfunc(list)
def load(filename):
    from rg import Fatgraph, Vertex, BoundaryCycle
    result = list()
    checksum = 0
    try:
        with open(filename, 'r') as checkpoint_file:
            for line in checkpoint_file:
                result.append(eval(line))
                checksum = adler32(line, checksum)
        checksum &= 0xffffffff

        with open(filename+'.sum', 'r') as checksum_file:
            saved_checksum = int(checksum_file.read(), 16)
            if checksum != saved_checksum:
                logging.warning("Computed checksum of file '%s' is 0x%x,"
                                " but saved checksum is 0x%x.  Ignoring checkpoint file.",
                                filename, checksum, saved_checksum)
                return None

        # checksum ok, return to caller
        return result
    except IOError, error:
        if error.errno == 2: # No such file or directory
            return None
        else:
            raise error
        

#@cython.locals(items=list, filename=str,
#               checksum=cython.long,
#               checkpoint_file=object, checksum_file=object,
#               line=str, error=Exception)
#@cython.cfunc
def save(items, filename):
    checksum = 0
    try:
        with open(filename, 'w') as checkpoint_file:
            for item in items:
                line = "%s\n" % repr(item)
                checksum = adler32(line, checksum)
                checkpoint_file.write(line)
        with open(filename+".sum", 'w') as checksum_file:
            checksum_file.write("0x%x\n" % (checksum & 0xffffffff))
    except Exception, error:
        # remove files to avoid corrupted restart
        try:
            if os.path.exists(filename):
                os.remove(filename)
        except Exception, ex:
            # log error and ignore it: must pass original exception to caller
            logging.warning("Error removing checkpoint file '%s': %s",
                            filename, str(ex))
        try:
            if os.path.exists(filename+".sum"):
                os.remove(filename+".sum")
        except Exception, ex:
            # log error and ignore it: must pass original exception to caller
            logging.warning("Error removing checksum file '%s': %s",
                            filename+".sum", str(ex))
        raise error



## main: run tests

if "__main__" == __name__:
    import doctest
    doctest.testmod(name="loadsave",
                    optionflags=doctest.NORMALIZE_WHITESPACE)
