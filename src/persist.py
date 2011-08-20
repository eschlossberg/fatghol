#! /usr/bin/env python
#
"""A simple-minded persistent framework.
"""
__docformat__ = 'reStructuredText'


## stdlib imports

import pickle


## local imports

from utils import Iterator



#: object-to-filename dict
_backing_store = dict()



class RecordingIterator(Iterator):
    """Records the results of another iterator, and lets caller rewind
    and replay values.

    Example::
      >>> r = RecordingIterator([1,2,3])
      >>> for x in r: print x
      1
      2
      3
      >>> r.rewind(1)
      >>> for x in r: print x
      2
      3
      
    """

    def __init__(self, iterable, checkpoint=10):
        self.iterable = iter(iterable)
        self.pos = -1
        self.history = []
        self.checkpoint = checkpoint

    def next(self):
        if (self.pos < len(self.history) - 1):
            self.pos += 1
            return self.history[self.pos]
        else:
            next = self.iterable.next()
            self.history.append(next)
            self.pos += 1
            if len(self.history) % self.checkpoint == 0:
                # checkpoint to file
                try:
                    checkpoint(self)
                except AttributeError:
                    pass
            return next

    def rewind(self, pos=0):
        """Rewind iterator at the given position.
        By default, rewinds iterator at start of the recorded values.
        """
        self.pos = pos-1

    def thaw(self):
        """Called after the iterator has been un-pickled."""
        try:
            self.iterable.thaw()
        except AttributeError:
            # rewind iterator at start
            self.rewind()


def PersistedIterator(factory):
    """Wrap a factory function: the wrapped version looks for and uses
    the saved copy of an object, and falls back to constructing a new
    one if no saved copy exists.

    The file name of the saved copy is determined by the factory
    function name, with the string representation of the invocation
    arguments appended.

    Example::
      >>> class _X(object):
      ...   def __init__(self, arg):
      ...     self.arg = arg
      >>> X = cached(_X)
      >>> x = X(1)
      >>> x.arg
      1
      >>> x.arg =2
      >>> checkpoint(x)
      >>> x = X(1) # uses saved copy
      >>> x.arg
      2

    *Note:* If a class constructor is wrapped, like in the example above,
    the wrapped constructor *must* be bound to a different Python name,
    or `pickle` won't be able to load the saved object.
    """

    def make(*args, **kwargs):
        """
        If there is a saved state file matching the given init args,
        the contents of that file are loaded and the pickled object is
        returned; otherwise, normal class initialization is run.
        """
        # compute backing store name from init args
        if len(args) > 0:
            filename = factory.__name__ \
                + "-" + str.join(",", [str(x) for x in args]) \
                + '.pickle'
        else:
            filename = factory.__name__ + ".pickle"
        try:
            # try loading it
            store = open(filename, 'r')
            instance = pickle.load(store)
            store.close()
            # tell instance it has been un-pickled
            instance.thaw()
        except IOError: # filename does not exist
            # run normal class initialization
            instance = RecordingIterator(factory(*args, **kwargs))

        # save persistent store ref into instance
        _backing_store[id(instance)] = filename

        return instance

    return make


def checkpoint(self, filename=None):
    """Save an object to the associated persistent store.
    """
    if filename is None:
        filename = _backing_store[id(self)]
    store = open(filename, 'w')
    pickle.dump(self, store)
    _backing_store[id(self)] = filename
    store.close()



## main: run tests

if "__main__" == __name__:
    import doctest
    doctest.testmod(name="persist",
                    optionflags=doctest.NORMALIZE_WHITESPACE)
