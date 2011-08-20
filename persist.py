#! /usr/bin/env python
#
"""A simple-minded persistent framework.
"""
__docformat__ = 'reStructuredText'


## stdlib imports

import pickle


## local imports

from iterators import Iterator



class CheckpointingIterator(Iterator):
    """Iterator class, periodically checkpointing state to persistent storage.

    Every `checkpoint_every` (initialization parameter) invocations of
    the `.next()` method, the iterator state is saved to disk.
    """

    def __init__(self, iterable, filename, checkpoint_every=10):
        self.iterable = iter(iterable)
        self.checkpoint_every = checkpoint_every
        self.__filename = filename
        self.__countdown = checkpoint_every
        self.__stop_iteration = False

    def checkpoint(self, filename=None, remember=False):
        """Save an object to the associated persistent store.
        """
        if filename is None:
            filename = self.__filename
        store = open(filename, 'w')
        self.pickle()
        pickle.dump(self, store)
        store.close()
        if remember:
            self.__filename = filename

    def next(self):
        if self.__stop_iteration:
            raise StopIteration
        try:
            next = self.iterable.next()
            self.__countdown -= 1
            if self.__countdown == 0:
                # checkpoint to file
                self.checkpoint()
                self.__countdown = self.checkpoint_every
            return next
        except StopIteration:
            self.__stop_iteration = True
            self.checkpoint()
            raise StopIteration

    def pickle(self):
        """Called *before* the iterator is pickled."""
        try:
            self.iterable.pickle()
        except AttributeError:
            pass

    def unpickle(self):
        """Called *after* the iterator has been un-pickled."""
        try:
            self.iterable.unpickle()
        except AttributeError:
            pass

    
class RecordingIterator(CheckpointingIterator):
    """Records the results of another iterator, and lets caller rewind
    and replay values.

    Example::
      >>> r = RecordingIterator([1,2,3])
      >>> for x in r: print x
      1
      2
      3
      >>> r.replay(1)
      >>> for x in r: print x
      2
      3
      
    """

    def __init__(self, iterable, filename, checkpoint_every=10, restart=True):
        self.pos = -1
        self.history = []
        self.restart = restart
        CheckpointingIterator.__init__(self, iterable,
                                       filename, checkpoint_every)

    def next(self):
        if (self.pos < len(self.history) - 1):
            self.pos += 1
            return self.history[self.pos]
        else:
            next = super(RecordingIterator, self).next()
            self.history.append(next)
            self.pos += 1
            return next

    def replay(self, pos=0):
        """Replay iterator starting at the given position.
        By default, replays iterator from start of the recorded values.
        """
        self.pos = pos-1

    def unpickle(self):
        """Called *after* the iterator has been un-pickled."""
        try:
            self.iterable.unpickle()
        except AttributeError:
            # replay iterator at start
            if self.restart:
                self.replay()


def PersistedIterator(factory, checkpoint_every=10, record=True, replay=True):
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
      ...   def __iter__(self):
      ...     return self
      ...   def next(self):
      ...     return self.arg
      >>> X = PersistedIterator(_X)
      >>> x = X(1)
      >>> x.next()
      1
      >>> x.arg =2
      >>> checkpoint(x)
      >>> x = X(1) # uses saved copy
      >>> x.next()
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
                + '.cache'
        else:
            filename = factory.__name__ + ".cache"
        try:
            # try loading it
            store = open(filename, 'r')
            instance = pickle.load(store)
            store.close()
            # tell instance it has been un-pickled
            instance.unpickle()
        except IOError: # filename does not exist
            # run normal class initialization
            if record:
                instance = RecordingIterator(factory(*args, **kwargs),
                                             filename, checkpoint_every,
                                             restart=replay)
            else:
                instance = CheckpointingIterator(factory(*args, **kwargs),
                                                 filename, checkpoint_every)

        return instance

    return make


## main: run tests

if "__main__" == __name__:
    import doctest
    doctest.testmod(name="persist",
                    optionflags=doctest.NORMALIZE_WHITESPACE)
