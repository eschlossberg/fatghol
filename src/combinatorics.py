#! /usr/bin/env python
#
"""Common combinatorics constructions.
"""
__docformat__ = 'reStructuredText'


import itertools


class SetProductIterator:
    """Iterate over all elements in a cartesian product.

    Argument `factors` is a sequence, all whose items are sequences
    themselves: the returned iterator will return -upon each
    successive invocation- a list `[t_1, t_2, ..., t_n]` where `t_k`
    is an item in the `k`-th sequence.

    Examples::
      >>> list(SetProductIterator([]))
      [[]]
      >>> list(SetProductIterator([[1]]))
      [[1]]
      >>> list(SetProductIterator([[1],[1]]))
      [[1, 1]]
      >>> list(SetProductIterator([[1,2],[]]))
      [[]]
      >>> list(SetProductIterator([[1,2],[1]]))
      [[1, 1], [2, 1]]
      >>> list(SetProductIterator([[1,2],[1,2]]))
      [[1, 1], [2, 1], [1, 2], [2, 2]]
    """
    def __init__(self, factors):
        self.__closed = False
        self.__factors = factors
        self.__L = len(factors)
        self.__M = [ len(s)-1 for s in factors ]
        self.__m = [0] * self.__L
        self.__i = 0

    def __iter__(self):
        return self

    def next(self):
        if self.__closed:
            raise StopIteration
        if (0 == self.__L) or (-1 in self.__M):
            # there are no factors, or one of them has no elements
            self.__closed = True
            return []
        else:
            if self.__i < self.__L:
                # will return element corresponding to current multi-index
                result = [ s[self.__m[i]]
                           for (i,s) in enumerate(self.__factors) ]
                # advance multi-index
                i = 0
                while (i < self.__L):
                    if self.__m[i] == self.__M[i]:
                        self.__m[i] = 0
                        i += 1
                    else:
                        self.__m[i] += 1
                        break
                self.__i = i
                # back to caller
                return result
            else:
                # at end of iteration
                self.__closed = True
                raise StopIteration


class PermutationList:
    """Simulate a (read-only) list of all permutations of a prescribed order.

    The list is ordered lexicographically; at position `0` and `n!`
    are the identity permutation and the order-reversal permutation,
    respectively.
    
    The code is a port of the one described in Wikipedia at:
      http://en.wikipedia.org/wiki/Permutation#Algorithm_to_generate_permutations

    Examples::
      >>> ps = PermutationList(3)
      >>> ps[0]
      [0, 1, 2]
      >>> ps[5]
      [2, 1, 0]
    """
    def __init__(self, order):
        self.__order = order
        # pre-compute factorial
        self.factorial = [1]
        for j in xrange(0, self.__order):
            self.factorial.append((j+1)*self.factorial[-1])
    def __getitem__(self, i):
        """Return permutation at `i`-th place."""
        if i >= self.factorial[-1]:
            raise StopIteration
        order = self.__order
        perm = [ n for n in xrange(0, order) ]
        for k in xrange(0, i+1):
            for j in xrange(1, order):
                jj = j - ((k / self.factorial[j]) % (j+1))
                perm[j], perm[jj] = perm[jj], perm[j]
        return perm
    def __len__(self):
        return self.factorial[-1]


class PermutationIterator:
    """Iterate over all permutations of a given sequence.

    The code is a port of the one described in Wikipedia at:
      http://en.wikipedia.org/wiki/Permutation#Algorithm_to_generate_permutations

    Examples::
      >>> p = PermutationIterator([1,2,3])
      >>> p.next()
      [1, 2, 3]
      >>> p.next()
      [2, 1, 3]
      >>> p.next()
      [2, 3, 1]
      >>> p.next()
      [3, 1, 2]
      >>> p.next()
      [2, 1, 3]
      >>> p.next()
      [3, 2, 1]
    """
    def __init__(self, seq, initial=0):
        self.seq = seq
        self.rank = initial
        # pre-compute factorial
        self.factorial = [1]
        for j in xrange(0, len(seq)):
            self.factorial.append((j+1)*self.factorial[-1])
    def __iter__(self):
        return self
    def next(self):
        """Return next permutation of initially given `sequence`."""
        if self.rank >= self.factorial[-1]:
            raise StopIteration
        def swap(seq, pos1, pos2):
            """Swap items at positions `pos1` and `pos2` in `seq`."""
            seq[pos1],seq[pos2] = seq[pos2],seq[pos1]
        for j in xrange(1, len(self.seq)):
            swap(self.seq, j, j - ((self.rank / self.factorial[j]) % (j+1)))
        self.rank += 1
        return self.seq


class InplacePermutationIterator:
    """Iterate (destructively) over all permutations of a given sequence.

    The given sequence `seq` is altered as new permutations are
    requested through the `next()` method.  In order to use the
    returned value after subsequent calls to
    `InplacePermutationIterator.next()`, you must make a copy of it.

    The code is a port of the C++ STL one, explained in:
      http://marknelson.us/2002/03/01/next-permutation

    Examples::
      >>> [ x[:] for x in InplacePermutationIterator([0]) ]
      [[0]]
      >>> [ x[:] for x in InplacePermutationIterator([0,1]) ]
      [[1, 0], [0, 1]]
      >>> [ x[:] for x in InplacePermutationIterator([0,1,2])]
      [[0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0], [0, 1, 2]]
    """
    def __init__(self, seq, start=0, end=None):
        self.seq = seq
        self.start = start
        if end is None:
            end = len(self.seq)
        self.end = end
        if (start == end):
            self.enumeration_finished = True
        else:
            self.enumeration_finished = False
    def __iter__(self):
        return self
    def next(self):
        """Return next permutation of initially given `sequence`."""
        if self.enumeration_finished:
            raise StopIteration
        i = self.end - 1
        while True:
            if (i == self.start):
                self.seq.reverse()
                self.enumeration_finished = True
                return self.seq
            j = i
            i -= 1
            if self.seq[i] < self.seq[j]:
                k = self.end-1
                while self.seq[i] >= self.seq[k]:
                    k -= 1
                # swap seq[i] and seq[k]
                self.seq[i],self.seq[k] = self.seq[k],self.seq[i]
                # reverse slice seq[j:] *in-place*
                if self.end-2 == j:
                    self.seq[j],self.seq[-1] = self.seq[-1],self.seq[j]
                else:
                    for l in xrange(0, (self.end - 1 - j) / 2):
                        self.seq[j+l],self.seq[-1-l] = self.seq[-1-l],self.seq[j+l]
                return self.seq


def sign_of_permutation(p):
    """Return sign of permutation `p`.

    A permutation of N elements is represented as a linear
    list `p`, whose items are the numbers 0,1,...,N-1: `p`
    maps `i` to `p[i]`.

    Examples::

      >>> sign_of_permutation([0,1,2])
      1
      >>> sign_of_permutation([0,2,1])
      -1
      >>> sign_of_permutation([2,0,1])
      1

    The trivial permutations mapping a single element into
    itself and the empty permutation are assigned sign +1::

      >>> sign_of_permutation([0])
      1

      >>> sign_of_permutation([])
      1

    This is an adaptation of the `perm_sign` code by John Burkardt
    (see it among the collection at
    http://orion.math.iastate.edu/burkardt/f_src/subset/subset.f90
    or http://www.scs.fsu.edu/~burkardt/math2071/perm_sign.m ); it
    computes the sign by counting the number of interchanges
    required to change the given permutation into the identity
    one.

    """
    n = len(p)
    s = +1
    # get elements back in their home positions
    for j in xrange(n):
        q = p[j]
        if q !=j :
            p[j],p[q] = p[q],q # interchange p[j] and p[p[j]]
            s = -s             # and account for the interchange
    # note that q is now in its home position
    # whether or not an interchange was required
    return s


class FixedLengthPartitionIterator:
    """Iterate over partitions of integer `N` into *exactly* `K`
    positive integers.

    Each returned partition is a list of positive integers in
    descending order, such that their sum is `N`.

    Arguments `min_` and `max_` bound the values in each partition.

    Examples::
      >>> list(FixedLengthPartitionIterator(3,1))
      [(3,)]
      >>> list(FixedLengthPartitionIterator(3,2))
      [(2, 1)]
      >>> list(FixedLengthPartitionIterator(3,3))
      [(1, 1, 1)]
      >>> list(FixedLengthPartitionIterator(6,2))
      [(5, 1), (4, 2), (3, 3)]
      >>> list(FixedLengthPartitionIterator(8,3))
      [(6, 1, 1), (5, 2, 1), (4, 3, 1), (4, 2, 2), (3, 3, 2)]
      >>> list(FixedLengthPartitionIterator(8,4,2))
      [(2, 2, 2, 2)]
      >>> list(FixedLengthPartitionIterator(8,3,2))
      [(4, 2, 2), (3, 3, 2)]
      >>> list(FixedLengthPartitionIterator(8,3,2,3))
      [(3, 3, 2)]      
    """
    def __init__(self, N, K, min_=1, max_=None):
        # `max_` really defaults to N
        if max_ is None:
            max_ = N

        # rule out trivial cases
        if (K*min_ > N) or (K*max_ < N):
            self.done = True
            return

        self._N = N  #: integer to be partitioned
        self._K = K  #: maximum number of nonzero parts
        self._k = 0  #: current number of nonzero parts
        self._min = min_  #: minimum value of each part
        self._max = max_  #: maximum value of each part
        self._p = [0] * K #: current partition
        self.done = False #: when `True`, enumeration is over

    def __iter__(self):
        return self

    def next(self):
        if self.done:
            raise StopIteration
        if self._N == self._K * self._min:
            self.done = True
            return tuple([self._min]*self._K)
        else:
            while self._k <= self._K:
                i = self._k - 1
                while i > 0:
                    if (self._p[i]+1 <= self._p[i-1]-1) \
                           and (self._p[i-1]-1 >= self._min) \
                           and (self._p[i]+1 <= self._max):
                        self._p[i-1] -= 1
                        self._p[i] += 1
                        return tuple(self._p)
                    else:
                        i -= 1
                # only change the first `k` parts
                self._k += 1
                head = (self._N - self._k - self._min*self._K + self._min + 1)
                if (head < self._min+1) or (head > self._max) \
                       or (self._k > self._K):
                    continue
                # advance to next partition:
                # [N-2*(k-1)-(K-k), 2, ..., 2, (k-1 times) 1, ..., 1 (K-k times)]
                self._p = [head] \
                          + ([self._min + 1] * (self._k - 1)) \
                          + ([self._min] * (self._K - self._k))
                return tuple(self._p)
            raise StopIteration


def PartitionIterator(N, K, min_=1, max_=None):
    """Iterate over partitions of integer `N` into *at most* `K`
    positive integers.

    Each returned partition is a list of positive integers in
    descending order, such that their sum is `N`.

    Optional arguments `min_` and `max_` bound the values in each
    partition.

    Examples::
      >>> list(PartitionIterator(2,1))
      [(2,)]
      >>> list(PartitionIterator(3,3))
      [(3,), (2, 1), (1, 1, 1)]
      >>> list(PartitionIterator(8,3,2))
      [(8,), (6, 2), (5, 3), (4, 4), (4, 2, 2), (3, 3, 2)]
      >>> list(PartitionIterator(8,3,2,3))
      [(3, 3, 2)]      
    """
    return itertools.chain(*[FixedLengthPartitionIterator(N,k,min_,max_)
                             for k in xrange(1,K+1)])


## main: run tests

if "__main__" == __name__:
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
