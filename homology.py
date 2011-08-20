#! /usr/bin/env python
#
"""Generic classes for homological algebra.
"""
__docformat__ = 'reStructuredText'


import types

# import the NZMATH `Rational` class
import sys; sys.path.append('./NZMATH-0.7.0')
from nzmath.rational import Rational


class VectorSpace(object):
    """Represent the vector space generated by the given `base` vectors.

    After construction, you can retrieve the base vectors set and the
    dimension from instance attributes `base` and `dimension`.
    
    The `base` elements are assumed to be *linearly independent*, so
    the `dimension` of the generated vector space equals the number of
    elements in the base set.
    """
    def __init__(self, base, typecast=Rational):
        """Constructor, taking explicit base and coefficient numeric type.

        First argument `base` is a sequence of base vectors; no
        requirement is placed on the type of base vectors.  The `base`
        object should support:
          - the `len` operator;
          - the `index` operator (with the same semantics of the `list` one)

        Optional second argument `typecast` is a factory function for
        the coefficients in the coordinate vectors.  The function
        `typecast` should accept a single integer argument and return
        an instance of a type that supports addition (Python operators
        `__add__` and `__iadd__`).
        """
        self.base = base
        self.dimension = len(base)
        self._typecast = typecast

    def __repr__(self):
        return "VectorSpace(%s)" % self.base
    def __str__(self):
        return "Vector space with base %s" % self.base
    
    def coordinates(self, element):
        """Return the coordinate vector of `element`.

        Argument `element` represents a linear combination as a list
        of pairs `(vector, coefficient)`, where `vector` is an item in
        the `base` (specified when constructing this object).
        """
        coordinates = [ self._typecast(0)
                        for i in xrange(self.dimension) ]
        for (vector, coefficient) in iter(element):
            coordinates[self.base.index(vector)] += coefficient
        return coordinates

    
class ChainComplex(object):
    """Represents a (finite-length) chain (homology) complex.

    A `ChainComplex` `C` of length `l` comprises vector spaces `C[i]`
    and differentials `C.differential[i]`; each `C[i]` represents the
    part of the graded vector space `C` having degree `i`.  The map
    `C.differential[i]` sends (linear combinations of) elements in
    vector space `C[i]` to linear combinations of vectors in `C[i-1]`;
    the `coordinates` method of `C[i-1]` will be used to obtain a
    numerical representation of the differentiated element.
    
    A `ChainComplex` instance must be initialized by assigning
    `VectorSpace` instances into each `C[i]` (for 0 <= `i` <
    `len(C)`), and appropriate maps into `C.differential[i]` (for 1 <=
    `i` < `len(C)`)::

      >>> # chain homology of a segment
      >>> C = ChainComplex(2)
      >>> C[1] = VectorSpace(['a'])
      >>> C[0] = VectorSpace(['b0', 'b1'])
      >>> C.differential[1] = lambda _: [('b0',1), ('b1', -1)]

    Indices of the slices `C[i]` run from 0 to `len(C)-1` (inclusive).
    The Python `len` operator returns the total length of the
    complex::
      
      >>> len(C)
      2

    At present, the only supported operation on chain complexes is
    computing the rank of homology groups::

      >>> C.compute_homology_ranks()
      [1, 0]
    """
    
    def __init__(self, length, modules=None, differentials=None):
        """Create a chain complex of specified length."""
        assert length > 0, \
                   "ChainComplex.__init__:"\
                   " argument `length` must be a positive integer," \
                   " but got `%s`." % length
        #: Total length of the complex.
        self.length = length
        #: Boundary operators; `differentials[i]` sends elements in
        #  `C[i]` to elements in `C[i+1]`.
        self.differential = [None]
        if differentials:
            assert len(differentials) == length-1, \
                   "ChainComplex.__init__:" \
                   " supplied `differentials` argument does not match" \
                   " supplied `length` argument."
            self.differential.extend(differentials)
        else:
            self.differential.extend([None] * length)
        #: The vector spaces supporting the differential complex.
        if modules:
            assert len(modules) == length, \
                   "ChainComplex.__init__:" \
                   " supplied `modules` argument does not match" \
                   " supplied `length` argument."
            self.module = modules
        else:
            self.module = [None] * length

    def __repr__(self):
        return "ChainComplex(%d, modules=%s, differentials=%s)" \
               % (self.length, self.module, self.differential)
    def __str__(self):
        return repr(self)

    ## list-like interface: support C[i] and len(C) syntax
    def __len__(self):
        return self.length
    def __getitem__(self, i):
        """Return the `i`-th pair (module, boundary operator)."""
        return (self.module[i], self.differential[i])
    def __setitem__(self, i, val):
        """Set the `i`-th support module and, optionally, boundary operator.

        ::
          C[i] = (module, differential)  # set `i`-th module and boundary op.
          C[i] = module                  # only set module
        """
        if (isinstance(val, tuple)):
            assert len(val) == 2, \
                   "ChainComplex.__setitem__:" \
                   " Need a 2-tuple (module, differential), but got `%s`" % val
            (self.module[i], self.differential[i]) = val
        else:
            self.module[i] = val

    def compute_homology_ranks(self):
        """Compute and return (list of) homology group ranks.

        Returns a list of integers: item at index `n` is the rank of
        the `n`-th homology group of this chain complex.  Since the
        chain complex has finite length, homology group indices can
        only run from 0 to the length of the complex (all other groups
        being, trivially, null).

        Examples::
        
          >>> # chain homology of a point
          >>> C_point = ChainComplex(1)
          >>> C_point[0] = VectorSpace(['a'])
          >>> C_point.compute_homology_ranks()
          [1]
          
          >>> # chain homology of a segment
          >>> C_segment = ChainComplex(2)
          >>> C_segment[1] = VectorSpace(['a'])
          >>> C_segment[0] = VectorSpace(['b0', 'b1'])
          >>> C_segment.differential[1] = lambda _: [('b0',1), ('b1', -1)]
          >>> C_segment.compute_homology_ranks()
          [1, 0]
          
          >>> # chain homology of a circle
          >>> C_circle = ChainComplex(2)
          >>> C_circle[1] = VectorSpace(['a'])
          >>> C_circle[0] = VectorSpace(['b'])
          >>> C_circle.differential[1] = lambda _: []
          >>> C_circle.compute_homology_ranks()
          [1, 1]
          
        """
        ## pass 1: compute boundary operators in matrix form
        ## 
        # FIXME: since we're only interested in computing the
        # rank, should we instanciate the matrix as row-major or
        # column-major, depending on which dimension is lesser?
        # (For doing Gaussian elimination on a smaller set.)

        #: Matrix form of boundary operators; the `i`-th differential
        #: is `dim C[i-1]` rows (range) by `dim C[i]` columns (domain),
        #: stored in column-major format.
        D = [ [ self.module[i-1].coordinates(self.differential[i](b))
                for b in self.module[i].base ]
              for i in xrange(1, self.length) ]
        
        ## pass 2: compute rank and nullity of boundary operators
        ##
        ## We reduce (destructively) every boundary operator matrix to
        ## column Echelon form by Gaussian elimination, computing the
        ## rank in the process.
        ##
        #: ranks of `D[n]` matrices, for 0 <= n < len(self); the differential
        #: `D[0]` is the null map.
        ranks = [ 0 ]  
        for n in xrange(len(D)):
            A = D[n]            # micro-optimization (saves a few lookups)
            columns = len(A)    #: number of columns
            rows = len(A[0])    #: number of rows
            i = 0  #: row index
            j = 0  #: column index
            rank = 0  #: computed rank of matrix `A`
            while (i < rows) and (j < columns):
              # find pivot in row i, starting at column j:
              pivot_column = j
              for jj in xrange(j+1, columns):
                if abs(A[jj][i]) > abs(A[pivot_column][i]):
                  pivot_column = jj
              if A[pivot_column][i] != 0:
                rank += 1
                # swap columns `j` and `pivot_column`
                self[pivot_column], self[j] = self[j], self[pivot_column]
                # divide each entry in column `j` by `A[j][i]`
                lead = A[j][i]
                for ii in xrange(rows):
                    A[j][ii] /= lead
                # A[-,u] -= A[-,j] * A[i,u]
                for u in xrange(j+1, columns):
                  for ii in xrange(rows):
                      A[u][ii] -= A[j][ii] * A[u][i]
                      # now A[u][i] will be 0, since:
                      # A[u][i] - A[j][i] * A[u][i] = A[u][i] - 1 * A[u][i] = 0.
                j += 1
              i += 1
            ranks.append(rank)

        ## pass 3: compute homology group ranks from rank and nullity
        ##         of boundary operators.
        ##
        ## By the rank-nullity theorem, if A:V-->W is a linear map,
        ## then null(A) =  dim(V) - rk(A), hence:
        ##   dim(Z_i) = null(D_i) = dim(C_i) - rk(D_i)
        ##   dim(B_i) = rk(D_{i+1})
        ## Therefore:
        ##   h_i = dim(H_i) = dim(Z_i / B_i) = dim(Z_i) - dim(B_i)
        ##       = dim(C_i) - rk(D_i) - rk(D_{i+1})
        ##
        ranks.append(0) # augment complex with the null map.
        return [ (self.module[i].dimension - ranks[i] - ranks[i+1])
                 for i in xrange(self.length) ]
    

        
## main: run tests

if "__main__" == __name__:
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
