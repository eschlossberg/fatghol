#! /usr/bin/env python
#
"""Classes for computing graph homology.
"""
__docformat__ = 'reStructuredText'


## logging subsystem

import itertools
import logging


## application-local imports

from aggregate import AggregateList
from combinatorics import (
    bernoulli,
    factorial,
    minus_one_exp,
    Permutation,
    PermutationIterator,
    InplacePermutationIterator,
    )
from copy import copy
from homology import *
from iterators import IndexedIterator
from rg import (
    Fatgraph,
    MgnGraphsIterator,
    NumberedFatgraph,
    Vertex, # for the doctests
    )
from rational import Rational
from valences import vertex_valences_for_given_g_and_n


def FatgraphComplex(g, n):
    """Return the fatgraph complex for given genus `g` and number of
    boundary components `n`.

    This is a factory method returning a `homology.ChainComplex`
    instance, populated with the correct vector spaces and
    differentials to compute the graph homology of the space
    `M_{g,n}`.
    """
    ## Minimum number of edges is attained when there's only one
    ## vertex; so, by Euler's formula `V - E + n = 2 - 2*g`, we get:
    ## `E = 2*g + n - 1`.
    min_edges = 2*g + n - 1
    logging.debug("  Minimum number of edges: %d", min_edges)
    
    ## Maximum number of edges is reached in graphs with all vertices
    ## tri-valent, so, combining Euler's formula with `3*V = 2*E`, we
    ## get: `E = 6*g + 3*n - 6`.  These are also graphs corresponding
    ## to top-dimensional cells.
    top_dimension = 6*g + 3*n - 6
    logging.debug("  Maximum number of edges: %d", top_dimension)

    #: list of primitive graphs, graded by number of edges
    #generators = [ AggregateList() for dummy in xrange(top_dimension) ]
    C = MgnChainComplex(top_dimension)

    # gather graphs
    chi = Rational(0)
    for graph in MgnGraphsIterator(g,n):
        grade = graph.num_edges - 1
        pool = NumberedFatgraphPool(graph)
        # compute orbifold Euler characteristics (needs to include *all* graphs)
        chi += Rational(minus_one_exp(grade-min_edges)*len(pool), pool.num_automorphisms)
        # discard non-orientable graphs
        if not pool.is_orientable:
            continue
        C.module[grade].aggregate(pool)
    C.orbifold_euler_characteristics = chi
        
    for i in xrange(top_dimension):
        logging.debug("  Initialized grade %d chain module (with dimension %d)",
                     i, len(C[i]))

    ## if __debug__:
    ##     outputs = [ None ]
    ##     for i in xrange(1, top_dimension):
    ##         output = file("M%d,%d-D%d.txt" % (g,n,i), 'w')
    ##         outputs.append(output)
    ##         # write out row vectors
    ##         for (seq, G) in enumerate(C.module[i-1].itervalues()):
    ##             output.write("R %04d %s\n" % (seq, repr(G)))
    ##         # dump column vectors
    ##         for (seq, G) in enumerate(C.module[i].itervalues()):
    ##             output.write("C %04d %s\n" % (seq, repr(G)))
    ##         # dump matrix
    ##     D = C.compute_boundary_operators(outputs)
    ##     for i in xrange(1, top_dimension):
    ##         if outputs[i]:
    ##             outputs[i].close()
        
    return C



class MgnChainComplex(ChainComplex):
    """A specialized `ChainComplex`.
    """
    def __init__(self, length):
        ChainComplex.__init__(self, length)
        for i in xrange(length):
            self.module[i] = AggregateList()

    def compute_boundary_operators(self, outputs=None):
        #: Matrix form of boundary operators; the `i`-th differential
        #: `D[i]` is `dim C[i-1]` rows (range) by `dim C[i]` columns
        #: (domain).
        m = self.module # micro-optimization
        D = DifferentialComplex()
        D.append(NullMatrix, 0, len(m[0]))
        for i in xrange(1, len(self)):
            p = len(m[i-1]) # == dim C[i-1]
            q = len(m[i])   # == dim C[i]
            d = SimpleMatrix(p, q)
            j0 = 0
            for pool1 in m[i].iterblocks():
                k0 = 0
                for pool2 in m[i-1].iterblocks():
                    for l in pool1.graph.edges():
                        if pool1.graph.is_loop(l):
                            continue # with next `l`
                        for (j, k, s) in NumberedFatgraphPool.facets(pool1, l, pool2):
                            assert k < len(pool2)
                            assert j < len(pool1)
                            assert k+k0 < p
                            assert j+j0 < q
                            d.addToEntry(k+k0, j+j0, s)
                            ## if __debug__:
                            ##     if outputs:
                            ##         outputs[i].write("M %04d %04d %d\n" % (k+k0, j+j0, s))
                    k0 += len(pool2)
                j0 += len(pool1)
            D.append(d, p, q)
            logging.info("  Computed %dx%d matrix D[%d]", p, q, i)
        return D



class NumberedFatgraphPool(object):
    """An immutable virtual collection of `NumberedFatgraph`s.
    Items are all distinct (up to isomorphism) decorations of a
    `Fatgraph` instance `graph` with a numbering of the boundary
    cycles.

    Implements object lookup by index and iteration over the whole list.

    Examples::
    
      >>> ug1 = Fatgraph([Vertex([2,0,0]), Vertex([2,1,1])])
      >>> p = NumberedFatgraphPool(ug1)
      >>> len(p)
      3
      >>> for g in p: print g
      NumberedFatgraph(Fatgraph([Vertex([2, 0, 0]), Vertex([2, 1, 1])]),
                       numbering={BoundaryCycle([(0, 2, 0), (1, 2, 0), (0, 0, 1), (1, 0, 1)]): 0,
                                  BoundaryCycle([(1, 1, 2)]): 1,
                                  BoundaryCycle([(0, 1, 2)]): 2})
      NumberedFatgraph(Fatgraph([Vertex([2, 0, 0]), Vertex([2, 1, 1])]),
                       numbering={BoundaryCycle([(0, 1, 2)]): 0,
                                  BoundaryCycle([(0, 2, 0), (1, 2, 0), (0, 0, 1), (1, 0, 1)]): 1,
                                  BoundaryCycle([(1, 1, 2)]): 2})
      NumberedFatgraph(Fatgraph([Vertex([2, 0, 0]), Vertex([2, 1, 1])]),
                       numbering={BoundaryCycle([(0, 1, 2)]): 0,
                                  BoundaryCycle([(1, 1, 2)]): 1,
                                  BoundaryCycle([(0, 2, 0), (1, 2, 0), (0, 0, 1), (1, 0, 1)]): 2})
       
    Note that, when only one numbering out of many possible ones is
    returned because of isomorphism, the returned numbering may not be
    the trivial one (it is infact the first permutation of 0..n
    returned by `InplacePermutationIterator`)::
      
      >>> ug2 = Fatgraph([Vertex([2,1,0]), Vertex([2,0,1])])
      >>> for g in NumberedFatgraphPool(ug2): print g
      NumberedFatgraph(Fatgraph([Vertex([2, 1, 0]), Vertex([2, 0, 1])]),
                        numbering={BoundaryCycle([(1, 2, 0), (0, 0, 1)]): 0,
                                   BoundaryCycle([(0, 2, 0), (1, 0, 1)]): 1,
                                   BoundaryCycle([(0, 1, 2), (1, 1, 2)]): 2})

    When the graph has only one boundary component, there is only one
    possible numbering, which is actually returned::
    
      >>> ug3 = Fatgraph([Vertex([1,0,1,0])])
      >>> for g in NumberedFatgraphPool(ug3): print g
      NumberedFatgraph(Fatgraph([Vertex([1, 0, 1, 0])]),
                        numbering={BoundaryCycle([(0, 3, 0), (0, 2, 3),
                                                  (0, 1, 2), (0, 0, 1)]): 0})

    Index lookup returns a `NumberedFatgraph` instance, produced
    on-the-fly::

      >>> ug4 = Fatgraph([Vertex([0, 1, 0, 1, 2, 2])])
      >>> pool = NumberedFatgraphPool(ug4)
      >>> pool[0]
      NumberedFatgraph(Fatgraph([Vertex([0, 1, 0, 1, 2, 2])]),
                       numbering={BoundaryCycle([(0, 4, 5)]): 0,
                                  BoundaryCycle([(0, 2, 3), (0, 3, 4),
                                                 (0, 1, 2), (0, 0, 1), (0, 5, 0)]): 1})
    """
    def __init__(self, graph):
        bc = graph.boundary_cycles
        n = len(bc) # == graph.num_boundary_cycles
        orientable = True

        ## Find out which automorphisms permute the boundary cycles among
        ## themselves.
        P = []  #: permutation of boundary cycles induced by `a \in Aut(G)`
        P_ = [] #: corresponding graph automorphisms: `P[i]` is induced by `P_[i]`
        automorphisms = [] #: `NumberedFatgraph` automorphisms
        for a in graph.automorphisms():
            p = Permutation()
            for src in xrange(n):
                dst_cy = bc[src].transform(a, graph)
                try:
                    dst = bc.index(dst_cy)
                except ValueError:
                    # `dst_cy` not in `bc`
                    break # continue with next `a`
                p[src] = dst
            if len(p) != n: # not all `src` were mapped to a `dst`
                continue # with next `a`
            if p.is_identity():
                # `a` preserves the boundary cycles pointwise,
                # so it induces an automorphism of the numbered graph
                automorphisms.append(a)
                if Fatgraph._cmp_orient(graph, graph, a) == -1:
                    orientable = False
            if (p not in P):
                # `a` induces permutation `p` on the set `bc`
                P.append(p)
                P_.append(a)
        assert len(P) > 0 # XXX: should verify that `P` is a group!

        ## There will be as many distinct numberings as there are cosets
        ## of `P` in `Sym(n)`.
        if len(P) > 1:
            def unseen(candidate, P, already):
                """Return `False` iff any of the images of `candidate` by an
                element of group `P` is contained in set `already`.
                """
                for p in P:
                    if p.rearranged(candidate) in already:
                        return False
                return True
            numberings = []
            for candidate in InplacePermutationIterator(range(n)):
                if unseen(candidate, P, numberings):
                    numberings.append(copy(candidate))
        else:
            # if `P` is the one-element group, then all orbits are trivial
            #numberings = PermutationList(n)
            #numberings = list(PermutationIterator(range(n)))
            numberings = [ copy(p) for p in InplacePermutationIterator(range(n)) ]

        # things to remember
        self.graph = graph
        self.is_orientable = orientable
        self.P = P
        self.P_ = P_
        self.numberings = numberings
        self.num_automorphisms = len(automorphisms)


    def __getitem__(self, pos):
        return NumberedFatgraph(self.graph, zip(self.graph.boundary_cycles, self.numberings[pos]))


    def __iter__(self):
        return IndexedIterator(self)
        

    def __len__(self):
        return len(self.numberings)


    def __repr__(self):
        if hasattr(self, 'graph'):
            return "NumberedFatgraphPool(%s)" % self.graph
        else:
            return object.__repr__(self)
        
    __str__ = __repr__
        

    def facets(self, edge, other):
        """Iterate over facets obtained by contracting `edge` and
        projecting onto `other`.

        Each returned item is a triple `(j, k, s)`, where:
          - `j` is the index of a `NumberedFatgraph` in `self`;
          - `k` is the index of a `NumberedFatgraph` in `other`;
          - `s` is the sign by which `self[j].contract(edge)` projects onto `other[k]`.
        Only triples for which `s != 0` are returned.

        Examples::
        
          >>> p0 = NumberedFatgraphPool(Fatgraph([Vertex([1, 2, 0, 1, 0]), Vertex([3, 3, 2])]))
          >>> p1 = NumberedFatgraphPool(Fatgraph([Vertex([0, 1, 0, 1, 2, 2])]))
          >>> list(NumberedFatgraphPool.facets(p0, 2, p1))
          [(0, 0, 1), (1, 1, 1)]
        """
        assert not self.graph.is_loop(edge)
        assert self.is_orientable
        assert other.is_orientable
        
        def push_fwd(f, g1, g2):
            """Return a `Permutation` instance corresponding to the
            map induced by `f` on the boundary cycles, or `None`
            if no such map exists.

            Indeed, there are cases (see examples in the
            `Fatgraph.__eq__` docstring), in which the
            `Fatgraph.isomorphisms` algorithm may find a valid
            mapping, changing from `g1` to an *alternate*
            representation of `g2` - these should fail as they don't
            preserve the boundary cycles.
            """
            result = Permutation() #: map induced by `f` on boundary cycles
            for (src, bc1) in enumerate(g1.boundary_cycles):
                bc2 = bc1.transform(f, g2)
                try:
                    result[src] = g2.boundary_cycles.index(bc2)
                except KeyError:
                    return None
            return result

        def compute_nb_map(pull_back_bcy_map, src, dst):
            """
            Return a pair of `dict` instances:
              - the first maps (indices of ) numberings in
                `NumberedFatgraphPool` instance `src` to corresponding
                (indices of) numberings in `NumberedFatgraphPool`
                instance `dst`.

              - the second maps the same indices into the
                automorphism `a` of `dst.graph` such that::

                  src.numberings[i] = pull_back(<permutation induced by a applied to> dst.numberings[i])

            Raises `KeyError` if no such mapping could be found.

            Indices must be `int`s here, as `list`s cannot be used as
            dictionary keys.
            """
            nb_map = dict()
            a_map = dict()
            for (i, nb) in enumerate(src.numberings):
                (i_, k) = dst._index(pull_back_bcy_map(nb))
                nb_map[i] = i_
                a_map[i] = k
            return (nb_map, a_map)

        g0 = self.graph
        g1 = g0.contract(edge)
        g2 = other.graph
        assert len(g1.boundary_cycles) == len(g2.boundary_cycles)
        
        ## 1. compute map induced on `g0.boundary_cycles` from the
        ##    graph map `f0` which contracts `edge`.
        (e1, e2) = g0.endpoints(edge)
        assert set(g1.boundary_cycles) == set([ bcy.contract(e1, e2, g1) for bcy in g0.boundary_cycles ]), \
               "NumberedFatgraphPool.facets():" \
               " Boundary cycles of contracted graph are not the same" \
               " as contracted boundary cycles of parent graph:" \
               " `%s` vs `%s`" % (g1.boundary_cycles, [ bcy.contract(e1, e2, g1)
                                                        for bcy in g0.boundary_cycles ])
        f0_push_fwd = Permutation(enumerate(
            g1.boundary_cycles.index(bcy.contract(e1, e2, g1))
            for bcy in g0.boundary_cycles
            ))
        assert len(f0_push_fwd) == len(g1.boundary_cycles)
        assert len(f0_push_fwd) == g0.num_boundary_cycles
        if __debug__:
            p1 = NumberedFatgraphPool(g1)
            for (i, nb) in enumerate(self.numberings):
                nb_ = f0_push_fwd.inverse().rearranged(nb)
                i_ = p1._index(nb_)[0]
                assert p1[i_] == self[i].contract(edge)

        ## 2. compute isomorphism map `f1` from `other.graph` to `self.graph.contract(edge)`;
        ##    if there is no such isomorphisms, then stop iteration.
        ##    The reason for having this map go in the opposite direction as `f0`
        ##    is to avoid inverting its push-forward.  (Thus saving a little complexity.)

        # if there is no isomorphism from `g2` to `g1`, the "for" loop
        # down here will not run, so we need a way to detect that a
        # valid isomorphism has been found.
        f1_found = False
        for f1 in Fatgraph.isomorphisms(g2, g1):
            f1_push_fwd = push_fwd(f1, g2, g1)
            if f1_push_fwd is None:
                continue # with next isomorphism
            # now check that every numbering of `g0` can be mapped to
            # a numbering on `g2` through the composite map `f1^(-1) * f0`
            def pull_back(nb): return f1_push_fwd.rearranged(f0_push_fwd.inverse().rearranged(nb))
            try:
                nb_map, a_map = compute_nb_map(pull_back, self, other)
            except KeyError:
                continue # with next isomorphism
            # if we got to this point, the `f1` is a graph
            # automorphism that induces a map on the boundary cycles
            f1_found = True
            break # correct `f1` found
        if not f1_found:
            # no isomorphism between `g1` and `g2`, stop here.
            raise StopIteration

        # check that non-orientable contractions have been caught by
        # the above code.
        if __debug__:
            p1 = NumberedFatgraphPool(g1)
            assert p1.is_orientable
        assert len(f1_push_fwd) == len(g1.boundary_cycles)
        assert len(f1_push_fwd) == len(g2.boundary_cycles)

        ## For every numbering `nb` on `g0`, compute the (index of)
        ## corresponding numbering on `g2` (under the composition map
        ## `f1^(-1) * f0`) and return the triple `(index of nb, index of
        ## push-forward, sign)`.
        for (j, nb) in enumerate(self.numberings):
            k = nb_map[j]
            a = a_map[j]
            ## there are three components to the sign `s`:
            ##   - the sign given by the ismorphism `f1`
            ##   - the sign of the automorphism of `g2` that transforms the
            ##     push-forward numbering into the chosen representative in the same orbit
            ##   - the alternating sign from the homology differential
            s = Fatgraph._cmp_orient(g2, g1, f1) \
                * Fatgraph._cmp_orient(g2, g2, a) \
                * minus_one_exp(g0.edge_numbering[edge])
            if __debug__:
                p2 = NumberedFatgraphPool(g2)
                assert p2[k] == self[j].contract(edge)
                assert s == self[j].contract(edge).projection(p2[k]) * minus_one_exp(g0.edge_numbering[edge])
            yield (j, k, s)

    def _index(self, numbering):
        """
        """
        for (i, p) in enumerate(self.P):
            try:
                j = self.numberings.index(p.rearranged(numbering))
                # once a `p` has matched, there's no reason to try others
                return (j, self.P_[i])
            except ValueError:
                pass
        assert False, \
               "%s._index(%s): No match found." % (self, numbering)
        


## main: run tests

if "__main__" == __name__:
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
