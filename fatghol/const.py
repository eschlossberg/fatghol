#! /usr/bin/env python
#
"""
Compute well-known invariants of `M_{g,n}`.
"""
#
#   Copyright (C) 2008-2012 Riccardo Murri <riccardo.murri@gmail.com>
#   All rights reserved.
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
__docformat__ = 'reStructuredText'


from fractions import Fraction
from fatghol.combinatorics import (
    bernoulli,
    factorial,
    minus_one_exp,
    )


def orbifold_euler_characteristics(g,n):
    """
    Return the orbifold/virtual Euler characteristics of `M_{g,n}`,
    computed according to Harer-Zagier.
    """
    if g==0:
        return factorial(n-3) * minus_one_exp(n-3)
    elif g==1:
        return Fraction(minus_one_exp(n), 12) * factorial(n-1)
    else: # g > 1
        return bernoulli(2*g) * factorial(2*g+n-3) / (factorial(2*g-2) * 2*g) * minus_one_exp(n)


def euler_characteristics(g,n):
    """
    Return Euler characteristics of `M_{g,n}`.

    The Euler characteristics is computed according to formulas and
    tables found in:
    * Bini-Gaiffi-Polito, arXiv:math/9806048, p.3
    * Bini-Harer, arXiv:math/0506083, p. 10
    """
    if g==0:
        # according to Bini-Gaiffi-Polito arXiv:math/9806048, p.3
        return factorial(n-3)*minus_one_exp(n-3)
    elif g==1:
        # according to Bini-Gaiffi-Polito, p. 15
        if n>4:
            return factorial(n-1)*Fraction(minus_one_exp(n-1),12)
        else:
            es = [1,1,0,0]
            return es[n-1] # no n==0 computed in [BGP]
    elif g==2:
        # according to Bini-Gaiffi-Polito, p. 14
        if n>6:
            return factorial(n+1)*Fraction(minus_one_exp(n+1),240)
        else:
            es = [1,2,2,0,-4,0,-24]
            return es[n]
    elif g>2:
        # according to Bini-Harer arXiv:math/0506083, p. 10
        es = [# n==1  n==2     n==3     n==4        n==5        n==6          n==7          n==8
              [    8,    6,       4,     -10,         30,       -660,         6540,        79200], # g==3
              [   -2,  -10,     -24,     -24,       -360,       2352,       -37296,       501984], # g==4
              [   12,   26,      92,     182,       1674,     -16716,       238980,     -3961440], # g==5
              [    0,  -46,    -206,     188,      -7512,     124296,     -2068392,     37108656], # g==6
              [   38,  120,     676,   -1862,      71866,   -1058676,     21391644,   -422727360], # g==7
              [ -166, -630,   -5362,   16108,    -680616,   12234600,   -259464240,   5719946400], # g==8
              [  748, 2132,   29632, -323546,    7462326, -164522628,   3771668220, -90553767840], # g==9
              [-1994, 6078, -213066, 4673496, -106844744, 2559934440, -64133209320,    1.664e+12], # g==10
              ]
        # g=0,1,2 already done above
        return es[g-3][n]
    else:
        raise ValueError("No Euler characteristics known for M_{g,n},"
                         " where g=%s and n=%s" % (g,n))



## main: run tests

if "__main__" == __name__:
    import doctest
    doctest.testmod(name="const",
                    optionflags=doctest.NORMALIZE_WHITESPACE)
