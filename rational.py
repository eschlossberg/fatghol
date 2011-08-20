#! /usr/bin/env python
#
# This is a straightforward port of the `boost::rational` C++ class, whose
# source file bears the following copyright notice:
#
#   (C) Copyright Paul Moore 1999. Permission to copy, use, modify, sell and
#   distribute this software is granted provided this copyright notice appears
#   in all copies. This software is provided "as is" without express or
#   implied warranty, and with no claim as to its suitability for any purpose.
#   
#   See http://www.boost.org/libs/rational for documentation.
#
"""Rational exact arithmetic in Python.

Implementation notes: Fractions are kept in normalized form at all
times. Normalized form is defined as `gcd(numerator, denominator) ==
1` and `self.denominator > 0`.  In particular, note that the
implementation of `abs()` below relies on `self.denominator` always
being positive.

This is a straightforward port of Boost's `boost::rational` C++ class,
which see for additional comments on the algorithms used.
"""
__docformat__ = 'reStructuredText'



def gcd(m, n):
    """Return the GCD of integers `m` and `n`."""
    # ensure `m` and `n` are positive
    if m < 0:
        m = -m
    if n < 0:
        n = -n
    while True:
        if 0 == m:
            return n
        n %= m
        if 0 == n:
            return m
        m %= n


def lcm(m, n):
    """Return the LCM of integers `m` and `n`."""
    if (0 == m) or (0 == n):
        return 0
    n /= _gcd(m,n)
    n *= m
    if n < 0:
        n = -n
    return n


class Rational(object):
    """A rational number.
    Rational(p,q) -> p/q
    Rational(p) -> p/1
    """

    __slots__ = [ 'numerator', 'denominator' ]
    
    def __init__(self, numerator, denominator=1,
                 _normalize=True):
        """Create a new rational number from a pair of integers."""
        self.numerator = numerator
        self.denominator = denominator

        if _normalize:
            self.normalize()


    def normalize(self):
        if (self.denominator == 0):
            raise ZeroDivisionError

        # handle the case of zero separately, to avoid division by zero
        if (self.numerator == 0):
            self.denominator = 1
            return self

        div = gcd(self.numerator, self.denominator)
        self.numerator /= div
        self.denominator /= div

        # ensure that the denominator is positive
        if (self.denominator < 0):
            self.numerator = -self.numerator;
            self.denominator = -self.denominator

        return self


    def __abs__(self):
        if self.numerator >= 0:
            return self
        
        return Rational(-self.numerator, self.denominator, _normalize=False)
        

    def __add__(self, other):
        new = Rational(self.numerator, self.denominator, _normalize=False)
        new += other
        return new


    def __cmp__(self, other):
        try:
            # If the two values have different signs, we don't need to do the
            # expensive calculations below. We take advantage here of the fact
            # that the denominator is always positive.
            if (self.numerator < 0 and other.numerator >= 0): # -ve < +ve
                return -1
            if (self.numerator > 0 and other.numerator < 0): # +ve or zero is not < -ve or zero
                return +1

            # Avoid overflow
            gcd1 = gcd(self.numerator, other.numerator)
            if gcd1 == 0:
                gcd1 = 1
            gcd2 = gcd(other.denominator, self.denominator)
            return cmp((self.numerator/gcd1) * (other.denominator/gcd2),
                       (self.denominator/gcd2) * (other.numerator/gcd1))

        except AttributeError: # assume `other` is a whole number
            # If the two values have different signs, we don't need to do the
            # expensive calculations below. We take advantage here of the fact
            # that the denominator is always positive.
            if (self.numerator < 0 and other >= 0): # -ve < +ve
                return -1
            if (self.numerator >= 0 and other < 0): # +ve or zero is not < -ve or zero
                return +1

            # Now, use the fact that n/d truncates towards zero as long as n and d
            # are both positive.
            # Divide instead of multiplying to avoid overflow issues. Of course,
            # division may be slower, but accuracy is more important than speed...
            if (self.numerator > 0):
                return cmp((self.numerator/self.denominator), other)
            else:
                return cmp(-other, (-self.numerator/self.denominator))


    def __coerce__(self, other):
        return (self, Rational(other))
    

    def __div__(self, other):
        new = Rational(self.numerator, self.denominator, _normalize=False)
        new /= other
        return new
    

    def __eq__(self, other):
        try:
            return (self.numerator == other.numerator) \
                   and (self.denominator == other.denominator)

        except AttributeError: # presume `other` is a whole number
            return (self.numerator == other) and (self.denominator == 1)
            
        
    def __hash__(self):
        return hash((self.numerator, self.denominator))


    def __iadd__(self, other):
        try:
            # see boost::rational for an explanation of this algorithm

            # protect against self-modification
            other_num = other.numerator
            other_den = other.denominator

            g = gcd(self.denominator, other_den)
            self.denominator /= g
            self.numerator = self.numerator * (other_den / g) \
                             + other_num * self.denominator
            g = gcd(self.numerator, g)
            self.numerator /= g
            self.denominator *= other_den/g
            
        except AttributeError: # assume `other` is a whole number
            # assume `other` is some kind of integer
            self.numerator += self.denominator * other

        return self


    def __idiv__(self, other):
        try:
            # protect against self-modification
            other_num = other.numerator
            other_den = other.denominator
            
            # trap division by zero
            if (other_num == 0):
                raise ZeroDivisionError
            
            # shortcut
            if (self.numerator == 0):
                return self
        
            # Avoid overflow and preserve normalization
            gcd1 = gcd(self.numerator, other_num)
            gcd2 = gcd(other_den, self.denominator)
            self.numerator = (self.numerator/gcd1) * (other_den/gcd2)
            self.denominator = (self.denominator/gcd2) * (other_num/gcd1)

            if (self.denominator < 0):
               self.numerator = -self.numerator
               self.denominator = -self.denominator

        except AttributeError:
            self.numerator /= other

            if (self.denominator < 0):
               self.numerator = -self.numerator
               self.denominator = -self.denominator
            
        return self

           
    def __imul__(self, other):
        try:
            # avoid overflow and preserve normalization
            gcd1 = gcd(self.numerator, other.denominator)
            gcd2 = gcd(other.numerator, self.denominator)
            self.numerator = (self.numerator/gcd1) * (other.numerator/gcd2)
            self.denominator = (self.denominator/gcd2) * (other.denominator/gcd1)

        except AttributeError: # assume `other` is a whole number
            self.numerator *= other

        return self


    def __int__(self):
        return self.numerator // self.denominator

    __long__ = __int__


    def __isub__(self, other):
        try:
            # see boost::rational for an explanation of this algorithm

            # protect against self-modification
            other_num = other.numerator
            other_den = other.denominator

            g = gcd(self.denominator, other_den)
            self.denominator /= g
            self.numerator = self.numerator * (other_den / g) - other_num * self.denominator
            g = gcd(self.numerator, g)
            self.numerator /= g
            self.denominator *= other_den/g
            
        except AttributeError: # assume `other` is a whole number
            self.numerator -= self.denominator * other

        return self


    def __mul__(self, other):
        new = Rational(self.numerator, self.denominator, _normalize=False)
        new *= other
        return new
    

    def __neg__(self):
        return Rational(-self.numerator, self.denominator, _normalize=False)


    def _nonzero__(self):
        if self.numerator != 0:
            return True
        else:
            return False

    def __pos__(self):
        return self


    def __repr__(self):
        return "Rational(%s, %s)" %  (self.numerator, self.denominator)


    def __str__(self):
        return "%s/%s" % (self.numerator, self.denominator)


    def __sub__(self, other):
        new = Rational(self.numerator, self.denominator, _normalize=False)
        new -= other
        return new


## main: run tests

if "__main__" == __name__:
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
