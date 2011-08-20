#! /usr/bin/env python

import math
import sys

from decorator import *

def getattr_(obj, name, default_thunk):
    "Similar to .setdefault in dictionaries."
    try:
        return getattr(obj, name)
    except AttributeError:
        default = default_thunk()
        setattr(obj, name, default)
        return default

@decorator
def memoize(func, *args):
    dic = getattr_(func, "memoize_dic", dict)
    # memoize_dic is created at the first call
    if args in dic:
        return dic[args]
    else:
        result = func(*args)
        dic[args] = result
        return result
    
def m_min(g,n):
    return (2*g + n - 1)

def m_max(g,n):
    return (6*g + 3*n -6)

def add_up_to(x, min=0):
    if x == 0 and min == 0:
        yield (0,0)
    elif x-min >= 0:
        for y in xrange(min, x-min+1):
            yield (y, x-y)

def choose2(n):
    return n*(n-1)/2

@memoize
def fact(n):
    if n <= 1:
        return 1
    else:
        return n*fact(n-1)

@memoize
def fact2(n):
    if n <= 1:
        return 1
    else:
        return n*fact(n-2)

def N2(g,n):
    """Number of combinations for graph generation through grafting trees into clovers."""
    m = m_min(g,n)
    #return 2**(l-2) * fact(2*l) / fact(l-1)
    return fact(4*m-2) / fact(2*m) / fact2(2*m-2)

def N3(g,n):
    """Number of combinations for generating graphs from all pairings."""
    m = m_max(g,n)
    #return fact(2*l) / fact(l) / 2**l
    k = 2*m
    r = fact2(k)
    while k > 0:
        r *= (k-1)*(k-2)
        k -= 3
    return r

@memoize
def N1max(g,n):
    """Max. number of combinations for graph generation by bridging lesser genus graphs."""
    if m_max(g,n) <= 0:
        return 0
    if g == 0:
        if n <= 2:
            return 0
        elif n == 3:
            return 2
    if (g,n) == (1,1):
        return 1
    result = 2 * m_max(g,n-1) * N1max(g,n-1)  # hang a circle to an existing edge
    result += 4 * (m_max(g,n-1)**2) * N1max(g,n-1)   # bridging edge separates boundary cycles
    result += 4 * (m_max(g-1,n+1)**2) * N1max(g-1,n+1) # bridging edge does not separate; 1 connected component
    #for gg in add_up_to(g, min=0): # bridging edge does not separate; 2 connected components
    #    for nn in add_up_to(n+1, min=1):
    #        (g1, g2) = gg
    #        (n1, n2) = nn
    #        if (g1, n1) == (g,n) or (g2,n2) == (g,n):
    #            continue
    #        result += 4 * m_max(g1,n1) * m_max(g2,n2) * N3(g1,n1) * N3(g2,n2) 
    return result


@memoize
def N1(g,n):
    """Number of combinations for graph generation by bridging lesser genus graphs (with recursion shortcut)."""
    def estimate(g,n):
        try:
            # use result from actual computation
            return {
                # no. of edges: 9   8   7   6    5    4    3    2
                # -----------------------------------------------
                (0,3): 2,  #                               2 +  1 
                (0,4): 6,  #                6 +  6 +  7 +  6       mmax=6  mmin=3
                (0,5): 26, #   26 +26 +72+103 + 65 + 21            mmax=9  mmin=4
                (1,1): 1,  #                               1 +  1  mmax=3  mmin=2
                (1,2): 5,  #                5 +  5 +  8 +  8       mmax=6  mmin=3
                (1,3): 46, #   46 +46+162+256 +198 + 72            mmax=9  mmin=4
                (2,1): 9,  #    9  +9 +29 +52  +45 + 21            mmax=9  mmin=4
                }[g,n]
        except KeyError:
            # recurse and use as estimate
            return N1(g,n)
    if m_max(g,n) <= 0:
        return 0
    if g == 0:
        if n <= 2:
            return 0
        elif n == 3:
            return 4
    if (g,n) == (1,1):
        return 2
    result = 2 * m_max(g,n-1) * estimate(g,n-1)  # hang a circle to an existing edge
    result += 4 * (m_max(g,n-1)**2) * estimate(g,n-1)   # bridging edge separates boundary cycles
    result += 4 * (m_max(g-1,n+1)**2) * estimate(g-1,n+1) # bridging edge does not separate; 1 connected component
    #for gg in add_up_to(g, min=0): # bridging edge does not separate; 2 connected components
    #    for nn in add_up_to(n+1, min=1):
    #        (g1, g2) = gg
    #        (n1, n2) = nn
    #        if (g1, n1) == (g,n) or (g2,n2) == (g,n):
    #            continue
    #        result += 4 * m_max(g1,n1) * m_max(g2,n2) * N3(g1,n1) * N3(g2,n2) 
    return result


def print_table(func=N3, chi_max=None, g_max=None, n_max=None, maxwidth=12):
    assert chi_max or (g_max and n_max), \
           "print_table: Need either `chi_max` or `g_max` and `n_max`."
    if maxwidth is None or maxwidth == 0:
        format_small = "%10d"
        format_large = "%d"
    else:
        format_small = "%%%dd" % maxwidth
        format_large = "%E"
    if g_max is None:
        g_max = chi_max - 1
    if n_max is None:
        n_max = chi_max
    print "      " + str.join(" ", [ (format_small % n) for n in xrange(1, n_max+1) ])
    print "-------" + ("-" * (maxwidth+1) * n_max)
    def fmt_entry(num):
        digits = 1 + math.log10(1 + num)
        if digits > maxwidth:
            return (format_large % num)
        else:
            return (format_small % num)
    for g in xrange(0,g_max+1):
        fmt = "g=" + ("%-4d " % g) + (" %s" * n_max)
        if chi_max is not None:
            this_n_max = max(chi_max - 2*g, 0)
        else:
            this_n_max = n_max
        r = [ fmt_entry(func(g,n)) for n in xrange(1, this_n_max+1) ]
        if this_n_max < n_max:
            r += [ "" ] * (n_max - this_n_max)
        print fmt % tuple(r)


def print_list(func=N3, chi_max=None, g_max=None, n_max=None, maxwidth=12):
    assert chi_max or (g_max and n_max), \
           "print_list: Need either `chi_max` or `g_max` and `n_max`."
    if maxwidth is None or maxwidth == 0:
        format_small = "%10d"
        format_large = "%d"
    else:
        format_small = "%%%dd" % maxwidth
        format_large = "%E"
    if g_max is None:
        g_max = chi_max - 1
    if n_max is None:
        n_max = chi_max
    print "-------" + ("-" * (maxwidth+1))
    def fmt_entry(num):
        digits = 1 + math.log10(1 + num)
        if digits > maxwidth:
            return (format_large % num)
        else:
            return (format_small % num)
    for g in xrange(0,g_max+1):
        fmt = "g=" + ("%-4d " % g) + (" %s" * n_max)
        if chi_max is not None:
            this_n_max = max(chi_max - 2*g, 0)
        else:
            this_n_max = n_max
        for n in xrange(1, this_n_max + 1):
            print ("g=%d,n=%d %s" % (g, n, fmt_entry(func(g,n))))


## main

if len(sys.argv) <= 1:
    max_chi = 8
    max_g = None
    max_n = None
elif len(sys.argv) == 2:
    max_chi = int(sys.argv[1])
    max_g = None
    max_n = None
else:
    max_chi = None
    max_g = int(sys.argv[1])
    max_n = int(sys.argv[2])

for fn in N1, N2, N3:
    print fn.__doc__
    #print_table(fn, max_chi, max_g, max_n)
    print_list(fn, max_chi, max_g, max_n)
    print

