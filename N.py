#! /usr/bin/env python

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
    
def l_min(g,n):
    return (2*g + n - 1)

def l_max(g,n):
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

def N1(g,n):
    """Number of combinations for graph generation through grafting trees into clovers."""
    l = l_min(g,n)
    return 2**(l-2) * fact(2*l) / fact(l-1)

def N2(g,n):
    """Number of combinations for generating graphs from all pairings."""
    l = l_max(g,n)
    return fact(2*l) / fact(l) / 2**l

@memoize
def N3(g,n):
    """Number of combinations for graph generation by bridging lesser genus graphs."""
    if l_max(g,n) <= 0:
        return 0
    if g == 0:
        if n <= 2:
            return 0
        elif n == 3:
            return 4
    if (g,n) == (1,1):
        return 1
    result = 2 * l_max(g,n-1) * N3(g,n-1)  # hang a circle to an existing edge
    result += 4 * choose2(l_max(g-1,n-1)) * N3(g,n-1)   # bridging edge separates boundary cycles
    result += 4 * choose2(l_max(g-1,n+1)) * N3(g-1,n+1) # bridging edge does not separate; 1 connected component
    for gg in add_up_to(g, min=0): # bridging edge does not separate; 2 connected components
        for nn in add_up_to(n+1, min=1):
            (g1, g2) = gg
            (n1, n2) = nn
            if (g1, n1) == (g,n) or (g2,n2) == (g,n):
                continue
            result += 4 * l_max(g1,n1) * l_max(g2,n2) * N3(g1,n1) * N3(g2,n2) 
    return result

def print_table(g_max, n_max, func=N3):
    print "      " + str.join(" ", [ ("%18d" % n) for n in xrange(1, n_max+1) ])
    print "-------" + ("-" * 19 * n_max)
    for g in xrange(0,g_max+1):
        fmt = "g=" + ("%-4d " % g) + (" %18d" * n_max)
        r = tuple([ func(g,n) for n in xrange(1, n_max+1) ])
        print fmt % r


## main

try:
    max_g = int(sys.argv[1])
except:
    max_g = 5

try:
    max_n = int(sys.argv[2])
except:
    max_n = 4

for fn in N1, N2, N3:
    print fn.__doc__
    print_table(max_g, max_n, fn)
    print

