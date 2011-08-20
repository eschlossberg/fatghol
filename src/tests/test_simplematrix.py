#! /usr/bin/env python

from simplematrix import SimpleMatrix

m = SimpleMatrix(100, 100)

i = 0
j = 0
for x in range(16):
    m.entry(i, j, x)
    i += 13; i %= 100
    j += 29; j %= 100

r = m.rank()
assert r == 15
print (r)
