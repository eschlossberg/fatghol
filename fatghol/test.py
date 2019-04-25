#! /usr/bin/env python

import homology_reps
import numpy as np

f = open("M2,2.txt", 'w')
# to_compute = [(0,3),(0,4),(1,1),(1,2),(1,3),(2,1),(0,5)]
to_compute = [(2,2)]
for (g,n) in to_compute:
    mgn = homology_reps.NullSpaceComplex(g,n)
    f.write("g:{} n:{}\n".format(g, n))
    print("g:", g, "n:", n)
    for i in range(len(mgn)):
        f.write("Degree {}:{}\n".format(i, mgn.homology_characters[i]))
        print(i, mgn.homology_characters[i])
    f.write("\n\n")

