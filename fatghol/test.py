#! /usr/bin/env python

import homology_reps
import numpy as np
import os

filename = "M2,2.txt"
f = open(filename, 'w')
# to_compute = [(0,3),(0,4),(1,1),(1,2),(1,3),(2,1),(0,5)]
to_compute = [(2,2)]
for (g,n) in to_compute:
    mgn = homology_reps.NullSpaceComplex(g,n)
    f.write("g:{} n:{}\n".format(g, n))
    f.close()
    #print("g:", g, "n:", n)
    for i in range(len(mgn)):
        mgn.compute_null_space(i)
        np.save(".null_spaces/d{}.mat".format(i), mgn.null_spaces[i][1])
    for i in range(len(mgn)):
        mgn.null_space_character(i)
        f = open('.M22nsc.txt', 'a')
        f.write(mgn.null_space_characters[i])
        f.close()
    mgn.compute_homology_characters()
    for i in range(len(mgn)):
        f = open(filename, 'a')
        f.write("Degree {}:{}\n".format(i, mgn.homology_characters[i]))
        #print(i, mgn.homology_characters[i])
        f.close()
    f = open(filename, 'a')
    f.write("\n\n")
f.close()
