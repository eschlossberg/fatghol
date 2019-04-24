#! /usr/bin/env python

import homology_reps
import numpy as np

f = open("/home/schlossberg2/homology_characters.txt", 'w')
for (g,n) in [(0,3),(0,4),(1,1),(1,2),(1,3),(2,1),(0,5)]:
    mgn = homology_reps.NullSpaceComplex(g,n)
    f.write("g:{} n:{}\n".format(g, n))
    print("g:", g, "n:", n)
    for i in range(len(mgn)):
        f.write("Degree {}:{}\n".format(i, mgn.homology_characters[i]))
        print(i, mgn.homology_characters[i])
    f.write("\n\n")
