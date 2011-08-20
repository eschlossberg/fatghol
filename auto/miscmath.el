(TeX-add-style-hook "miscmath"
 (lambda ()
    (TeX-add-symbols
     '("Hermitian" ["argument"] 0)
     '("EndOp" ["argument"] 0)
     '("catVect" ["argument"] 0)
     '("avg" ["argument"] 1)
     '("pairing" ["argument"] 2)
     '("ipang" ["argument"] 2)
     '("ipsq" ["argument"] 2)
     '("ip" ["argument"] 2)
     '("inner" ["argument"] 2)
     '("correlator" 1)
     '("card" 1)
     '("norm" 1)
     '("abs" 1)
     '("third" 1)
     '("half" 1)
     '("oneof" 1)
     '("Perm" 1)
     '("sbtext" 1)
     '("sptext" 1)
     '("xp" 1)
     '("tp" 1)
     '("ldl" 1)
     '("rdl" 1)
     '("cplx" 1)
     '("lspan" 1)
     '("lie" 1)
     '("sheaf" 1)
     '("field" 1)
     '("operad" 1)
     '("category" 1)
     "setC"
     "setN"
     "setQ"
     "setZ"
     "setR"
     "fk"
     "k"
     "del"
     "delbar"
     "ud"
     "Htop"
     "Stab"
     "Diff"
     "Def"
     "isoarrow"
     "opp"
     "onehalf"
     "onethird"
     "twothirds"
     "I"
     "U"
     "E"
     "inv"
     "geq"
     "leq"
     "oo"
     "catSet"
     "catMod"
     "catBraid"
     "M"
     "Mbar"
     "Mcomb"
     "Mbarcomb")
    (TeX-run-style-hooks
     "eucal"
     "mathrsfs"
     "amsfonts"
     "amsxtra"
     "amssymb"
     "amsmath")))

