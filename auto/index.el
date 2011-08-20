(TeX-add-style-hook "index"
 (lambda ()
    (LaTeX-add-bibliographies
     "math")
    (LaTeX-add-environments
     "theorem"
     "lemma"
     "proposition"
     "claim"
     "conjecture"
     "corollary"
     "remark"
     "definition"
     "example")
    (TeX-add-symbols
     '("ginttriv" ["argument"] 1)
     '("gint" ["argument"] 1)
     '("biggint" ["argument"] 1)
     '("Tcomb" ["argument"] 0)
     '("Strebel" ["argument"] 0)
     '("SG" ["argument"] 0)
     '("ERG" ["argument"] 0)
     '("RG" ["argument"] 0)
     '("RTG" ["argument"] 0)
     '("RTE" ["argument"] 0)
     '("RTD" ["argument"] 0)
     '("Oo" ["argument"] 0)
     '("vv" 1)
     '("Holes" 1)
     '("Legs" 1)
     '("Edges" 1)
     '("Vertices" 1)
     '("btpc" 1)
     '("tpc" 1)
     '("freemsc" 1)
     '("freemc" 1)
     '("textsuperscript" 1)
     '("textsubscript" 1)
     "theoremname"
     "ksz"
     "aksz"
     "A"
     "B"
     "rev"
     "catN"
     "Ao"
     "Conf"
     "Lb"
     "Hn"
     "Harmonic"
     "ad"
     "comb"
     "T"
     "rel"
     "xyc")
    (TeX-run-style-hooks
     "miscmath"
     "amsthm"
     "xytree"
     "rg"
     "xy"
     "ps"
     "dvips"
     "xdvi"
     "colour"
     "all"
     "arc"
     "knot"
     "poly"
     "paralist"
     "newitem"
     "newenum"
     "url"
     "hyperref"
     "graphicx"
     "final"
     "csref"
     "debug"
     "latex2e"
     "amsbook10"
     "amsbook"
     "a4paper"
     "reqno"
     "openany"
     "draft"
     "title"
     "prelim"
     "graphcomplex"
     "gc"
     "fd"
     "ainfty"
     "kontsevich"
     "arrows"
     "btc"
     "rt"
     "operads")))

