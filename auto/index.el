(TeX-add-style-hook "index"
 (lambda ()
    (LaTeX-add-bibliographies
     "math"
     "tech")
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
     '("Hermitian" ["argument"] 0)
     '("Strebel" ["argument"] 0)
     '("SG" ["argument"] 0)
     '("RTG" ["argument"] 0)
     '("RTE" ["argument"] 0)
     '("RTD" ["argument"] 0)
     '("Oo" ["argument"] 0)
     '("vv" 1)
     '("expval" 1)
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
     '("FIXME" 1)
     "theoremname"
     "ksz"
     "aksz"
     "A"
     "B"
     "rev"
     "catN"
     "Ao"
     "RG"
     "ERG"
     "Conf"
     "Lb"
     "Hn"
     "Harmonic"
     "ad"
     "comb"
     "T"
     "Tcomb"
     "M"
     "Mbar"
     "Mcomb"
     "Mbarcomb"
     "R"
     "negquad"
     "negqquad"
     "X"
     "dec"
     "res"
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
     "csref"
     "debug"
     "url"
     "paralist"
     "newitem"
     "newenum"
     "listings"
     "hyperref"
     "graphicx"
     "final"
     "mdwlist"
     "fontenc"
     "T1"
     "color"
     "prelim2e"
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
     "algorithm"
     "gc"
     "fd"
     "ainfty"
     "kontsevich"
     "dfiz"
     "arrows"
     "btc"
     "rt"
     "operads"
     "python")))

