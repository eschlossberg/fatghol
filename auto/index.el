(TeX-add-style-hook "index"
 (function
  (lambda ()
    (LaTeX-add-bibliographies
     "math")
    (LaTeX-add-environments
     "theorem"
     "thmsec"
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
     '("Strebel" ["argument"] 0)
     '("SG" ["argument"] 0)
     '("RG" ["argument"] 0)
     '("RTG" ["argument"] 0)
     '("RTE" ["argument"] 0)
     '("RTD" ["argument"] 0)
     '("Oo" ["argument"] 0)
     '("vv" 1)
     '("Holes" 1)
     '("Edges" 1)
     '("Vertices" 1)
     '("ftpc" 1)
     '("btpc" 1)
     '("tpc" 1)
     '("freemsc" 1)
     '("freemc" 1)
     "theoremname"
     "A"
     "B"
     "rev"
     "catN"
     "Conf"
     "T"
     "Lb"
     "Hn"
     "xyc"
     "nx")
    (TeX-run-style-hooks
     "graphicx"
     "final"
     "csref"
     "nomencl"
     "refpage"
     "hyperref"
     "url"
     "mdwtab"
     "mdwlist"
     "paralist"
     "newitem"
     "newenum"
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
     "miscmath"
     "amsthm"
     "rcs"
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
     "gc"
     "fd"
     "kontsevich"
     "arrows"
     "btc"
     "rt"))))

