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
     '("ginttriv" 2)
     '("gint" 2)
     '("biggint" 2)
     '("vv" 1)
     '("Holes" 1)
     '("Edges" 1)
     '("Vertices" 1)
     '("Strebel" 1)
     '("ftpc" 1)
     '("btpc" 1)
     '("tpc" 1)
     '("freemsc" 1)
     '("freemc" 1)
     '("SG" 1)
     '("RG" 1)
     '("RTG" 1)
     '("RTE" 1)
     '("RTD" 1)
     '("Oo" 1)
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
     "csref"
     "nomencl"
     "hyperref"
     "url"
     "mdwtab"
     "mdwlist"
     "paralist"
     "rg"
     "xy"
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

