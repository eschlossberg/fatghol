(TeX-add-style-hook "figures"
 (function
  (lambda ()
    (LaTeX-add-environments
     "eps")
    (TeX-add-symbols
     '("ginttriv" 2)
     '("gint" 2)
     '("biggint" 2)
     '("vv" 1)
     '("Holes" 1)
     '("Edges" 1)
     '("Vertices" 1)
     '("Ribbon" 1)
     '("freemsc" 1)
     '("freemc" 1)
     '("SG" 1)
     '("RG" 1)
     '("RT" 1)
     '("RTD" 1)
     '("Oo" 1)
     "A"
     "B"
     "rev"
     "catN")
    (TeX-run-style-hooks
     "miscmath"
     "rg"
     "xy"
     "color"
     "latex2e"
     "amsart10"
     "amsart"
     "a4paper"
     "fleqn"
     "reqno"))))

