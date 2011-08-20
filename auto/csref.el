(TeX-add-style-hook "csref"
 (lambda ()
    (TeX-add-symbols
     '("csref" 1)
     '("cslabel" 1)
     '("renewcsref" 1)
     '("newcsref" 1)
     "filename"
     "label"
     "refsectionname")
    (TeX-run-style-hooks
     "prettyref"
     "fancyref"
     "varioref")))

