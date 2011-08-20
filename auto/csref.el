(TeX-add-style-hook "csref"
 (lambda ()
    (TeX-add-symbols
     '("csref" 1)
     '("cslabel" 1)
     '("renewcsref" 1)
     '("newcsref" 1)
     '("defcsref" 1)
     "filename"
     "label"
     "sectionrefname")
    (TeX-run-style-hooks
     "prettyref"
     "fancyref"
     "varioref")))

