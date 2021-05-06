(TeX-add-style-hook
 "phys4004_assignment_2_mpi"
 (lambda ()
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art10"
    "base"
    "pgfplotstable"
    "color"
    "listings")
   (LaTeX-add-labels
    "sec:overview"
    "sec:serial"
    "sec:static"
    "fig:static-load-balance"
    "sec:master-worker"
    "fig:master-load-balance"
    "fig:scaling-chunksize"
    "sec:cyclic"
    "fig:cyclic-load-balance"
    "sec:appendix"
    "sec:serial-code"
    "sec:static-code"
    "sec:master-worker-code"
    "sec:cyclic-code"))
 :latex)

