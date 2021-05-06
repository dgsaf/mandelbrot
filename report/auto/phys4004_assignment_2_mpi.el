(TeX-add-style-hook
 "phys4004_assignment_2_mpi"
 (lambda ()
   (add-to-list 'LaTeX-verbatim-environments-local "lstlisting")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "lstinline")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "lstinline")
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
    "fig:master-load-balance-100000"
    "fig:scaling-chunksize"
    "fig:master-load-balance-10"
    "fig:master-load-balance-10000000"
    "sec:cyclic"
    "fig:cyclic-load-balance"
    "sec:appendix"
    "sec:serial-code"
    "sec:static-code"
    "sec:master-worker-code"
    "sec:cyclic-code")
   (LaTeX-add-listings-lstdefinestyles
    "ff"))
 :latex)

