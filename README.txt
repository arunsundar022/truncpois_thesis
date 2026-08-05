truncpois - Thesis Folder
=====================================

This folder contains the thesis document (Truncpois_Thesis.Rmd / .pdf) for
the truncpois R package, along with the code and package files needed to
reproduce its results.

1. Installing the package
--------------------------
The truncpois package can be installed in either of two ways:

  a) From the provided source tarball in this folder:
       install.packages("truncpois_0.1.0.tar.gz", repos = NULL, type = "source")

  b) From GitHub, via pak:
       # install.packages("pak")
       pak::pak("arunsundar022/truncpois")

2. Reproducing the thesis results
-----------------------------------
"Thesis code.R" reproduces every experiment, figure, and table reported in
the thesis. It requires the truncpois package (see above) plus the
LaplacesDemon and microbenchmark packages to be installed. Running it
end-to-end regenerates all of the PNG figures into the Figures/ folder.

3. Further examples
---------------------
Beyond what is reproduced in "Thesis code.R", the truncpois package itself
ships with additional worked examples in its function documentation
(e.g. help(dtruncpois)) and in its vignette:
  vignette("truncpois", package = "truncpois")

The complete package source, tests, and documentation are also available on
GitHub at https://github.com/arunsundar022/truncpois (see Appendix A of the
thesis).
