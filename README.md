# hdpExtra
Extension to the package hdp for mutational signature discovery with the Hierarchical Dirichlet Process

Prior to installing hdpExtra, make sure that dependencies listed in the `DESCRIPTION` file are installed.
To install, use the `devtools` package
```R
devtools::install_github("victorvelasco/hdpExtra")
```

To infer mutational signatures, use the `hdp` package of Nicola Roberts to initialise an HDPM model
with the appropriate hierarchy of Dirichlet Process. To draw samples from the posterior distribution, 
use our function `hdpExtra::hdpExtra_posterior`. This function has been adapted from `hdp::hdp_posterior`,
and saves the clustering allocations to memory so that it is possible to use our post-processing procedure. 
Note that `hdpExtra` keeps full compatibility with `hdp`. Therefore, using `hdp`'s post-processing procedure
with `hdpExtra`'s output is also possible.

Our post-processing procedure allows one to find a representative solution of the vast MCMC output, and 
most importantly to quantify the uncertainty around that solution. 
