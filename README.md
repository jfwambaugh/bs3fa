# About bs3fa

`bs3fa` is a package for Bayesian sparse and smooth supervised factor analysis. This model is appropriate for data in which you observe functional Y and numeric (continuous, binary, count) data X. The model assumes all variation in Y is explained by some low-dimensional factors eta, and these factors also explain part (but not all) of the variation in X.

# Navigating this repository

The [R](R) directory contains R functions available in this package, including `run_sampler()`, the main model sampler, and the post-processing functions used to resolve label and sign switches.

The [src](src) directory contains cpp source code used for sampling specific parameters in the Gibbs sampler (i.e., this directory contains the sampling functions used within the sampler loop of `run_sampler()`).

The [demo](demo) directory contains a demo of the method using realistically simulated active chemicals, and continuous chemical features with sparsit in the toxicity-relevant loadings.

The [data](data) directory contains a cleaned set of samples from the ToxCast Attagene PXR assay measured at a common grid of dose values. It is used in the `simulate_data()` function if the user desires more realistic loading vectors be simulated.
