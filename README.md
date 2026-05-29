# GWnorm

Simulation code and examples to accompany 

_A new way to evaluate G-Wishart normalising constants via Fourier analysis_

Ching Wong, Giusi Moffa and Jack Kuipers (2024), [doi:10.48550/arXiv.2404.06803](https://arxiv.org/abs/2404.06803)

### GWnorm package

The code relies on the [_GWnorm_](https://CRAN.R-project.org/package=GWnorm) R package.

### Examples

Examples for the first figure can be run with the `examples_figure1.R` file, while the **iris_example** directory hosts the results for the Fisher Iris data example.

### Simulations

The main simulations run from the files `comparison_prime.R` and `comparison_density.R`, with the simulations results stored in the **sim_results** folder and plotted in the **plots** folder. 

### Special cases

Code to check and count the special cases is in the **special_case_count** folder, with the `.Rdata` file storing the results for Table 1. This is superceded by functionality in the [_GWnorm_](https://CRAN.R-project.org/package=GWnorm) package to evaluate the normalising constant for those cases, as in the Examples above.
