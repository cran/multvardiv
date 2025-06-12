# multvardiv

The **multvardiv** packages replaces the [mggd](https://forge.inrae.fr/imhorphen/mggd), [mcauchyd](https://forge.inrae.fr/imhorphen/mcauchyd) and [mstudentd](https://forge.inrae.fr/imhorphen/mstudentd) packages.

This package provides tools for multivariate probability distributions:

* Multivariate Gaussian distribution (MGGD)
* Multivariate Cauchy distribution (MCD)
* Multivariate $t$ distribution (MTD)

For each of these probability distributions, some functions are available to compute the divergence between two distributions:

* Renyi, Bhattacharyya and Hellinger divergence
* Kullback-Leibler divergence

And some functions for the manipulation of probability distributions:
  + Probability density
  + Parameter estimation
  + Sample simulation
  + Density plot (bivariate distribution)

# Installation

Install the package from CRAN:
```
install.packages("multvardiv")
```

Or install the development version from the repository, using the [`devtools`](https://CRAN.R-project.org/package=devtools) package:

```
install.packages("devtools")
devtools::install_git("https://forge.inrae.fr/imhorphen/multvardiv")
```

## Authors

[Pierre Santagostini](mailto:pierre.santagostini@institut-agro.fr) and [Nizar Bouhlel](mailto:nizar.bouhlel@institut-agro.fr)

## Reference

N. Bouhlel, A. Dziri,
Kullback-Leibler Divergence Between Multivariate Generalized Gaussian Distributions.
IEEE Signal Processing Letters, vol. 26 no. 7, July 2019.
<https://doi.org/10.1109/LSP.2019.2915000>

N. Bouhlel, D. Rousseau,
A Generic Formula and Some Special Cases for the Kullback–Leibler Divergence between Central Multivariate Cauchy Distributions.
Entropy, 24, 838, July 2022.
<https://doi.org/10.3390/e24060838>

N. Bouhlel and D. Rousseau (2023), Exact Rényi and Kullback-Leibler Divergences Between Multivariate t-Distributions, IEEE Signal Processing Letters.
<https://doi.org/10.1109/LSP.2023.3324594>
