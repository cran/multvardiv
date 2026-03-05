# multvardiv 1.0.15

* The title of the package is modified, to highlight the three
  probability distributions covered by the 'multvardiv' package.
  Now "Multivariate Generalized Gaussian Distribution,
  Multivariate t Distribution, Multivariate Cauchy Distribution,
  Statistical Divergence"

# multvardiv 1.0.13

* The diststudent function now returns `NaN`, with a warning, if `delta1 + delta2 - p/2 <= 0`
  When $\delta_1 + \delta_2 - \frac{p}{2} \leq 0$, the Rényi divergence
  (and, therefore, the Bhattacharyya or Hellinger divergence)
  cannot be computed using the expression implemented in diststudent.
* Added a `NEWS.md` file to track changes to the package.
