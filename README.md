# Stata implementation of clubSandwich

`clubSandwich` provides bias-reduced linearization variance estimators (i.e.,
sandwich estimators or cluster-robust variance estimators) for ordinary and
weighted least squares linear regression models. Several adjustments are
incorporated to improve small-sample performance. The package includes functions
for estimating the variance-covariance matrix and for testing single- and
multiple-contrast hypotheses based on Wald test statistics. Tests of single
regression coefficients use Satterthwaite corrections. Tests of
multiple-contrast hypotheses use an approximation to Hotelling's T-squared
distribution. 

# Installing clubSandwich

The package is available on the SSC Archive. To install it, type 
```{r}
install.packages("clubSandwich")
```

To install the latest development version directly from Github, type:
```{r}
net install github, from("https://haghish.github.io/github/")
github install jepusto/clubSandwich-Stata
```