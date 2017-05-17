# Stata implementation of clubSandwich

`clubSandwich` provides cluster-robust variance estimators (i.e., sandwich
estimators) for ordinary and weighted least squares linear regression models. 
Several adjustments are incorporated to improve small-sample performance. The
package includes functions for estimating linear regression models with
cluster-robust variance-covariance matrices and for testing single- and
multiple-contrast hypotheses based on Wald test statistics. Variance-covariance
estimators are based on a version of the bias-reduced linearization estimator
proposed by Bell and McCaffrey (2002) and further developed by Tipton and
Pustejovsky (2015) and Pustejovsky and Tipton (2016). Tests of single regression
coefficients use Satterthwaite corrections. Tests of multiple-contrast
hypotheses use an approximation to Hotelling's T-squared distribution.

### Requirements

* Stata version 14.2 
* matsqrt from http://www.stata.com/users/jpitblado

### Authors

Marcelo Tyszler. Sustainable Economic Development and Gender, Royal Tropical Institute, Netherlands. m.tyszler@kit.nl

James E. Pustejovsky (maintainer). Department of Education Psychology, University of Texas at Austin. pusto@austin.utexas.edu

Elizabeth Tipton. Department of Human Development, Teachers College, Columbia University. tipton@tc.columbia.edu

# Installing clubSandwich

The package is available on the SSC Archive. To install it, type 
```
ssc install clubSandwich 
```

To install the latest development version directly from Github, type: 
```
net install github, from("https://haghish.github.io/github/") 
github install jepusto/clubSandwich-Stata 
```