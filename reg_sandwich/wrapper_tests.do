/* ****************************************************************
reg_sandwich wrapper caller test - 17/Jan/2017
github repository: https://github.com/jepusto/clubSandwich-Stata
Updated by Marcelo Tyszler (tyszler.jobs@gmail.com):


Tests if the wrapper reg_sandwich is calling the appropriate functions
******************************************************************* */



/**** Initial Set-up ****/
set more off
clear

capture log close
log using "wrapper_tests_log.txt", replace text
***********************************************
display "Tests if the wrapper reg_sandwich is calling the appropriate functions"
display "WARNING: absorb is not fully implemented yet"
display "The goal of theses tests is to check if the coefficients, se and dfs are being properly calculated"
display "To match aweights with R, R pop needs to be divided by 100"
display "DateTime: $S_DATE $S_TIME"

use "MortalityRates", replace

* filter: cause=="Motor Vehicle"
label list cause
keep if cause == 2

* model specification

local specification = "mrate legal beertaxa beerpercap winepercap i.year"
  
* unweighted
* ols_pooled <- lm(specification, data = MV_Mortality)
*coef_test(ols_pooled, vcov = "CR2", cluster = MV_Mortality$state)
* xi, noomit: reg `specification', nocons cluster(state)
xi, noomit: reg_sandwich `specification', nocons cluster(state)

* a-weighted
* wls_pooled <- lm(specification, weights = pop, data = MV_Mortality)
* coef_test(wls_pooled, vcov = "CR2", cluster = MV_Mortality$state, inverse_var = TRUE)
*xi, noomit: reg `specification' [aweight=pop], nocons cluster(state)
xi, noomit: reg_sandwich `specification' [aweight=pop], nocons cluster(state)

* p-weighted
* coef_test(wls_pooled, vcov = "CR2", cluster = MV_Mortality$state, inverse_var = FALSE)
*xi, noomit: reg `specification' [pweight=pop], nocons cluster(state)
xi, noomit: reg_sandwich `specification' [pweight=pop], nocons cluster(state)

/**
** with absorption: 
disp "WARNING: areg has no 'noconstant' option, therefore the values for the dummies are in a different level"

* unweighted
* ols_within <- plm(update(specification, . ~ . - 0), data = MV_Mortality, effect = "individual", index = c("state","year"))
*xi: areg `specification', absorb(state) cluster(state)
xi: reg_sandwich `specification', absorb(state) cluster(state)

* a-weighted
* wls_within <- lm(update(specification, . ~ . + factor(state)), weights = pop, data = MV_Mortality)
* coef_test(wls_within, vcov = "CR2", cluster = MV_Mortality$state, inverse_var = TRUE)
*xi, noomit: areg `specification' [aweight=pop], absorb(state) cluster(state)
xi, noomit: reg_sandwich `specification' [aweight=pop], absorb(state) cluster(state)

* p-weighted
* coef_test(wls_within, vcov = "CR2", cluster = MV_Mortality$state, inverse_var = FALSE, ignore_FE = TRUE)
*xi, noomit: areg `specification' [pweight=pop], absorb(state) cluster(state)
xi, noomit: reg_sandwich `specification' [pweight=pop], absorb(state) cluster(state)
**/
****************
log close
