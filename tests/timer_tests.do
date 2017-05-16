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
log using "timer_tests_log.txt", replace text
***********************************************
display "Tests timer"

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
/** Ftests:
test_sandwich legal beertaxa
test_sandwich legal beerpercap
test_sandwich legal winepercap

test_sandwich beertaxa beerpercap
test_sandwich beertaxa winepercap

test_sandwich beerpercap winepercap

test_sandwich legal beertaxa beerpercap
test_sandwich legal beertaxa winepercap
test_sandwich legal beerpercap winepercap
test_sandwich beertaxa beerpercap winepercap

test_sandwich legal beertaxa beerpercap winepercap


**/

* a-weighted
* wls_pooled <- lm(specification, weights = pop, data = MV_Mortality)
* coef_test(wls_pooled, vcov = "CR2", cluster = MV_Mortality$state, inverse_var = TRUE)
*xi, noomit: reg `specification' [aweight=pop], nocons cluster(state)
xi, noomit: reg_sandwich `specification' [aweight=pop], nocons cluster(state)
/**test_sandwich legal beertaxa
test_sandwich legal beerpercap
test_sandwich legal winepercap

test_sandwich beertaxa beerpercap
test_sandwich beertaxa winepercap

test_sandwich beerpercap winepercap

test_sandwich legal beertaxa beerpercap
test_sandwich legal beertaxa winepercap
test_sandwich legal beerpercap winepercap
test_sandwich beertaxa beerpercap winepercap

test_sandwich legal beertaxa beerpercap winepercap
**/

* p-weighted
* coef_test(wls_pooled, vcov = "CR2", cluster = MV_Mortality$state, inverse_var = FALSE)
*xi, noomit: reg `specification' [pweight=pop], nocons cluster(state)
xi, noomit: reg_sandwich `specification' [pweight=pop], nocons cluster(state)
/**test_sandwich legal beertaxa
test_sandwich legal beerpercap
test_sandwich legal winepercap

test_sandwich beertaxa beerpercap
test_sandwich beertaxa winepercap

test_sandwich beerpercap winepercap

test_sandwich legal beertaxa beerpercap
test_sandwich legal beertaxa winepercap
test_sandwich legal beerpercap winepercap
test_sandwich beertaxa beerpercap winepercap

test_sandwich legal beertaxa beerpercap winepercap

**/
** with absorption: 
disp "WARNING: areg has no 'noconstant' option, therefore the values for the dummies are in a different level"

* unweighted
* ols_within <- plm(update(specification, . ~ . - 0), data = MV_Mortality, effect = "individual", index = c("state","year"))
*xi: areg `specification', absorb(state) cluster(state)

* compare results with explicit dummies
*xi: reg_sandwich `specification' i.state,  cluster(state)
xi: reg_sandwich `specification', absorb(state) cluster(state)
/**
test_sandwich legal beertaxa
test_sandwich legal beerpercap
test_sandwich legal winepercap

test_sandwich beertaxa beerpercap
test_sandwich beertaxa winepercap

test_sandwich beerpercap winepercap

test_sandwich legal beertaxa beerpercap
test_sandwich legal beertaxa winepercap
test_sandwich legal beerpercap winepercap
test_sandwich beertaxa beerpercap winepercap

test_sandwich legal beertaxa beerpercap winepercap
**/

* a-weighted
* wls_within <- lm(update(specification, . ~ . + factor(state)), weights = pop, data = MV_Mortality)
* coef_test(wls_within, vcov = "CR2", cluster = MV_Mortality$state, inverse_var = TRUE)
*xi, noomit: areg `specification' [aweight=pop], absorb(state) cluster(state)

* compare results with explicit dummies
*xi, noomit: reg_sandwich `specification' i.state [aweight=pop], cluster(state)
xi, noomit: reg_sandwich `specification' [aweight=pop], absorb(state) cluster(state)
/**
test_sandwich legal beertaxa
test_sandwich legal beerpercap
test_sandwich legal winepercap

test_sandwich beertaxa beerpercap
test_sandwich beertaxa winepercap

test_sandwich beerpercap winepercap

test_sandwich legal beertaxa beerpercap
test_sandwich legal beertaxa winepercap
test_sandwich legal beerpercap winepercap
test_sandwich beertaxa beerpercap winepercap

test_sandwich legal beertaxa beerpercap winepercap
**/
* p-weighted
* coef_test(wls_within, vcov = "CR2", cluster = MV_Mortality$state, inverse_var = FALSE, ignore_FE = TRUE)
*xi, noomit: areg `specification' [pweight=pop], absorb(state) cluster(state)

* compare results with explicit dummies
*xi, noomit: reg_sandwich `specification' i.state [pweight=pop], cluster(state)
xi, noomit: reg_sandwich `specification' [pweight=pop], absorb(state) cluster(state)
/**
test_sandwich legal beertaxa
test_sandwich legal beerpercap
test_sandwich legal winepercap

test_sandwich beertaxa beerpercap
test_sandwich beertaxa winepercap

test_sandwich beerpercap winepercap

test_sandwich legal beertaxa beerpercap
test_sandwich legal beertaxa winepercap
test_sandwich legal beerpercap winepercap
test_sandwich beertaxa beerpercap winepercap

test_sandwich legal beertaxa beerpercap winepercap
**/
****************
log close
