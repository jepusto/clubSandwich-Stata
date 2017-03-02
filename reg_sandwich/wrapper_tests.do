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
display "The goal of theses tests is to check if the coefficients, se and dfs are being properly calculated"
display "To match aweights with R, R pop needs to be divided by 100"
display "DateTime: $S_DATE $S_TIME"

use "MortalityRates", replace

* filter: cause=="Motor Vehicle"
label list cause
keep if cause == 2

* model specification

local specification = "mrate legal beertaxa beerpercap winepercap i.year"

* timer
timer clear
local i = 0 
  
* unweighted
* ols_pooled <- lm(specification, data = MV_Mortality)
*coef_test(ols_pooled, vcov = "CR2", cluster = MV_Mortality$state)
* xi, noomit: reg `specification', nocons cluster(state)

local i = `i'+1
timer on `i'
xi, noomit: reg_sandwich `specification', nocons cluster(state)
ereturn list
** Ftests:
timer off `i'
local i = `i'+1
timer on `i'
test_sandwich legal beertaxa

timer off `i'
local i = `i'+1
timer on `i'
test_sandwich legal beerpercap
timer off `i'
local i = `i'+1
timer on `i'
test_sandwich legal winepercap

timer off `i'
local i = `i'+1
timer on `i'
test_sandwich beertaxa beerpercap
timer off `i'
local i = `i'+1
timer on `i'
test_sandwich beertaxa winepercap

timer off `i'
local i = `i'+1
timer on `i'
test_sandwich beerpercap winepercap

timer off `i'
local i = `i'+1
timer on `i'
test_sandwich legal beertaxa beerpercap
timer off `i'
local i = `i'+1
timer on `i'
test_sandwich legal beertaxa winepercap
timer off `i'
local i = `i'+1
timer on `i'
test_sandwich legal beerpercap winepercap
timer off `i'
local i = `i'+1
timer on `i'
test_sandwich beertaxa beerpercap winepercap

timer off `i'
local i = `i'+1
timer on `i'
test_sandwich legal beertaxa beerpercap winepercap
timer off `i'
disp "Timers:"
timer list
timer clear
local i = 0 

* a-weighted
* wls_pooled <- lm(specification, weights = pop, data = MV_Mortality)
* coef_test(wls_pooled, vcov = "CR2", cluster = MV_Mortality$state, inverse_var = TRUE)
*xi, noomit: reg `specification' [aweight=pop], nocons cluster(state)
local i = `i'+1
timer on `i'
xi, noomit: reg_sandwich `specification' [aweight=pop], nocons cluster(state)
ereturn list

timer off `i'
local i = `i'+1
timer on `i'
test_sandwich legal beertaxa
timer off `i'
local i = `i'+1
timer on `i'
test_sandwich legal beerpercap
timer off `i'
local i = `i'+1
timer on `i'
test_sandwich legal winepercap

timer off `i'
local i = `i'+1
timer on `i'
test_sandwich beertaxa beerpercap
timer off `i'
local i = `i'+1
timer on `i'
test_sandwich beertaxa winepercap

timer off `i'
local i = `i'+1
timer on `i'
test_sandwich beerpercap winepercap

timer off `i'
local i = `i'+1
timer on `i'
test_sandwich legal beertaxa beerpercap
timer off `i'
local i = `i'+1
timer on `i'
test_sandwich legal beertaxa winepercap
timer off `i'
local i = `i'+1
timer on `i'
test_sandwich legal beerpercap winepercap
timer off `i'
local i = `i'+1
timer on `i'
test_sandwich beertaxa beerpercap winepercap

timer off `i'
local i = `i'+1
timer on `i'
test_sandwich legal beertaxa beerpercap winepercap
timer off `i'
disp "Timers:"
timer list
timer clear
local i = 0 

* p-weighted
* coef_test(wls_pooled, vcov = "CR2", cluster = MV_Mortality$state, inverse_var = FALSE)
*xi, noomit: reg `specification' [pweight=pop], nocons cluster(state)

local i = `i'+1
timer on `i'
xi, noomit: reg_sandwich `specification' [pweight=pop], nocons cluster(state)
ereturn list

timer off `i'
local i = `i'+1
timer on `i'
test_sandwich legal beertaxa
timer off `i'
local i = `i'+1
timer on `i'
test_sandwich legal beerpercap
timer off `i'
local i = `i'+1
timer on `i'
test_sandwich legal winepercap

timer off `i'
local i = `i'+1
timer on `i'
test_sandwich beertaxa beerpercap
timer off `i'
local i = `i'+1
timer on `i'
test_sandwich beertaxa winepercap

timer off `i'
local i = `i'+1
timer on `i'
test_sandwich beerpercap winepercap

timer off `i'
local i = `i'+1
timer on `i'
test_sandwich legal beertaxa beerpercap
timer off `i'
local i = `i'+1
timer on `i'
test_sandwich legal beertaxa winepercap
timer off `i'
local i = `i'+1
timer on `i'
test_sandwich legal beerpercap winepercap
timer off `i'
local i = `i'+1
timer on `i'
test_sandwich beertaxa beerpercap winepercap

timer off `i'
local i = `i'+1
timer on `i'
test_sandwich legal beertaxa beerpercap winepercap
timer off `i'
disp "Timers:"
timer list
timer clear
local i = 0 

** with absorption: 
disp "WARNING: areg has no 'noconstant' option, therefore the values for the dummies are in a different level"

* unweighted
* ols_within <- plm(update(specification, . ~ . - 0), data = MV_Mortality, effect = "individual", index = c("state","year"))
*xi: areg `specification', absorb(state) cluster(state)

* compare results with explicit dummies

local i = `i'+1

xi: reg_sandwich `specification' i.state,  cluster(state)
timer on `i'
xi: reg_sandwich `specification', absorb(state) cluster(state)
ereturn list

timer off `i'
local i = `i'+1
timer on `i'
test_sandwich legal beertaxa
timer off `i'
local i = `i'+1
timer on `i'
test_sandwich legal beerpercap
timer off `i'
local i = `i'+1
timer on `i'
test_sandwich legal winepercap

timer off `i'
local i = `i'+1
timer on `i'
test_sandwich beertaxa beerpercap
timer off `i'
local i = `i'+1
timer on `i'
test_sandwich beertaxa winepercap

timer off `i'
local i = `i'+1
timer on `i'
test_sandwich beerpercap winepercap

timer off `i'
local i = `i'+1
timer on `i'
test_sandwich legal beertaxa beerpercap
timer off `i'
local i = `i'+1
timer on `i'
test_sandwich legal beertaxa winepercap
timer off `i'
local i = `i'+1
timer on `i'
test_sandwich legal beerpercap winepercap
timer off `i'
local i = `i'+1
timer on `i'
test_sandwich beertaxa beerpercap winepercap

timer off `i'
local i = `i'+1
timer on `i'
test_sandwich legal beertaxa beerpercap winepercap
timer off `i'
disp "Timers:"
timer list
timer clear
local i = 0 

* a-weighted
* wls_within <- lm(update(specification, . ~ . + factor(state)), weights = pop, data = MV_Mortality)
* coef_test(wls_within, vcov = "CR2", cluster = MV_Mortality$state, inverse_var = TRUE)
*xi, noomit: areg `specification' [aweight=pop], absorb(state) cluster(state)

* compare results with explicit dummies

local i = `i'+1

xi, noomit: reg_sandwich `specification' i.state [aweight=pop], cluster(state)
timer on `i'
xi, noomit: reg_sandwich `specification' [aweight=pop], absorb(state) cluster(state)
ereturn list

timer off `i'
local i = `i'+1
timer on `i'
test_sandwich legal beertaxa
timer off `i'
local i = `i'+1
timer on `i'
test_sandwich legal beerpercap
timer off `i'
local i = `i'+1
timer on `i'
test_sandwich legal winepercap

timer off `i'
local i = `i'+1
timer on `i'
test_sandwich beertaxa beerpercap
timer off `i'
local i = `i'+1
timer on `i'
test_sandwich beertaxa winepercap

timer off `i'
local i = `i'+1
timer on `i'
test_sandwich beerpercap winepercap

timer off `i'
local i = `i'+1
timer on `i'
test_sandwich legal beertaxa beerpercap
timer off `i'
local i = `i'+1
timer on `i'
test_sandwich legal beertaxa winepercap
timer off `i'
local i = `i'+1
timer on `i'
test_sandwich legal beerpercap winepercap
timer off `i'
local i = `i'+1
timer on `i'
test_sandwich beertaxa beerpercap winepercap

timer off `i'
local i = `i'+1
timer on `i'
test_sandwich legal beertaxa beerpercap winepercap
timer off `i'
disp "Timers:"
timer list
timer clear
local i = 0 
* p-weighted
* coef_test(wls_within, vcov = "CR2", cluster = MV_Mortality$state, inverse_var = FALSE, ignore_FE = TRUE)
*xi, noomit: areg `specification' [pweight=pop], absorb(state) cluster(state)

* compare results with explicit dummies

local i = `i'+1

xi, noomit: reg_sandwich `specification' i.state [pweight=pop], cluster(state)
timer on `i'
xi, noomit: reg_sandwich `specification' [pweight=pop], absorb(state) cluster(state)
ereturn list

timer off `i'
local i = `i'+1
timer on `i'
test_sandwich legal beertaxa
timer off `i'
local i = `i'+1
timer on `i'
test_sandwich legal beerpercap
timer off `i'
local i = `i'+1
timer on `i'
test_sandwich legal winepercap

timer off `i'
local i = `i'+1
timer on `i'
test_sandwich beertaxa beerpercap
timer off `i'
local i = `i'+1
timer on `i'
test_sandwich beertaxa winepercap

timer off `i'
local i = `i'+1
timer on `i'
test_sandwich beerpercap winepercap

timer off `i'
local i = `i'+1
timer on `i'
test_sandwich legal beertaxa beerpercap
timer off `i'
local i = `i'+1
timer on `i'
test_sandwich legal beertaxa winepercap
timer off `i'
local i = `i'+1
timer on `i'
test_sandwich legal beerpercap winepercap
timer off `i'
local i = `i'+1
timer on `i'
test_sandwich beertaxa beerpercap winepercap

timer off `i'
local i = `i'+1
timer on `i'
test_sandwich legal beertaxa beerpercap winepercap
timer off `i'
local i = `i'+1
timer off `i'
disp "Timers:"
timer list
timer clear

****************
log close
