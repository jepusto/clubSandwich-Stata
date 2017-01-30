// github repository: https://github.com/jepusto/clubSandwich-Stata
//
*! version 0.0 updated 30-Jan-2017
// Updated by Marcelo Tyszler (tyszler.jobs@gmail.com):
//
// Initialize the wrapper function for reg_sandwich:
//
// allow for factor variables via xi:.
// have an option to absorb fixed effects, as in the areg command.
// work with pweights and aweights.
// disregard absorbed fixed effects when calculating the adjustment matrices, degrees of freedom, etc.
// save the abjustment matrices and any other required information for calculating F-tests based on the fitted model.
//
// Suggested syntax is as follows:
// reg_sandwich depvar [indepvars] [if] [in] [weight], cluster(varname) [absorb(varname)]
// 
// If the absorb option is present, then the model should be fit as in areg. 
// If the option is absent, then the model should be fit as in reg. 
// If weights are included, the model should be fit via WLS. 

capture program drop reg_sandwich
program define reg_sandwich, eclass sortpreserve

	version 14.2 
	syntax varlist(min=1 numeric fv) [if] [in] ///
	[aweight pweight],  ///
    cluster(varlist max=1 numeric) ///
	[absorb(varlist max=1 numeric)] ///
	[noCONstant]
	
	*mark sample
    marksample touse
	
	*create macros of variables
    tokenize `varlist'
    local t `1'
    macro shift
    local x `*'
	** determine main function
	capture confirm existence `absorb'
	if _rc == 6{
		* no absorb, use default reg
		local main_function = "reg"
		local absorb_call = ""
	} 
	else {
		* absorb, use areg
		local main_function = "areg"
		local absorb_call = "absorb(`absorb')"
		if "`constant'"!="" {
			di as error "absorb and noconstant cannot be used simultaneously"
			exit 198
		}
	} 
    
	** determine weights
	capture confirm existence `weight'
	if _rc == 6{
		* no weights
		local weight_call = ""
	} 
	else {
		* weights, use areg
		local weight_call = "[`weight'`exp']"
	}
	
	*collinearity 
    local olist "`x'"
	if "`constant'"=="" {
		_rmcoll `x' if `touse'
	}
	else {
		_rmcoll `x' if `touse', noconstant
	}
	
    local x = r(varlist)
    foreach v in `olist' {
        local x = regexr("`x'","o\.`v'","")
    }

    if "`x'" == "." local x ""
	
	
	
	** call regression:
	*disp "`main_function' `t' `x' `weight_call' if `touse', `constant' cluster(`cluster') `absorb_call'"
	noisily capture: `main_function' `t' `x'  `weight_call' if `touse', `constant' cluster(`cluster') `absorb_call'
	
	
end
