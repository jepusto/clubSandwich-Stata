// github repository: https://github.com/jepusto/clubSandwich-Stata
//
*! version 0.0 updated 17-Jan-2017
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
	
	
	** call regression:
	disp "`main_function' `varlist' `weight_call', `constant' cluster(`cluster') `absorb_call'"
	`main_function' `varlist' `weight_call', `constant' cluster(`cluster') `absorb_call'
	
end
