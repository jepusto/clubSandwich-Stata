{smcl}
{* *! version 0  18jan2017}{...}

{title:Title}

{phang}
{bf:reg_sandwich} {hline 2}  Linear regression with clubSandwich standard errors, and small-sample t-tests for each coefficient. 

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd: reg_sandwich}
{depvar} [{indepvars}] {ifin} [{it:{help weight:weight}}] {cmd:,}
cluster({varname}) 
[absorb({varname}) | {cmdab:nocon:stant}] 



{synoptset 26 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Mandatory}

{synopt :{opth cluster:(varname)}} {varname} for clustered sandwich estimator. {p_end}

{syntab:Optional}

{synopt :{opth absorb:(varname)}} categorical variable to be absorbed. {p_end}
{synopt :{cmdab:nocon:stant}} suppress constant term. {p_end}
{synopt :*} {it:absorb and noconstant cannot be used simultaneously}. {p_end}

{synoptline}
INCLUDE help fvvarlist
{p 4 6 2}
{cmd:aweight}s and {cmd:pweight}s are
allowed; see {help weight}.{p_end}





{title:Description}

{pstd}
HERE COMES THE DETAILED DESCRIPTION.
{cmd: reg_sandwich} fits a linear regression using {help regress}, optionally passing aweights or pweights (see {help weight}). If {it:absorb} is provided regression are fitted using {help areg}.
{p_end}

{title:Arguments}

{dlgtab:Mandatory}

{pmore}
{opth cluster:(varname)} {varname} for clustered sandwich estimator. See {helpb vce_option:[R] {it:vce_option}}.
{p_end}

{dlgtab:Optional}
{pmore}
({it:These two optional cannot be used simultaneously})

{pmore}
{opth absorb:(varname)}} specifies the categorical variable, which is to be included in the regression as if it were specified by dummy variables. See {help areg}. {p_end}

{pmore}
{cmdab:nocon:stant} suppresses the constant term (intercept) in the model.{p_end}

{title:Examples}

{phang}{cmd:. use MortalityRates.dta}{p_end}
{phang}{cmd:. keep if cause == 2} {p_end}
{phang}{cmd:. xi, noomit: reg_sandwich mrate legal beertaxa beerpercap winepercap i.year, nocon cluster(state)}{p_end}
{phang}{cmd:. xi, noomit: reg_sandwich mrate legal beertaxa beerpercap winepercap i.year [aweight=pop], nocon cluster(state)}{p_end}
{phang}{cmd:. xi, noomit: reg_sandwich mrate legal beertaxa beerpercap winepercap i.year [pweight=pop], nocon cluster(state)}{p_end}
{phang}{cmd:. xi: reg_sandwich mrate legal beertaxa beerpercap winepercap i.year, nocon cluster(state) absorb(state)}{p_end}
{phang}{cmd:. xi, noomit: reg_sandwich mrate legal beertaxa beerpercap winepercap i.year [aweight=pop], nocon cluster(state) absorb(state)}{p_end}
{phang}{cmd:. xi, noomit: reg_sandwich mrate legal beertaxa beerpercap winepercap i.year [pweight=pop], nocon cluster(state) absorb(state)}{p_end}



{title:References}

{phang}
Github repository:  {browse "https://github.com/jepusto/clubSandwich-Stata"} {p_end}

{phang}
Elizabeth Tipton and James E. Pustejovsky, 2015. Small-sample adjustments for tests of moderators and model fit 
using robust variance estimation in meta-regression. Journal of Educational and Behavioral Statistics December 2015 vol. 40 no. 6 604-634. 
DOI: 10.3102/1076998615606099
{p_end}


