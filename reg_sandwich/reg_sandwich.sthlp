{smcl}
{* *! version 0  18jan2017}{...}

{title:Title}

{phang}
{bf:reg_sandwich} {hline 2}  Linear regression with clubSandwich standard errors and small-sample t-tests for each coefficient. 

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd: reg_sandwich}
{depvar} [{indepvars}] {ifin} [{it:{help weight:weight}}]{cmd:,}
cluster({varname}) 
[absorb({varname}) | {cmdab:nocon:stant}] 
[level({#})]




{synoptset 26 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Mandatory}

{synopt :{opth cluster:(varname)}} {varname} for clustered sandwich estimator. {p_end}

{syntab:Optional}

{synopt :{opth absorb:(varname)}} categorical variable to be absorbed. {p_end}
{synopt :{cmdab:nocon:stant}} suppress constant term. {p_end}
{synopt :*} {it:absorb and noconstant cannot be used simultaneously}. {p_end}
{synopt : level(#)} set confidence level; default is level(95). {p_end}

{synoptline}

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

{dlgtab:Reporting}

{phang}
{opt level(#)}; see {helpb estimation options##level():[R] estimation options}.

{title:Examples}

{phang}{cmd:. use MortalityRates.dta}{p_end}
{phang}{cmd:. keep if cause == 2} {p_end}
{phang}{cmd:. xi, noomit: reg_sandwich mrate legal beertaxa beerpercap winepercap i.year, nocon cluster(state)}{p_end}
{phang}{cmd:. xi, noomit: reg_sandwich mrate legal beertaxa beerpercap winepercap i.year [aweight=pop], nocon cluster(state)}{p_end}
{phang}{cmd:. xi, noomit: reg_sandwich mrate legal beertaxa beerpercap winepercap i.year [pweight=pop], nocon cluster(state)}{p_end}
{phang}{cmd:. xi: reg_sandwich mrate legal beertaxa beerpercap winepercap i.year, cluster(state) absorb(state)}{p_end}
{phang}{cmd:. xi, noomit: reg_sandwich mrate legal beertaxa beerpercap winepercap i.year [aweight=pop], cluster(state) absorb(state)}{p_end}
{phang}{cmd:. xi, noomit: reg_sandwich mrate legal beertaxa beerpercap winepercap i.year [pweight=pop], cluster(state) absorb(state)}{p_end}

{title:Saved results}

{pstd}
{cmd:re_sandwich} saves the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(N_g)}}number of studies{p_end}
{synopt:{cmd:e(df_r)}}model degrees of freedom{p_end}

{synopt:{cmd:e(tau2)}} method-of-moments tau-square estimate {p_end}
{synopt:{cmd:e(tau2o)}} observed tau-square if estimate is negative{p_end}

{synopt:{cmd:e(omega2)}} method-of-moments omega-square estimate (used in hierarchical model) {p_end}
{synopt:{cmd:e(omega2o)}} observed omega-square if estimate is negative{p_end}

{synopt:{cmd:e(QE)}}QE used for estimating tau-square {p_end}
{synopt:{cmd:e(QR)}}QR used for estimating omega-square {p_end}

{synopt:{cmd:e(rho)}} in correlated effects models, use specified ICC{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:robumeta}{p_end}
{synopt:{cmd:e(depvar)}}{depvar}{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of the estimators{p_end}
{synopt:{cmd:e(dfs)}}Degrees of freedom for effects{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}



{title:References}

{phang}
Github repository:  {browse "https://github.com/jepusto/clubSandwich-Stata"} {p_end}

{phang}
Elizabeth Tipton and James E. Pustejovsky, 2015. Small-sample adjustments for tests of moderators and model fit 
using robust variance estimation in meta-regression. Journal of Educational and Behavioral Statistics December 2015 vol. 40 no. 6 604-634. 
DOI: 10.3102/1076998615606099
{p_end}


