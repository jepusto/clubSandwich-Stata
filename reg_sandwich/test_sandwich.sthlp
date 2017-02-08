{smcl}
{* *! version 1.0.0  12sep2016}{...}
{findalias asfradohelp}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "[R] help" "help help"}{...}
{viewerjumpto "Syntax" "examplehelpfile##syntax"}{...}
{viewerjumpto "Description" "examplehelpfile##description"}{...}
{viewerjumpto "Options" "examplehelpfile##options"}{...}
{viewerjumpto "Remarks" "examplehelpfile##remarks"}{...}
{viewerjumpto "Examples" "examplehelpfile##examples"}{...}
{title:Title}

{phang}
{bf:test_sandwich} {hline 2} Computes small sample corrected F-test of parameters estimated by {help reg_sandwich}


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd: test_sandwich}
{varlist}
[{cmd:,} {it:cons}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt cons}} includes the constant term{p_end}
{synoptline}
{p2colreset}{...}

{marker description}{...}
{title:Description}

{pstd}
{cmd:test_sandwich} is a post-estimation command to be run after {help reg_sandwich}. 
It computes a small sample corrected F-test of parameters included in {varlist}, using the same sample and weights as in the {it:reg_sandwich} estimation


{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt cons} includes the constant term in the estimation.

{marker remarks}{...}
{title:Remarks}

{pstd}
Parameters included in {varlist} need to have been estimated by {it:reg_sandwich}.



{marker examples}{...}
{title:Examples}
{phang}{cmd:. use MortalityRates.dta}{p_end}
{phang}{cmd:. keep if cause == 2} {p_end}
{phang}{cmd:. xi, noomit: reg_sandwich mrate legal beertaxa beerpercap winepercap i.year, nocon cluster(state)}{p_end}
{phang}{cmd:. test_sandwich legal beertaxa}{p_end}


{title:Saved results}
{pstd}
{cmd:test_sandwich} saves the following in {cmd:e()}:
{synoptset 20 tabbed}{...}
{p_end}
{synopt:{cmd:e(F_stat)}} F statistic{p_end}
{synopt:{cmd:e(F_df1)}} degrees of freedom from the numerator{p_end}
{synopt:{cmd:e(F_df1)}} degrees of freedom from the denominator{p_end}
{synopt:{cmd:e(F_pvalue)}} p-value of the F statistic{p_end}
{synopt:{cmd:e(F_eta)}} {it:eta} computed for the F statistic{p_end}
