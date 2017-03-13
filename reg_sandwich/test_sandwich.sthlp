{smcl}
{* *! version 1.0.0  12sep2016}{...}

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

{p 4 6 2}
{cmd:test_sandwich} uses the supporting mata function
test_sandwich_ftests.mo{p_end}

{p 4 6 2}
see also {help reg_sandwich}{p_end}

{marker description}{...}
{title:Description}

{pstd}
{cmd:test_sandwich} is a post-estimation command to be run after {help reg_sandwich}. 
It computes a small sample corrected F-test of parameters included in {varlist}, using the same sample and weights as in the {it:reg_sandwich} estimation


{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt cons} includes the constant term in the F-test.

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

{title:Citation}
{phang}
{cmd:test_sandwich} is not an official Stata command. It is a free contribution to the research community, like a paper.
Please cite it as such:{p_end}

{phang}
Tyszler, M., Pustejovsky, J.E., Tipton, E. 2017. clubSandwich: Cluster-robust variance estimation and hypothesis testing with small-sample corrections for linear regression. 
URL: {browse "https://github.com/jepusto/clubSandwich-Stata"}
{p_end}

{title:Authors}
{phang} Marcelo Tyszler. Sustainable Economic Development and Gender, Royal Tropical Institute, Netherlands. {browse "mailto:m.tyszler@kit.nl":m.tyszler@kit.nl} {p_end}

{phang} James E. Pustejovsky {bf:{it: (Package maintainer)}}. Department of Education Psychology, University of Texas at Austin. {browse "mailto:pusto@austin.utexas.edu":pusto@austin.utexas.edu}{p_end}

{phang} Elizabeth Tipton. Department of Human Development, Teachers College, Columbia University. {browse "mailto:tipton@tc.columbia.edu":tipton@tc.columbia.edu} {p_end}


{title:References}
{phang}
Pustejovsky, James E. & Elizabeth Tipton (2016). 
Small sample methods for cluster-robust variance estimation and hypothesis testing in fixed effects models. 
Journal of Business and Economic Statistics. In Press. DOI: 10.1080/07350015.2016.1247004
{p_end}

{phang}
Github repository:  {browse "https://github.com/jepusto/clubSandwich-Stata"} {p_end}

{phang}
Elizabeth Tipton and James E. Pustejovsky, 2015. Small-sample adjustments for tests of moderators and model fit 
using robust variance estimation in meta-regression. Journal of Educational and Behavioral Statistics December 2015 vol. 40 no. 6 604-634. 
DOI: 10.3102/1076998615606099
{p_end}

{phang}
Bell, R. M., & McCaffrey, D. F. (2002). Bias reduction in standard errors for linear regression with multi-stage samples. 
Survey Methodology, 28(2), 169–181. 
Retrieved from {browse "http://www.statcan.gc.ca/pub/12-001-x/2002002/article/9058-eng.pdf"}
{p_end}


