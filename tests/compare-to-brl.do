set more off
clear
use  http://stats.idre.ucla.edu/stat/stata/seminars/svy_stata_intro/srs, clear

regress api00 growth emer yr_rnd, cluster(dnum)
brl api00 growth emer yr_rnd, cluster(dnum)
reg_sandwich api00 growth emer yr_rnd, cluster(dnum)

clear
use http://masteringmetrics.com/wp-content/uploads/2015/01/deaths.dta
keep if dtype == 2 & agegr == 2

* no fixed effects
xi, noomit: regress mrate legal beertaxa beerpercap winepercap, cluster(state)
xi, noomit: reg_sandwich mrate legal beertaxa beerpercap winepercap, cluster(state)
xi, noomit: brl mrate legal beertaxa beerpercap winepercap, cluster(state)

* year fixed effects
xi, noomit: regress mrate legal beertaxa beerpercap winepercap i.year, cluster(state)
xi, noomit: reg_sandwich mrate legal beertaxa beerpercap winepercap i.year, cluster(state)
xi, noomit: brl mrate legal beertaxa beerpercap winepercap i.year, cluster(state)

* year and state fixed effects
xi, noomit: regress mrate legal beertaxa beerpercap winepercap i.year i.state, cluster(state)
xi, noomit: reg_sandwich mrate legal beertaxa beerpercap winepercap i.year i.state, cluster(state)
xi, noomit: brl mrate legal beertaxa beerpercap winepercap i.year i.state, cluster(state)

