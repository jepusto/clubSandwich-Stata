set more off
clear
use http://masteringmetrics.com/wp-content/uploads/2015/01/deaths.dta
keep if dtype == 2 & agegr == 2
xi, noomit: reg_sandwich mrate legal beertaxa beerpercap winepercap i.year, nocon cluster(state)
xi, noomit: reg_sandwich mrate legal beertaxa beerpercap winepercap i.year [aweight=pop], nocon cluster(state)
xi, noomit: reg_sandwich mrate legal beertaxa beerpercap winepercap i.year [pweight=pop], nocon cluster(state)
xi: reg_sandwich mrate legal beertaxa beerpercap winepercap i.year, cluster(state) absorb(state)
xi, noomit: reg_sandwich mrate legal beertaxa beerpercap winepercap i.year [aweight=pop], cluster(state) absorb(state)
xi, noomit: reg_sandwich mrate legal beertaxa beerpercap winepercap i.year [pweight=pop], cluster(state) absorb(state)
