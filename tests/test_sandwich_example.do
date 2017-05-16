set more off
clear
use http://masteringmetrics.com/wp-content/uploads/2015/01/deaths.dta
keep if dtype == 2 & agegr == 2
xi: reg_sandwich mrate legal beertaxa beerpercap winepercap i.year, cluster(state) absorb(state)
test_sandwich legal beertaxa
