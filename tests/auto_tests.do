set more off
clear

use http://www.stata-press.com/data/r13/auto7


regress mpg weight foreign, cluster(manufacturer)

reg_sandwich mpg weight foreign, cluster(manufacturer)

test_sandwich weight foreign
