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
	[noCONstant] ///
	[Level(cilevel)]
	
	*mark sample
    marksample touse
	
	*create macros of variables
    tokenize `varlist'
    local t `1'
    macro shift
    local x `*'
	
	* Count valid observations and check matsize
	qui count if `touse'
    local nobs = r(N)
	if c(matsize)<`nobs' {
		set matsize `nobs'
	}
	
   *specify the temporary variables used in the program
    tempvar keepers cons w wh wfinal prelim_hat prelim_resid  ///
    prime_hat prime_resid v_mean v_n clusterweight clusternumber variance
   
    *specifiy the temporary scalars and matrixes
    tempname A1  A2  b   B1  B2  C1  C2  D   e   E   F   I   J    ///
        k   kw  kXJWX   kXWJX   omega_squared   omega_squared_o      ///
        Q1  QE  QR  sigmahat  min_n max_n sumk    sumk2   T   T_XB    T_XBJT_XB    ///
        tau_squared tau_squared_o trW trW_1  TWT     TWX     V   VkXJWX  VkXWJX   ///
        Vw2XJX  VwkXJX_XX   VwkXX   VXJWX   VXJX    VXJXVXW2X    ///
        VXJXVXWJWX  VXW2X   VXW2X   VXWJWX  VXWJWX  VXWJX   W    ///
        w_k     w2  w2XJX   wkXJX_XX    wkXX    wXX     X   XB   ///
        XJWX    XJX     XJX     XW2X    XW2X    XWeeWX  XWJWX   XWJX     ///
        XWT     XWX     XX ///
		M ///
		MXWVWXM ///
		XWAeeAWX Big_B_relevant Big_BSigmaB_relevant Big_VV ///
		Aj Wj Xj middle_Aj ej ///
		sq_Wj ///
		Vj ///
		make_it_fit ///
		_dfs middle_Omega ///
		cluster_list ///
		middle_BSigmaB BSigmaB ///
		C_ttest ///
		Omega_ttest matrix_ttest ///
		temp_calc ///
		prob 

	
	*generate constant term
	if "`constant'"=="" {
		quietly : gen double `cons' = 1 if `touse'
	}

	
	** determine main function
	capture confirm existence `absorb'
	if _rc == 6{
		* no absorb, use default reg
		local main_function = "reg"
		local absorb_call = ""
		local absorb_call_display = ""
	} 
	else {
		* absorb, use areg
		local main_function = "areg"
		local absorb_call = "absorb(`absorb')"
		local absorb_call_display = ", with absorb option"
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
		local main_call_display = "OLS"
	} 
	else {
		* weights
		local weight_call = "[`weight'`exp']"
		local main_call_display = "WLS (`weight's)"
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
	
	
	** prep for small sample reduced t-test:
	* (based on robumetaMt12.ado:
	

	
	matrix p = rowsof(e(V))
	local p = p[1,1]
	
	
	*capture ids
    quietly : gen double `clusternumber' = `cluster' if `touse'
    quietly sort `clusternumber'
    quietly : levelsof `clusternumber', local(idlist)
	
	*count ids, create macro m, number of studies
    local m = 0
    foreach j in `idlist' {
        local ++m
    }
	
	* double check set matsize:
	if `nobs'<`m'*`p'{
		local temp = `m'*`p'
		set matsize `temp'
	}
	
	
	*weights & variance:
	capture confirm existence `weight'
	if _rc == 6{
		* no weights = OLS
		quietly : gen double `wfinal' = 1 if `touse'
		* Working variance is I
		qui: gen double `variance' = 1 if `touse'
	} 
	else {
		* weights = WLS
		local model_weights = substr("`exp'",2,.)
		quietly : gen double `wfinal' = `model_weights' if `touse'
		
		if "`weight'"=="aweight"{			
			qui: gen double `variance' = 1/`model_weights'  if `touse'
		}
		else{
		* p-weights
		* Working variance is I
			qui: gen double `variance' = 1 if `touse'
		}
		
	
	}
	
	*cluster average variance
	*******
	quietly : by `clusternumber', sort rc0: egen double `v_mean' = mean(`variance') if `touse'
	*number of cases per cluster
	quietly : by `clusternumber', sort rc0: egen double `v_n' = count(`variance') if `touse'
				quietly sum `v_n' if `touse'
				scalar `min_n' = r(min) 
				scalar `max_n' = r(max)
	

	
	if "`constant'"=="" {
		mkmat `x' `cons' if `touse', matrix(`X')
		matrix colnames `X' = `x' _cons
	}
	else{
		mkmat `x'  if `touse', matrix(`X')
		matrix colnames `X' = `x'
	}
	
	mkmat `wfinal' if `touse', matrix(`W')
	matrix `W' = diag(`W')
	
	matrix `M' = invsym(`X'' * `W' * `X')
	
	mkmat `v_mean'   if `touse', matrix(`V')
	
	matrix `V' = diag(`V')
	
	matrix `MXWVWXM' =  `M'*`X''*`W'*`V'*`W'*`X'*`M'
	matrix drop `V' 
	
	/********************************************************************/
    /*    Variance covariance matrix estimation for standard errors     */
	/* 																    */
	/*    And F-test												    */
    /********************************************************************/
	qui: predict `prime_resid', residuals
	
	matrix `XWAeeAWX' = J(`p', `p', 0)
		
	local current_jcountFtest = 0
	
    foreach j in `idlist' {
        *Aj parameter 
		
		mkmat `wfinal' if `touse' & `clusternumber' == `j', matrix(`Wj')
		matrix `Wj' = diag(`Wj')  
		
		if "`constant'"=="" {
			mkmat `x' `cons' if `touse' & `clusternumber' == `j', matrix(`Xj')
			matrix colnames `Xj' = `x' _cons
		} 
		else {
			mkmat `x'  if `touse' & `clusternumber' == `j', matrix(`Xj')
			matrix colnames `Xj' = `x'
		}
				
		
		* we use that 
		* (I-X*M*X'*W)j*V*(I-X*M*X'*W)j' = 
		*
		* Vj - Vj*(Wj*Xj*M*Xj') - (Xj*M*Xj'*W)*Vj + Xj*(M*X'*W*V*W*X*M)*Xj'
			
		sum `v_mean' if `touse' & `clusternumber' == `j', meanonly
		local meanweight = r(mean)
		
		mkmat `v_mean' if `touse' & `clusternumber' == `j', matrix(`Vj')
		matrix `Vj' = diag(`Vj')  
		
		matrix `middle_Aj'=`Vj'-`Vj'*`Wj'*`Xj'*`M'*`Xj''-`Xj'*`M'*`Xj''*`Wj'*`Vj'+ `Xj'*`MXWVWXM'*`Xj''	
		
		matsqrt `middle_Aj'
											
		matrix `Aj' = (`meanweight'^.5)*inv(sq_`middle_Aj')
		matrix drop sq_`middle_Aj'												
		
			
        mkmat `prime_resid' if `touse' & `clusternumber' == `j', matrix(`ej')
		matrix `XWAeeAWX' = (`Xj'' * `Wj' * `Aj' * `ej' * `ej'' * `Aj' * `Wj' * `Xj') + `XWAeeAWX'

		
		* F-test:
		* Bj are defined in equation (9):
		* Bj = Omega^(-1/2)*C*M*Xj*Wj*Aj*(Ik - X*M*X'*W)j
		*
		* and Aj = Wj^(-1/2)*[Wj^(-1/2)*(inv(Wj) - Xj*M*Xj')*Wj^(-1/2)]^(-1/2)*Wj^(-1/2)
		*
		* Therefore to compute Bj we need C, Wj and Xj
		* We use also that:
		*
		*(I-X*M*X'*W)j*inv(W)*(I-X*M*X'*W)j' = inv(Wj) - Xj*M*Xj' if i == j
		*(I-X*M*X'*W)i*inv(W)*(I-X*M*X'*W)j' = -Xi*M*Xj' if i != j
		* 
		*[ Notice we get to: (I)i*inv(W)*(j)
		
		*
		* So we need: 
		* Bi*Sigma*Bj'=
		* inv(sq_Omega)*C*M*Xi'*Wi*Ai*(I-X*M*X'*W)i*inv(W)*(I-X*M*X'*W)j'*Aj'*Wj*Xj*M'*C*inv(sq_Omega)
		* 
		* if i == j:
		* = inv(sq_Omega)*C*M*Xj'*Wj*Aj*(inv(Wj) - Xj*M*Xj')*Aj'*Wj*Xj*M'*C*inv(sq_Omega) 
		* We save just the "middle" portion, which is independent of C:
		* M*Xj'*Wj*Aj*(inv(Wj) - Xi*M*Xj')*Aj'*Wj*Xj*M'
		* CALL IT Bj_Sigma_Bj_relevant
		*
		* if 1!=j:
		* inv(sq_Omega)*C*M*Xi'*Wi*Ai*(-Xi*M*Xj')*Aj'*Wj'*Xj'*M'*C'*inv(sq_Omega)' if i!=j
		* We save just the "middle" portion, which is independent of C:
		* Call Bj_relevant: 
		* M*Xj'*Wj*Aj*Xj (and ignore the (min) sign, since it will be cancelled out)
		
		
		local current_jcountFtest = `current_jcountFtest'+1
		
		tempname B`current_jcountFtest'_relevant  B`current_jcountFtest'_Sigma_B`current_jcountFtest'_relevant
		
		* for this weighitng type, we use V instead of inv(W) for Sigma
		* we use the fact that 
		* (I-X*M*X'*W)i*V*(I-X*M*X'*W)j' = 
		*
		* if i==j
		* Vj - Vj*(Wj*Xj*M*Xj') - (Xj*M*Xj'*W)*Vj + Xj*(M*X'*W*V*W*X*M)*Xj' = middle_Aj
		*
		* if i!=j
		* - Vi*Wi*Xi*M*Xj'   - Xi*M*Xj'*Wj*Vj     + Xi*(M*X'*W*V*W*X*M)*Xj'
		* we call VVj = Vj*Wj*Xj*M
		tempname VV`current_jcountFtest' 
		
		matrix `VV`current_jcountFtest'' = `Vj'*`Wj'*`Xj'*`M' // kj x p
		
		matrix `B`current_jcountFtest'_relevant' =  `M'*`Xj''*`Wj'*`Aj' // p x kj
		
		matrix `B`current_jcountFtest'_Sigma_B`current_jcountFtest'_relevant' = ///
													`M'*`Xj''*`Wj'*`Aj'* ///
													(`middle_Aj')* ///
													`Aj''*`Wj'*`Xj'*`M''
						
		matrix drop `Vj'
	
		
		* save for later
		if `current_jcountFtest'==1 {
			matrix `Big_BSigmaB_relevant' = `B`current_jcountFtest'_Sigma_B`current_jcountFtest'_relevant'

			if 	colsof(`B`current_jcountFtest'_relevant')<`max_n'{
				
				matrix `make_it_fit' = I(`max_n')
				matrix `make_it_fit' = `make_it_fit'[1..colsof(`B`current_jcountFtest'_relevant' ), 1..`max_n']
				
				matrix `Big_VV' = `VV`current_jcountFtest'''*`make_it_fit' // store it transposed
				matrix `Big_B_relevant' = `B`current_jcountFtest'_relevant'*`make_it_fit'
								
				matrix drop `make_it_fit'
			}
			else {
				matrix `Big_VV' = `VV`current_jcountFtest''' // store it transposed
				matrix `Big_B_relevant' = `B`current_jcountFtest'_relevant'
			}
			
		}
		else {
			matrix `Big_BSigmaB_relevant' = [`Big_BSigmaB_relevant' \ `B`current_jcountFtest'_Sigma_B`current_jcountFtest'_relevant']
			if 	colsof(`B`current_jcountFtest'_relevant')<`max_n'{
			
				matrix `make_it_fit' = I(`max_n')
				matrix `make_it_fit' = `make_it_fit'[1..colsof(`B`current_jcountFtest'_relevant' ), 1..`max_n']
									
				matrix `Big_VV' = [`Big_VV' \ `VV`current_jcountFtest'''*`make_it_fit'] // store it transposed
				matrix `Big_B_relevant' = [`Big_B_relevant' \ `B`current_jcountFtest'_relevant'*`make_it_fit']
				
				matrix drop `make_it_fit'
			}
			else {
				matrix `Big_VV' = [`Big_VV' \ `VV`current_jcountFtest'''] // store it transposed
				matrix `Big_B_relevant' = [`Big_B_relevant' \ `B`current_jcountFtest'_relevant']
				
			}
			 
			
		}	
    }
	matrix drop `Aj' `Wj' `Xj' `middle_Aj'  `ej'
		
	* RVE estimator
	matrix `V' = `M' * `XWAeeAWX' * `M'
	
	* T-test, using as a special case of an F-test:

	matrix `_dfs' =  J(1,`p', 0) 
	
	qui: tab `clusternumber' if `touse', matrow(`cluster_list')
	matrix `middle_Omega' = `MXWVWXM' // Full Variance
	
	forvalues i = 1/`m'{
		tempname X`i' 
		if "`constant'"=="" {
			mkmat `x' `cons' if `touse' & `clusternumber' == `cluster_list'[`i',1], matrix(`X`i'')
		}
		else {
			mkmat `x'  if `touse' & `clusternumber' == `cluster_list'[`i',1], matrix(`X`i'')
		}
	}
		
	

	
	
	
	forvalues i = 1/`m'{
		* We use the symmetry here, since that temp(i,j) =temp(j,i)
		forvalues j = `i'/`m'{

			if `i' == `j'{
	
				matrix  `BSigmaB' = `B`i'_Sigma_B`i'_relevant'
					
			}
			else {
				
				
			* for this weighitng type, we need to use V for Sigma
			* inv(sq_Omega)*C*M*Xi'*Wi*Ai*(I-X*M*X'*W)i*V*(I-X*M*X'*W)j'*Aj'*Wj*Xj*M'*C*inv(sq_Omega)
			* we use the fact that 
			* (I-X*M*X'*W)i*V*(I-X*M*X'*W)j' = 
			*
			* if i!=j
			* - Vi*Wi*Xi*M*Xj'   - Xi*M*Xj'*Wj*Vj     + Xi*(M*X'*W*V*W*X*M)*Xj'
			* we call VVj = Vj*Wj*Xj*M

				matrix `middle_BSigmaB' = -`VV`i''*`X`j'''-`X`i''*`VV`j''' + `X`i''*`middle_Omega'*`X`j'''
				
				
				
				
				matrix  `BSigmaB' = `B`i'_relevant'*`middle_BSigmaB'*`B`j'_relevant''
				matrix drop `middle_BSigmaB'
			}
			
			
			forvalues coefficient = 1/`p' {
				
				* prep C
					
				* * Define C
				matrix `C_ttest' = J(1, `p', 0) // C is initialized as a 1 x p matrix of zeros
				
				matrix `C_ttest'[1,`coefficient'] = 1
				
				
				
				matrix `Omega_ttest' = `C_ttest'*`middle_Omega'*`C_ttest''	
				
				local temp_val = `Omega_ttest'[1,1]
				matrix `matrix_ttest' = (1/sqrt(`temp_val'))*`C_ttest'
				matrix `temp_calc' = `matrix_ttest'*`BSigmaB'*`matrix_ttest''
				local  temp_calc2 = 2*((`temp_calc'[1,1])^2)
				
				* We use the symmetry here, since that temp(i,j) =temp(j,i)
				if `i'==`j'{
					matrix `_dfs'[1,`coefficient'] = `_dfs'[1,`coefficient'] + `temp_calc2'
					
				}
				else {
					matrix `_dfs'[1,`coefficient'] = `_dfs'[1,`coefficient'] + 2*`temp_calc2'
				}
				
			}
			matrix drop `BSigmaB'
			
		} 
		matrix drop `X`i''
		
	}
	
	forvalues coefficient = 1/`p' {
		matrix `_dfs'[1,`coefficient'] = 2/`_dfs'[1,`coefficient']
	}
		
	
     
    display _newline
    display as text "Robust standard error estimation using " as result "`main_call_display'`absorb_call_display'"


    *name the rows and columns of the matrixes

	if "`constant'"=="" {	
		matrix colnames `V' = `x' _cons
		matrix rownames `V' = `x' _cons
		matrix colnames `_dfs' = `x' _cons
	}
	else {
		matrix colnames `V' = `x' 
		matrix rownames `V' = `x' 
		matrix colnames `_dfs' = `x' 
	}
	

    
    /*********************/
    /*  Display results  */
    /*********************/

    display _col(55) as text "N Level 1" _col(69) "=" _col(69) as result %9.0f `nobs'
    display _col(55) as text "N Level 2" _col(69) "=" _col(69) as result %9.0f `m'
    display _col(55) as text "Min Level 1 n" _col(69) "=" _col(69) as result %9.0f `min_n'
    display _col(55) as text "Max Level 1 n" _col(69) "=" _col(69) as result %9.0f `max_n'
    display _col(55) as text _col (5) "Average" _col(69) "=" _col(69) as result  %9.2f `nobs' / `m'

    
    display as text  "{hline 13}" "{c TT}" "{hline 64}"

    display %12s abbrev("`t'",12)   _col(14) "{c |}" ///
                                    _col(21) "Coef." ///
                                    _col(29) "Std. Err." ///
                                    _col(40) "dfs" ///
                                    _col(50) "p-value" ///
                                    _col(60) "[" `level' "%Conf. Interval]"

    display as text  "{hline 13}" "{c +}" "{hline 64}"                            

    scalar `prob' = 0
	tempname effect variance dof
    local i = 1
	matrix `b' = e(b)
    foreach v in `x' {
        scalar `effect' = `b'[1,`i']
        scalar `variance' = `V'[`i',`i']
        scalar `dof' = `_dfs'[1,`i']

        if `dof' < 4 {
            local problem "!"
            scalar `prob' = 1
        }
        else {
            local problem ""
        }

        display %12s abbrev("`v'",12)   _col(14) "{c |}" ///
                                        _col(16) "`problem'" ///
                                        _col(21) %5.3f `effect' ///
                                        _col(29) %5.2f sqrt(`variance') ///
                                        _col(40) %5.2f `dof' ///
                                        _col(50) %5.4f 2*ttail(`dof',abs(`effect'/sqrt(`variance'))) ///
                                        _col(60) %5.4f `effect' - invttail(`dof',((100-`level')/100)/2)*sqrt(`variance') ///
                                        _col(70) %5.4f `effect' + invttail(`dof',((100-`level')/100)/2)*sqrt(`variance')
        local ++i
    }
	
	if "`constant'"=="" {
		local v = "_cons"
	    scalar `effect' = `b'[1,`i']
        scalar `variance' = `V'[`i',`i']
        scalar `dof' = `_dfs'[1,`i']

        if `dof' < 4 {
            local problem "!"
            scalar `prob' = 1
        }
        else {
            local problem ""
        }

        display %12s abbrev("`v'",12)   _col(14) "{c |}" ///
                                        _col(16) "`problem'" ///
                                        _col(21) %5.3f `effect' ///
                                        _col(29) %5.2f sqrt(`variance') ///
                                        _col(40) %5.2f `dof' ///
                                        _col(50) %5.4f 2*ttail(`dof',abs(`effect'/sqrt(`variance'))) ///
                                        _col(60) %5.4f `effect' - invttail(`dof',((100-`level')/100)/2)*sqrt(`variance') ///
                                        _col(70) %5.4f `effect' + invttail(`dof',((100-`level')/100)/2)*sqrt(`variance')
        local ++i
    }

    display as text  "{hline 13}" "{c BT}" "{hline 64}" 

    if `prob' == 1 {
        di as error "! dof is less than 4, p-value untrustworthy"
        di as error "see Tipton, E. (in press) Small sample adjustments for robust variance"
        di as error "estimation with meta-regression. Forthcoming in Psychological Methods."
    }


	
	
end
