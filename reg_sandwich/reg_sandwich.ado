// github repository: https://github.com/jepusto/clubSandwich-Stata
//
*! version 0.0 updated 02-March-2017
// Updated by Marcelo Tyszler (tyszler.jobs@gmail.com):
//
// Wrapper function for reg_sandwich:
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
	
	*timer clear 
	*timer on 5
	
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
    tempvar cons w wh wfinal    ///
     prime_resid  v_n  clusternumber  ///
	theta
	   
    *specifiy the temporary scalars and matrixes
    tempname ///
		X V b W ///
		min_n max_n ///
		M ///
		MXWTWXM ///
		XWAeeAWX Big_P_relevant Big_PThetaP_relevant Big_PP ///
		PThetaP ///
		Aj Wj Xj ej Bj ///
		Tj ///
		_dfs  ///
		cluster_list ///
		C_ttest ///
		Omega_ttest matrix_ttest ///
		temp_calc ///
		prob  ///
		T Dj inv_Bj ///
		Xi Wi Ti Bi Ai inv_Bi


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
        local x = regexr("`x' ","o\.`v' ","")
    }

    if "`x'" == "." local x ""
	
	
	** for absorb:
	if "`main_function'" == "areg" {
		* predict:
		local new_x = ""
		foreach xr of varlist `x' {
			tempvar _RS`xr'
			local new_x = trim("`new_x'") + " " + "`_RS`xr''"
			noisily capture: reg `xr' i.`absorb' `weight_call' if `touse'
			qui: predict `_RS`xr'' if `touse' , residuals 
		}
	}
	** call main regression:
	*disp "timer 1 start: main regression"
	*timer on 1
	*disp "`main_function' `t' `x' `weight_call' if `touse', `constant' cluster(`cluster') `absorb_call'"
	noisily capture: `main_function' `t' `x'  `weight_call' if `touse', `constant' cluster(`cluster') `absorb_call'
	*timer off 1
	*disp "timer 1 off"
	if "`main_function'" == "areg" {
		local old_x = "`x'"
		local x = "`new_x'"
	}
	
	
	** prep for small sample reduced t-test:
	* (based on robumetaMt12.ado:
		
	matrix p = rowsof(e(V))
	local p = p[1,1]
	matrix drop p
	
	if "`main_function'" == "areg" {
		*ignore constant
		local --p
	}
	
	*capture ids
    quietly : gen double `clusternumber' = `cluster' if `touse'
    quietly sort `clusternumber' `_sortindex'
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
		qui: gen double `theta' = 1 if `touse'
		local type_VCR = "OLS"
	} 
	else {
		* weights = WLS
		local model_weights = substr("`exp'",2,.)
		quietly : gen double `wfinal' = `model_weights' if `touse' 
		
		if "`weight'"=="aweight"{				
			qui: gen double `theta' = 1/`model_weights'  if `touse'
			local type_VCR = "WLSa"
		}
		else{
		* p-weights
		* Working variance is I
			qui: gen double `theta' = 1 if `touse'
			local type_VCR = "WLSp"
		}
		
	
	}
	
	* count
	*******
	quietly : by `clusternumber': egen double `v_n' = count(`theta') if `touse'
				quietly sum `v_n' if `touse'				
				scalar `min_n' = r(min) 
				scalar `max_n' = r(max)
	

	*disp "timer 2 start: MXWTWXM"
	*timer on 2
	if "`constant'"=="" & "`main_function'" != "areg" {
		mkmat `x' `cons' if `touse', matrix(`X')
		matrix colnames `X' = `x' _cons
	}
	else{
		mkmat `x'  if `touse', matrix(`X')
		matrix colnames `X' = `x'
	}
	
	if "`type_VCR'" =="OLS"{
		*matrix `M' = invsym(`X'' * `X')
		mata: st_matrix("`M'",invsym(st_matrix("`X'")' * st_matrix("`X'")))
		matrix `MXWTWXM' =  `M'
	} 
	else {
	
		mkmat `wfinal' if `touse', matrix(`W')
		matrix `W' = diag(`W')
	
		*matrix `M' = invsym(`X'' * `W' * `X')
		mata: st_matrix("`M'",invsym(st_matrix("`X'")' * st_matrix("`W'") *st_matrix("`X'")))
		
	}
				
	
	
	if "`type_VCR'" == "WLSp" {	
		*matrix `MXWTWXM' =  `M'*`X''*`W'*`W'*`X'*`M'
		mata: st_matrix("`MXWTWXM'", st_matrix("`M'")*st_matrix("`X'")'*st_matrix("`W'")*st_matrix("`W'")*st_matrix("`X'")*st_matrix("`M'"))
	}
	else if "`type_VCR'" == "WLSa" {
		
		mkmat `theta'   if `touse', matrix(`T')
		matrix `T' = diag(`T')
		matrix `MXWTWXM' =  `M'
		matrix drop `T'
	}
	
	if "`type_VCR'" ~="OLS"{
		matrix drop `W' 
	}
	
	matrix drop `X'
	*timer off 2
	*disp "timer 2 off"
	
	/********************************************************************/
    /*    Variance covariance matrix estimation for standard errors     */
	/* 																    */
	/*    And F-test												    */
    /********************************************************************/
	qui: predict `prime_resid', residuals
	
	matrix `XWAeeAWX' = J(`p', `p', 0)
		
	local current_jcountFtest = 0
	local first_cluster = 1
	*local endj = 0
	tempname Pj_relevant  Pj_Theta_Pj_relevant
	tempname evecs evals sq_inv_Bj
	tempname PPj
	*disp "timer 3 start: Ajs and save per cluster"
	
    foreach j in `idlist' {
	*timer on 3
	*disp "cluster `j'"
		
		if "`constant'"=="" & "`main_function'" != "areg" {
			mkmat `x' `cons' if `touse' & `clusternumber' == `j', matrix(`Xj')
			matrix colnames `Xj' = `x' _cons
		} 
		else {
			mkmat `x'  if `touse' & `clusternumber' == `j', matrix(`Xj')
			matrix colnames `Xj' = `x'
		}
		
		mkmat `theta' if `touse' & `clusternumber' == `j', matrix(`Tj')
		matrix `Tj' = diag(`Tj')  
		
		****Adjustment matrix
		* 
		* we use that Bj = 
		* Dj*[(I-X*M*X'*W)j*T*(I-X*M*X'*W)j']Dj = 
		*
		* Dj*[Tj - Tj*(Wj*Xj*M*Xj') - (Xj*M*Xj'*Wj)*Tj + Xj*(M*X'*W*V*W*X*M)*Xj']*Dj
		*
		* For OLS this simplifies to (Dj = I):
		* Tj - Xj*M*Xj'
		*
		* For WLSp, this simplifies to (Dj = I):
		* Tj - Wj*Xj*M*Xj' - Xj*M*Xj'Wj + Xj'MXWWXM*Xj'
		*
		* For WLSa, this simplified to (Wj*Tj = Ij):
		* Dj*[Tj-Xj*M*Xj]*Dj
				
		if "`type_VCR'" == "OLS" {
			*matrix `Bj'=`Tj'-`Xj'*`M'*`Xj''
			mata: st_matrix("`Bj'",st_matrix("`Tj'")-st_matrix("`Xj'")*st_matrix("`M'")*st_matrix("`Xj'")')
		}
		else if "`type_VCR'" == "WLSp" {
			mkmat `wfinal' if `touse' & `clusternumber' == `j', matrix(`Wj')
			matrix `Wj' = diag(`Wj')  
			*matrix `Bj'=`Tj'-`Wj'*`Xj'*`M'*`Xj''-`Xj'*`M'*`Xj''*`Wj'+ `Xj'*`MXWTWXM'*`Xj''
			mata: st_matrix("`Bj'", st_matrix("`Tj'")-st_matrix("`Wj'")*st_matrix("`Xj'")*st_matrix("`M'")*st_matrix("`Xj'")'-st_matrix("`Xj'")*st_matrix("`M'")*st_matrix("`Xj'")'*st_matrix("`Wj'")+ st_matrix("`Xj'")*st_matrix("`MXWTWXM'")*st_matrix("`Xj'")')
		}
		else if "`type_VCR'" == "WLSa" {
			mkmat `wfinal' if `touse' & `clusternumber' == `j', matrix(`Wj')
			matrix `Wj' = diag(`Wj')
			matrix `Dj' = cholesky(`Tj')
			*matrix `Bj' = `Dj''*(`Tj'-`Xj'*`M'*`Xj'')*`Dj'
			mata: st_matrix("`Bj'", st_matrix("`Dj'")'*(st_matrix("`Tj'")-st_matrix("`Xj'")*st_matrix("`M'")*st_matrix("`Xj'")')*st_matrix("`Dj'"))
		}
		
		* Symmetric square root of the Moore-Penrose inverse of Bj
		
		mat symeigen `evecs' `evals' = `Bj'
		mata: st_matrix( "`sq_inv_Bj'", st_matrix( "`evecs'")*diag(editmissing(st_matrix( "`evals'"):^(-1/2),0))*st_matrix( "`evecs'")')
											
		if "`type_VCR'" == "WLSa" {
			*matrix `Aj' = `Dj'*(`sq_inv_Bj')*`Dj''
			mata: st_matrix("`Aj'" , st_matrix("`Dj'")*(st_matrix("`sq_inv_Bj'"))*st_matrix("`Dj'")')
		}
		else {
			matrix `Aj' = (`sq_inv_Bj')
		}
		matrix drop `sq_inv_Bj'	
		
        mkmat `prime_resid' if `touse' & `clusternumber' == `j', matrix(`ej')
		
		if "`type_VCR'" == "OLS" {
			*matrix `XWAeeAWX' = (`Xj'' * `Aj' * `ej' * `ej'' * `Aj' * `Xj') + `XWAeeAWX'
			mata: st_matrix("`XWAeeAWX'",   (st_matrix("`Xj'")'* st_matrix("`Aj'") * st_matrix("`ej'") * st_matrix("`ej'")' * st_matrix("`Aj'") * st_matrix("`Xj'")) + st_matrix("`XWAeeAWX'"))
		}
		else {
		
			*matrix `XWAeeAWX' = (`Xj'' * `Wj' * `Aj' * `ej' * `ej'' * `Aj' * `Wj' * `Xj') + `XWAeeAWX'
			mata: st_matrix("`XWAeeAWX'",   (st_matrix("`Xj'")'* st_matrix("`Wj'")* st_matrix("`Aj'") * st_matrix("`ej'") * st_matrix("`ej'")' * st_matrix("`Aj'") * st_matrix("`Wj'")* st_matrix("`Xj'")) + st_matrix("`XWAeeAWX'"))
		}

		
		**** F-test:
		* 
		* To compute the degress of freedom we need P:
		* Psi = (I-Hx)i'*Ai*Wi*Xi*M*C*gs
		*
		* These matrices are needed to compute the terms Psi'*Theta*Ptj:
		*  gs'*C'*M*Xi'*Wi*Ai*(I-Hx)i*Theta*(I-Hx)j'*Aj*Wj*Xj*M*C*gt
		*
		*
		* We save just the "middle" portion, which is independent of C and gs:
		* 
		* We use the fact that Hx = X*M*X'W and
		* (I-X*M*X'*W)i*T*(I-X*M*X'*W)j' = 
		*
		* if i==j
		* Tj - Tj*(Wj*Xj*M*Xj') - (Xj*M*Xj'*W)*Tj + Xj*(M*X'*W*V*W*X*M)*Xj' 
		*
		* For OLS this simplifies to:
		* Tj - Xj*M*Xj'
		*
		* For WLSp, this simplifies to (Dj = I):
		* Tj - Wj*Xj*M*Xj' - Xj*M*Xj'Wj + Xj'MXWWXM*Xj'
		*
		* For WLSa, this simplified to:
		* Tj - Xj*M*Xj
		* 
		*
		* and we call M*Xi'*Wi*Ai*(I-Hx)i*Theta*(I-Hx)j'*Aj*Wj*Xj*M:
		* Pi_Theta_Pi_relevant
		*
		*
		* if i!=j
		* - Ti*Wi*Xi*M*Xj'   - Xi*M*Xj'*Wj*Tj     + Xi*(M*X'*W*T*W*X*M)*Xj'
		*
		* For OLS this simplifies to:
		* - Xi*M*Xj'
		*
		* For WLSp, this simplifies to:
		* - Wi*Xi*M*Xj'   - Xi*M*Xj'*Wj     + Xi*(M*X'*W*W*X*M)*Xj'
		*
		* For WLSa, this simplified to:
		* - Xi*M*Xj' 
		*  
		* For OLS and WLSa we call M*Xi'*Wi*Ai*Xi:
		* Pi_relevant (and ignore the (min) sign, since it will be cancelled out after multiplication)
		*
		* For WLSp we call  M*Xi'*Wi*Ai
		* Pi_Pj_relevant, (this is more efficient to save)
		* 
		* and additionally save M*Xi'*Wi*Ai as PPi
		
		local current_jcountFtest = `current_jcountFtest'+1
        
		
		
		if "`type_VCR'" == "OLS" {
			/*
			matrix `Pj_Theta_Pj_relevant' = ///
													`M'*`Xj''*`Aj'* ///
													(`Bj')* ///
													`Aj''*`Xj'*`M'' // p x p
			*/
			mata: st_matrix("`Pj_Theta_Pj_relevant'", st_matrix("`M'")*st_matrix("`Xj'")'*st_matrix("`Aj'")* (st_matrix("`Bj'"))* st_matrix("`Aj'")'*st_matrix("`Xj'")*st_matrix("`M'")') 
			// p x p
														
			*matrix `Pj_relevant' =  `M'*`Xj''*`Aj'*`Xj' // p x p
			mata: st_matrix("`Pj_relevant'",  st_matrix("`M'")*st_matrix("`Xj'")'*st_matrix("`Aj'")*st_matrix("`Xj'") )
		
													

		}
		else if "`type_VCR'" == "WLSp" {	
			/*
			matrix `Pj_Theta_Pj_relevant' = ///
													`M'*`Xj''*`Wj'*`Aj'* ///
													(`Bj')* ///
													`Aj''*`Wj'*`Xj'*`M'' // p x p
			*/
			mata: st_matrix("`Pj_Theta_Pj_relevant'", st_matrix("`M'")*st_matrix("`Xj'")'* st_matrix("`Wj'")* st_matrix("`Aj'")* (st_matrix("`Bj'"))* st_matrix("`Aj'")'* st_matrix("`Wj'")* st_matrix("`Xj'")*st_matrix("`M'")') 
			
													
			*matrix `Pj_relevant' =  `M'*`Xj''*`Wj'*`Aj' // p x kj
			mata: st_matrix("`Pj_relevant'",  st_matrix("`M'")*st_matrix("`Xj'")'*st_matrix("`Wj'")*st_matrix("`Aj'") )
			matrix P`current_jcountFtest'_relevant = `Pj_relevant'
			
			
			*matrix `PPj' = `Wj'*`Xj'*`M' // kj x p
			mata: st_matrix("`PPj'", st_matrix("`Wj'")*st_matrix("`Xj'")*st_matrix("`M'"))
			matrix PP`current_jcountFtest' = `PPj'
		}
		else if "`type_VCR'" == "WLSa" {
			/*
			matrix `Pj_Theta_Pj_relevant' = ///
													`M'*`Xj''*`Wj'*`Aj'* ///
													(`Tj'-`Xj'*`M'*`Xj'')* ///
													`Aj''*`Wj'*`Xj'*`M'' //p x p
			*/
			mata: st_matrix("`Pj_Theta_Pj_relevant'", st_matrix("`M'")*st_matrix("`Xj'")'* st_matrix("`Wj'")* st_matrix("`Aj'")* (st_matrix("`Tj'")-st_matrix("`Xj'")*st_matrix("`M'")*st_matrix("`Xj'")')* st_matrix("`Aj'")'* st_matrix("`Wj'")* st_matrix("`Xj'")*st_matrix("`M'")') 
			
													
			*matrix `Pj_relevant' =  `M'*`Xj''*`Wj'*`Aj'*`Xj' // p x p
			mata: st_matrix("`Pj_relevant'",  st_matrix("`M'")*st_matrix("`Xj'")'*st_matrix("`Wj'")*st_matrix("`Aj'")*st_matrix("`Xj'") )
		}
	
	
		
		* save for later
		if `first_cluster'==1 {
			matrix `Big_PThetaP_relevant' = `Pj_Theta_Pj_relevant'
			matrix `Big_P_relevant' = `Pj_relevant''
			if "`type_VCR'" == "WLSp" {
				matrix `Big_PP' = `PPj'
			}
			local first_cluster = 0
		}
		else {	
			matrix `Big_PThetaP_relevant' = [`Big_PThetaP_relevant' \ `Pj_Theta_Pj_relevant']
			matrix `Big_P_relevant' = [`Big_P_relevant' \ `Pj_relevant'']
			
			if "`type_VCR'" == "WLSp" {
				matrix `Big_PP' = [`Big_PP' \ `PPj']
			
			}
		}
		
		
			
		*timer off 3
	
    }
	*disp "timer 3 off"

	
	*matrix drop `Aj' `Wj' `Xj' `middle_Aj'  `ej'
		
	* RVE estimator
	matrix `V' = `M' * `XWAeeAWX' * `M'
	
	
	* T-test, using as a special case of an F-test:
	if "`type_VCR'" == "WLSp" {
		
			qui: tab `clusternumber' if `touse', matrow(`cluster_list')
			forvalues i = 1/`m'{
		
				*tempname X`i'
				
				if "`constant'"=="" & "`main_function'" != "areg" {
					mkmat `x' `cons' if `touse' & `clusternumber' == `cluster_list'[`i',1], matrix(X`i')
				} 
				else {
					mkmat `x'  if `touse' & `clusternumber' == `cluster_list'[`i',1], matrix(X`i')
				}
			}
	}
	
	*matrix `_dfs' =  J(1,`p', 0) 
	*disp "timer 4 start: T-tests i:j"
	*timer on 4
	mata: st_matrix("`_dfs'", reg_sandwich_ttests("`type_VCR'", `m', `p', st_matrix("`Big_PThetaP_relevant'"),  st_matrix("`Big_P_relevant'"), st_matrix("`M'"),  st_matrix("`MXWTWXM'")))
	if "`type_VCR'" == "WLSp" {
	
			forvalues i = 1/`m'{
		
				matrix drop X`i'
				matrix drop PP`i'
				matrix drop P`i'_relevant
				
			}
	}

	
	/*
	forvalues i = 1/`m'{
		* We use the symmetry here, since that temp(i,j) =temp(j,i)
		forvalues j = `i'/`m'{
		
		*disp "`i':`j'"
			if `i' == `j'{
	
				matrix  `PThetaP' = `P`i'_Theta_P`i'_relevant'
	
			}
			else {
				* i ~= j 
				* M*Xi'*Wi*Ai*(I-Hx)i*Theta*(I-Hx)j'*Aj*Wj*Xj*M
				* and we call M*Xi'*Wi*Ai*(I-Hx)i*Theta*(I-Hx)j'*Aj*Wj*Xj*M:
				* Pi_Theta_Pi_relevant
				*
				*
				* if i!=j
				* (I-Hx)i*Theta*(I-Hx)j'
				*
				* For OLS this simplifies to:
				* - Xi*M*Xj'
				*
				* For WLSp, this simplifies to:
				* - Wi*Xi*M*Xj'   - Xi*M*Xj'*Wj     + Xi*(M*X'*W*W*X*M)*Xj'
				*
				* For WLSa, this simplified to:
				* - Xi*M*Xj' 
				* 
				* For OLS and WLSa we call M*Xi'*Wi*Ai*Xi:
				* Pi_relevant (and ignore the (min) sign, since it will be cancelled out after multiplication)
				*
				* For WLSp we call  M*Xi'*Wi*Ai
				* Pi_relevant,
				* and M*Xi'*Wi*Ai as PPi
				
				if "`type_VCR'" == "OLS" {
					matrix  `PThetaP' = `P`i'_relevant'*`M'*`P`j'_relevant''														
				}
				else if "`type_VCR'" == "WLSp" {
				
					matrix  `PThetaP' = `P`i'_relevant'* ///
										(-`PP`i''*`X`j'''-`X`i''*`PP`j''' + `X`i''*`MXWTWXM'*`X`j''')* ///
										`P`j'_relevant''
					
					
				}
				else if "`type_VCR'" == "WLSa" {
					matrix  `PThetaP' = `P`i'_relevant'*`M'*`P`j'_relevant''
				}
			}
		
			forvalues coefficient = 1/`p' {
				
				* prep C
					
				* * Define C
				matrix `C_ttest' = J(1, `p', 0) // C is initialized as a 1 x p matrix of zeros
				
				matrix `C_ttest'[1,`coefficient'] = 1
				
				
				matrix `Omega_ttest' = `C_ttest'*`MXWTWXM'*`C_ttest''	
				
				local temp_val = `Omega_ttest'[1,1]
				matrix `matrix_ttest' = (1/sqrt(`temp_val'))*`C_ttest'
				matrix `temp_calc' = `matrix_ttest'*`PThetaP'*`matrix_ttest''
				local  temp_calc2 = 2*((`temp_calc'[1,1])^2)
				
				* We use the symmetry here, since that temp(i,j) =temp(j,i)
				if `i'==`j'{
					matrix `_dfs'[1,`coefficient'] = `_dfs'[1,`coefficient'] + `temp_calc2'
					
				}
				else {
					matrix `_dfs'[1,`coefficient'] = `_dfs'[1,`coefficient'] + 2*`temp_calc2'
				}
				
			}
			matrix drop `PThetaP'
		
		} 
		
		
	}
	*/
	*timer off 4	
	*disp "timer 4 off"
	forvalues coefficient = 1/`p' {
		matrix `_dfs'[1,`coefficient'] = 2/`_dfs'[1,`coefficient']
	}
	
     
    display _newline
    display as text "Robust standard error estimation using " as result "`main_call_display'`absorb_call_display'"



	if "`main_function'" == "areg" {
		if "`type_VCR'" == "WLSp" {
			mkmat `x'  if `touse', matrix(`X')
			matrix colnames `X' = `old_x'
		}
		local x = "`old_x'"
			
	}
	
	*name the rows and columns of the matrixes
	if "`constant'"=="" & "`main_function'" != "areg" {	
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
	tempname b_temp
	matrix `b_temp' = e(b)
	matrix `b' = `b_temp'[1, 1..`p']
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
	
	if "`constant'"=="" & "`main_function'" != "areg" {
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

   /*********************/
    /*  post results  */
    /*********************/
	ereturn post `b' `V', obs(`nobs') depname(`t') esample(`touse')
	
	ereturn local type_VCR "`type_VCR'"
	ereturn local vce "cluster"
	ereturn local vcetype "Cluster Robust Small Sample Corrected"
	ereturn scalar N_clusters = `m'
	ereturn matrix dfs = `_dfs'
	
	ereturn local clustvar = "`cluster'"
	
	if "`type_VCR'" ~= "OLS" {
		ereturn local wtype = "`weight'"
		ereturn local wexp = "`exp'"
	}
	
	
	
	if "`main_function'" == "areg" {
			ereturn local absvar = "`absorb'"
	}
    
	* Ftest	
	ereturn matrix P_relevant = `Big_P_relevant' 
	ereturn matrix PThetaP_relevant = `Big_PThetaP_relevant'
	if "`type_VCR'" == "WLSp" {
		ereturn matrix PP = `Big_PP'
		if "`main_function'" == "areg" {
			ereturn matrix Ur = `X'
		}

	}
	
	ereturn matrix MXWTWXM = `MXWTWXM'			
	ereturn local indepvars `x'
	
	if "`constant'"=="" & "`main_function'" != "areg" {
		ereturn local constant_used = 1
	}
	else {
		ereturn local constant_used = 0
	}
	
	ereturn local cmdline "reg_sandwich `0'"
	ereturn local cmd "reg_sandwich"
	*timer off 5
	*disp "timers:"
	*timer list
end
	
