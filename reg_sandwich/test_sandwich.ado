*! version 3.0 updated 15-Aug-2016
// Update by Marcelo Tyszler (tyszler.jobs@gmail.com):
// - A few typos in the comments
// - small sample F-test
//
// This procedure estimates the robust standard errors for meta regressions by way of various weighting schemes
// Original program by Eric Hedberg (eric.hedberg@me.com) based on
// Hedges, Larry V., Elizabeth Tipton, and Matthew C. Johnson. 2010. Robust variance estimation
//     in meta-regression with dependent effect size estimates. Research Synthesis Methods.
//     (www.interscience.wiley.com) DOI: 10.1002/jrsm.5
//and
//Tipton, E. (in press) Small sample adjustments for robust variance estimation with meta-regression. Forthcoming in Psychological Methods.
//
// Added by Marcelo Tyszler (tyszler.jobs@gmail.com)
// Implements the small sample f-test based on:
// Tipton and Pustejovsky (2015). Small-sample adjustments for tests of moderators and model fit using robust variance estimation in meta-regression
// 
// The implemented test is the AHZ test: Hotelling's T2 approximation, using Zhang's approach to estimate eta. The key equations are (10) and (14).
// The test is computed under the modelled approach, where W = inv(Theta) (or Theta = inv(W)) and therefore Var(b) = M = inv(X'*W*X)


capture program drop test_sandwich
program define test_sandwich, eclass byable(recall) sortpreserve
	*version 13 // Marcelo Tyszler
	version 14.1 // Marcelo Tyszler
	
	* verify matqrt is installed:
	capture which matsqrt
	if _rc!=0 {
		display as text "Package {it:matsqrt} is missing. Trying to install..."
		capture net install matsqrt.pkg, from (http://www.stata.com/users/jpitblado)
		if _rc!=0 {
			display as error "You need the package {it:matsqrt}, but an error occurred while trying to install it:"
			error _rc
		}
		display as text "Install successful" _newline
	}



    set type double
    syntax [varlist(default=none)], [cons]
	
	tempname  C_Ftest ///
			  js jt Var_d_st temp_calc Sum_Var_d_st ///
			 Omega_Ftest matrix_Ftest middle_Omega ///
			 Big_P_relevant Big_PThetaP_relevant Pi_Theta_Pi Pi_relevant Pj_relevant Big_PP middle_PThetaP ///
			 Fconstant ///
			 clusternumber ///
			 cluster_list ///
			 cluster ///
			 eta_Ftest ///
			 Q_Ftest z_Ftest D_Ftest ///
			 b V ///
			 F_stat F_df1 F_df2 F_pvalue ///
			 MXWTWXM

	*verify that this is run after robumeta:

	if e(cmd) !="reg_sandwich" {
		display as error "{it:test_sandwich} can only be used after {it:reg_sandwich}"
		error 301
		exit
	}
	
	
	* F-test:
	* 
	* The Q statistic needs the definition of C matrix and a c vector
	*
	* p is the number of coefficients (including the constant)
	* q is the number of coefficients to be tested (1 <= q <= p)
	*
	* C is the contrast matrix (q x p) such that:
	*  H0: Cb = c
	*  For example for the test beta_s = 0, would have q=1,  c = 0 and C = [ 0 .. 1 .. 0], where entry s (between 1 and p) is 1 and 0 otherwise
	* 
	*  For the F-test, C will be a (q x p) matrix where entry rs (1 <= r <= q and 1 <= s <=p) will be 1 and 0 otherwise
	*  and c will be a (q x 1) vector of 0s
	*
	*  For example if p = 3, an F test of beta_1 = beta_3 = 0 would have c = [0 0]' and C = [1 0 0; 0 0 1]
		
	
	*count coefficients to be tested, 
	* Simultaneously check if they belong to the list of original coefficients
	*load original variables:
	local x = e(indepvars) 
	local type_VCR = e(type_VCR)
	local constant_used = e(constant_used)
	
	local q_Ftest = 0
	capture confirm existence `varlist' 
	if _rc != 6 {
		foreach current_x in `varlist' {
			* verify coefficient is in the original list
			if strpos("`x'","`current_x'")==0{
				display as error "F-test error: {it:`current_x'} does not belong to the list of coefficients from the {it:robumeta} estimation"
				error 101
				exit
			}
			
			local ++q_Ftest
		}
	}
	* verify constant term:
	* remove leading and trailing whitespace
	if "`cons'" != "" {
		if "`constant_used'" == "0" {
			display as error "Constant was not included in the estimation. It cannot be included in the tests."
			error 101
			exit
		}
		* increment q:
		local ++q_Ftest
	}
	
	* Define C
	local p = 0
    foreach v in `x' {
        local ++p
    }
	
	if "`constant_used'" == "1" {
		* increment q:
		local ++p
	}
	
	matrix `C_Ftest' = J(`q_Ftest', `p', 0) // C is initialized as a q x p matrix of zeros
	
	if "`constant_used'" == "1" {
		matrix colnames `C_Ftest' = `x' _cons
	}
	else {
		matrix colnames `C_Ftest' = `x' 	
	}
	
	local current_row = 1
	capture confirm existence `varlist' 
	if _rc != 6 {
		* for each var listed in ftest, check which column it corresponds to
		foreach current_q in `varlist'{
			matrix `C_Ftest'[`current_row', colnumb(`C_Ftest',"`current_q'")] = 1
			local ++current_row
		}
	}
	
	* If option constant is active, last column needs to be active:
	if "`cons'" != "" {
		matrix `C_Ftest'[`current_row', `p'] = 1
		local ++current_row
	}
	mata : st_matrix("`C_Ftest'", sort(st_matrix("`C_Ftest'"), -1..-`p'))
	
    * F-test:
	* We need to get sum (s = 1 to q) sum (t = 1 to q) Var(d_st) (denominator of eq(14))
	*
	* Var(d_st) can be computed using equation (8), assuming Theta = inv(W):
	* Cov(u1'*D*u2,u3'*D*u4) = sum (i=1 to m) sum(j = 1 to m) u1'*Bi*inv(W)*Bj'u4*u2'*Bi*inv(W)*Bj'*u3 + u1'*Bi*inv(W)*Bj'u3*u2'*Bi*inv(W)*Bj'*u4
	* * where m = is the numberf of studies
	*
	* we additionally set u1 = u3 = js and u2 = u4 = jt, where vector js has entry s = 1 and 0 otherwise and jt is similarly defined
	*
	* Therefore
	* Cov(js'*D*jt,js'*D*jt) =  sum (i=1 to m) sum(j = 1 to m) js'*Bi*Theta*Bj'*jt*jt'*Bi*Theta*Bj'*js + js'*Bi*Theta*Bj'*js*jt'*Bi*Theta*Bj'*jt
	* Var(d_st) = Var(js'*D*jt) = Cov(js'*D*jt,js'*D*jt) 
	*
	*
	
	local m = e(N_clusters)
	
	matrix `Big_P_relevant' = e(P_relevant)
	matrix `Big_PThetaP_relevant' = e(PThetaP_relevant)
	
	if "`type_VCR'" == "WLSp" {

		if "`constant_used'" == "1" {
			quietly : gen double `Fconstant' = 1 if e(sample)
		}
			
		
		matrix `Big_PP' = e(PP)
		
		local cluster = e(cluster)
		quietly : gen double `clusternumber' = `cluster' if e(sample)
		quietly sort `clusternumber'
		qui: tab `clusternumber' if e(sample), matrow(`cluster_list')

		
		local starti = 0
		
		forvalues i = 1/`m'{
 

						
			tempname X`i' PP`i' P`i'_relevant
			if "`constant_used'" == "1" {
				mkmat `x' `Fconstant' if e(sample) & `clusternumber' == `cluster_list'[`i',1], matrix(`X`i'')		
			}
			else {
				mkmat `x' if e(sample) & `clusternumber' == `cluster_list'[`i',1], matrix(`X`i'')	
			}
			
			local endi = `starti'+1
			local starti =  `endi' + rowsof(`X`i'')-1
			
			matrix `PP`i'' = `Big_PP'[1.. rowsof(`X`i''),1..`p']
			matrix `P`i'_relevant' = `Big_P_relevant'[1.. rowsof(`X`i''),1..`p']'
		}
		
	} 
	else {
		forvalues i = 1/`m'{
			local starti = (`i'-1)*`p'+1
			local endi = `starti'+`p'-1
						
			tempname P`i'_relevant
			
			matrix `P`i'_relevant' = `Big_P_relevant'[`starti'..`endi',1..`p']'
		}
	}
	
	matrix `MXWTWXM' = e(MXWTWXM)
	matrix `Omega_Ftest' = `C_Ftest'*`MXWTWXM'*`C_Ftest''
	matsqrt `Omega_Ftest'
	matrix `matrix_Ftest' = invsym(sq_`Omega_Ftest')*`C_Ftest'
	
	local Sum_Var_d_st = 0
	
	forvalues s = 1/`q_Ftest'{
		
		matrix `js' = J(`q_Ftest',1,0)
		matrix `js'[`s',1]=1
		
		* We use the symmetry here, since that Var_d_st = Var_d_ts
		forvalues t = `s'/`q_Ftest'{
			
			matrix `jt' = J(`q_Ftest',1,0)
			matrix `jt'[`t',1]=1
			
			* get Var(d_st)
			local Var_d_st = 0
			forvalues i = 1/`m'{
				* We use the symmetry here, since that temp(i,j) =temp(j,i)
				forvalues j = `i'/`m'{
					
					if `i' == `j'{
						local start = (`i'-1)*`p'+1
						local end = `start'+`p'-1
						
						matrix `Pi_Theta_Pi' = `Big_PThetaP_relevant'[`start'..`end',1..`p']
						
						matrix `temp_calc' = `js''*`matrix_Ftest'*`Pi_Theta_Pi'*`matrix_Ftest''*`jt'*`jt''*`matrix_Ftest'*`Pi_Theta_Pi'*`matrix_Ftest''*`js'+ ///
											 `js''*`matrix_Ftest'*`Pi_Theta_Pi'*`matrix_Ftest''*`js'*`jt''*`matrix_Ftest'*`Pi_Theta_Pi'*`matrix_Ftest''*`jt'

					}
					else {
						
						if "`type_VCR'" == "WLSp" {
						
							* for this weighitng type, we need to use V for Theta
							* inv(sq_Omega)*C*M*Xi'*Wi*Ai*(I-X*M*X'*W)i*V*(I-X*M*X'*W)j'*Aj'*Wj*Xj*M'*C*inv(sq_Omega)
							* we use the fact that 
							* (I-X*M*X'*W)i*V*(I-X*M*X'*W)j' = 
							*
							* if i!=j
							* - Vi*Wi*Xi*M*Xj'   - Xi*M*Xj'*Wj*Vj     + Xi*(M*X'*W*V*W*X*M)*Xj'
							* we call PPj = Vj*Wj*Xj*M
											
							
							matrix `middle_PThetaP' = -`PP`i''*`X`j'''-`X`i''*`PP`j''' + `X`i''*`MXWTWXM'*`X`j'''

							
						}
						else {
							
							matrix `middle_PThetaP' = `MXWTWXM'
							
						}
						
						
						matrix `temp_calc' = `js''*`matrix_Ftest'*`P`i'_relevant'*`middle_PThetaP'*`P`j'_relevant''*`matrix_Ftest''*`jt'*`jt''*`matrix_Ftest'*`P`i'_relevant'*`middle_PThetaP'*`P`j'_relevant''*`matrix_Ftest''*`js' + ///
											 `js''*`matrix_Ftest'*`P`i'_relevant'*`middle_PThetaP'*`P`j'_relevant''*`matrix_Ftest''*`js'*`jt''*`matrix_Ftest'*`P`i'_relevant'*`middle_PThetaP'*`P`j'_relevant''*`matrix_Ftest''*`jt'
					
						matrix drop `middle_PThetaP'
					}
					* We use the symmetry here, since that temp(i,j) =temp(j,i)
					if `i'==`j'{
						local Var_d_st = `Var_d_st' + `temp_calc'[1,1]	
					}
					else {
						local Var_d_st = `Var_d_st' + 2*`temp_calc'[1,1]
					}
					
					
				} 
			}
			
			* update SumSum Var(d_st)
			* We use the symmetry here, since that Var_d_st = Var_d_ts
			if `s'==`t' {
				local Sum_Var_d_st = `Sum_Var_d_st' + `Var_d_st'
			}
			else {
				local Sum_Var_d_st = `Sum_Var_d_st' + 2*`Var_d_st'
			}
			
		}
	}
			
		
			
	* eta needs to be computed according to equation (14):
	* eta = q*(q+1) / [sum(s=1 to q) sum(t=1 to q) Var(d_st)]
	
	local eta_Ftest = (`q_Ftest'*(`q_Ftest'+1))/`Sum_Var_d_st'
	
	* z = Omega^(-1/2)(Cb-c)
	* D = Omega^(-1/2)*C*VR*C'*Omega^(-1/2)
	* Q = z'inv(D)z (equation 6)
	
	matrix `b' = e(b)
	matrix `V' = e(V)
	
	matrix `z_Ftest' = inv(sq_`Omega_Ftest')*(`C_Ftest'*`b'')
	matrix `D_Ftest' = inv(sq_`Omega_Ftest')*`C_Ftest'*`V'*`C_Ftest''*inv(sq_`Omega_Ftest')
	matrix `Q_Ftest' = `z_Ftest''*inv(`D_Ftest')*`z_Ftest'
	
	* Now we can compute the F-statistic:
	* (eta - q + 1)/(eta*q) * Q  follows F(q, eta - q + 1) distribution
	local F_stat = ((`eta_Ftest' - `q_Ftest' + 1)/(`eta_Ftest'*`q_Ftest'))* `Q_Ftest'[1,1]
	local F_df1 = `q_Ftest'
	local F_df2 = `eta_Ftest' - `q_Ftest' + 1
	
	local F_pvalue = Ftail(`F_df1',`F_df2',`F_stat')
	
	
	* F-test:
	* display some results
	display in b _newline
	display in b  "Small Sample Corrected F-test:" 
	display _col(10) in b  "F(" as result %5.4f `F_df1' "," as result %5.4f `F_df2' ")" _col(30) "=" _col(35) as result  %5.4f `F_stat'
	display _col(10)  "Prob > F" _col(30) "=" _col(35) as result  %5.4f `F_pvalue'
	
	/* Naive-F */
	/*
	local m_min_p = `m'-`p'
	local Naive_F = `Q_Ftest'[1,1]/`q_Ftest'
	local Naive_F_p = Ftail(`q_Ftest',`m_min_p',`Naive_F')
	disp "Naive F = `Naive_F', q: `q_Ftest', m: `m', p: `p', (m-p) = `m_min_p', p-value: `Naive_F_p'"
	*/
	
	
	* F-test:
	* save some results
	ereturn scalar F_stat = `F_stat'
	ereturn scalar F_df1 = `F_df1'
	ereturn scalar F_df2 = `F_df2'
	ereturn scalar F_pvalue = `F_pvalue'
	ereturn scalar F_eta = `eta_Ftest'	
	
	
** Clean:
	matrix drop sq_`Omega_Ftest'

end
