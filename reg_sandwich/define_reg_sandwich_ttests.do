// Created by Marcelo Tyszler (tyszler.jobs@gmail.com):
// 
// This do file simply creates and compiles the mata supporting function
// reg_sandwich_ttests
// 
// which calculates all ttests associated with reg_sandwich
//
//
// after the function defition the file reg_sandwich_ttest.mo is save
// and this file needs to be supplied together with reg_sandwich_ttests.mo


capture mata: mata drop reg_sandwich_ttests()
mata:

real vector reg_sandwich_ttests(string scalar type_VCR, real scalar m, real scalar p, matrix Big_PThetaP_relevant, matrix Big_P_relevant, matrix M, matrix MXWTWXM){
	
	_dfs =  J(1,p, 0) 
	endi = 0
	for (i=1; i<=m; i++) {
		// We use the symmetry here, since that temp(i,j) =temp(j,i)
		starti = endi+1
		endi  =  starti + p - 1
		
		// initialize startj and endj as same position as i
		startj = starti
		endj = endi
		for (j=i; j<=m; j++) {
			
			if (i == j) {
	
				PThetaP = Big_PThetaP_relevant[(starti .. endi) , (1 .. p)]
	
			}
			else {
				/* i ~= j 
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
				*/
				
				if (type_VCR == "OLS") {
					PThetaP = Big_P_relevant[(starti .. endi),(1 .. p)]'*M*Big_P_relevant[(startj .. endj),(1 .. p)]														
				}
				else if (type_VCR == "WLSp") {
					/*
					PThetaP = Big_P_relevant* ///
										(-`PP`i''*`X`j'''-`X`i''*`PP`j''' + `X`i''*`MXWTWXM'*`X`j''')* ///
										Big_P_relevant'
					*/
					
				}
				else if (type_VCR == "WLSa") {
					PThetaP = Big_P_relevant[(starti .. endi),(1 .. p)]'*M*Big_P_relevant[(startj .. endj),(1 .. p)]
				}
			}
		
			for (coefficient=1; coefficient<=p; coefficient++) {
				
				// prep C
				C_ttest = J(1, p, 0) // C is initialized as a 1 x p matrix of zeros
				
				// Define C				
				C_ttest[1,coefficient] = 1
				
				
				Omega_ttest = C_ttest*MXWTWXM*C_ttest'	
				
				temp_val = Omega_ttest[1,1]
				matrix_ttest = (1/sqrt(temp_val))*C_ttest
				temp_calc = matrix_ttest*PThetaP*matrix_ttest'
				temp_calc2 = 2*((temp_calc[1,1])^2)
				
				// We use the symmetry here, since that temp(i,j) =temp(j,i)
				if (i==j){
					_dfs[1,coefficient] = _dfs[1,coefficient] + temp_calc2
					
				}
				else {
					_dfs[1,coefficient] = _dfs[1,coefficient] + 2*temp_calc2
				}
				
			}
		
			// Update startj and endj
			startj = endj+1
			endj  =  startj + p - 1

		} 
		
		
	}
	
	return(_dfs)
}

/* Save and clear */
mata mosave reg_sandwich_ttests(), replace
mata drop reg_sandwich_ttests()

end
