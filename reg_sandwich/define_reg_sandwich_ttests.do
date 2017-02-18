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

real vector reg_sandwich_t(string scalar type_VCR, real scalar m, real scalar p, matrix Big_PThetaP_relevant, matrix Big_P_relevant, matrix M, matrix MXWTWXM){

	for (i=1; i<=m; i++) {
		// We use the symmetry here, since that temp(i,j) =temp(j,i)
		for (j=i; j<=m; j++) {
		
			if (i == j) {
	
				matrix  PThetaP = `P`i'_Theta_P`i'_relevant'
	
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
					matrix  PThetaP = `P`i'_relevant'*M*`P`j'_relevant''														
				}
				else if (type_VCR == "WLSp") {
				
					matrix  PThetaP = `P`i'_relevant'* ///
										(-`PP`i''*`X`j'''-`X`i''*`PP`j''' + `X`i''*`MXWTWXM'*`X`j''')* ///
										`P`j'_relevant''
					
					
				}
				else if (type_VCR == "WLSa") {
					matrix  PThetaP = `P`i'_relevant'*M*`P`j'_relevant''
				}
			}
		
			for (coefficient=1; coefficient<=p; coefficient++) {
				
				// prep C
				matrix C_ttest = J(1, p, 0) // C is initialized as a 1 x p matrix of zeros
				
				// Define C				
				matrix C_ttest[1,coefficient] = 1
				
				
				matrix Omega_ttest = C_ttest*MXWTWXM*C_ttest'	
				
				local temp_val = Omega_ttest[1,1]
				matrix matrix_ttest = (1/sqrt(temp_val))*C_ttest
				matrix temp_calc = matrix_ttest*PThetaP*matrix_ttest'
				local  temp_calc2 = 2*((temp_calc[1,1])^2)
				
				// We use the symmetry here, since that temp(i,j) =temp(j,i)
				if (i==j){
					matrix _dfs[1,coefficient] = _dfs[1,coefficient] + temp_calc2
					
				}
				else {
					matrix _dfs[1,coefficient] = _dfs[1,coefficient] + 2*temp_calc2
				}
				
			}
			matrix drop `PThetaP'
		
		} 
		
		
	}
}

/* Save and clear */
mata mosave reg_sandwich_ttests(), replace
mata drop reg_sandwich_ttests()

end
