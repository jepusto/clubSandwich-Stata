// Created by Marcelo Tyszler (tyszler.jobs@gmail.com):
// 
// This do file simply creates and compiles the mata supporting function
// test_sandwich_ftests
// 
// which calculates all ttests associated with test_sandwich
//
//
// after the function defition the file test_sandwich_ftest.mo is save
// and this file needs to be supplied together with test_sandwich_ftests.mo


capture mata: mata drop test_sandwich_ftests()
mata:

string scalar test_sandwich_ftests(string scalar type_VCR, real scalar q_Ftest,  real scalar m, real scalar p, Big_PThetaP_relevant, Big_P_relevant, MXWTWXM, matrix_Ftest, C_Ftest){
	
	Sum_temp_calc2 = 0
	
	for (s=1; s<=q_Ftest; s++) {
		
		gs = matrix_Ftest[1..q_Ftest,s]
		
		// We use the symmetry here, since that temp_calc2 = Var_d_ts
		for (t=s; t<=q_Ftest; t++) {
			
			gt = matrix_Ftest[1..q_Ftest,t]
			
			// get Var(d_st)
			temp_calc2 = 0
			endi = 0
			
			for (i=1; i<=m; i++) {
				starti = endi+1
				endi  =  starti + p - 1
				// We use the symmetry here, since that temp(i,j) =temp(j,i)
				startj = starti
				endj = endi
				
				// Read auxiliary variables for WLSp
				if (type_VCR == "WLSp") {
					Pi_relevant=valofexternal("P" + strofreal(i) + "_relevant")		
					PPi = valofexternal("PP" + strofreal(i))
					Xi = valofexternal("X" + strofreal(i))
				}
				for (j=i; j<=m; j++) {
				// Read auxiliary variables for WLSp
				if (type_VCR == "WLSp") {
					Pj_relevant=valofexternal("P" + strofreal(j) + "_relevant")		
					PPj = valofexternal("PP" + strofreal(j))
					Xj = valofexternal("X" + strofreal(j))
				}
					
					if (i == j){

						temp_calc = gs'*C_Ftest*Big_PThetaP_relevant[starti..endi,1..p]*C_Ftest'*gt*gt'*C_Ftest*Big_PThetaP_relevant[starti..endi,1..p]*C_Ftest'*gs+ gs'*C_Ftest*Big_PThetaP_relevant[starti..endi,1..p]*C_Ftest'*gs*gt'*C_Ftest*Big_PThetaP_relevant[starti..endi,1..p]*C_Ftest'*gt

					}
					else {
						
						if (type_VCR == "WLSp") {
							/*
							* for this weighitng type, we need to use V for Theta
							* inv(sq_Omega)*C*M*Xi'*Wi*Ai*(I-X*M*X'*W)i*V*(I-X*M*X'*W)j'*Aj'*Wj*Xj*M'*C*inv(sq_Omega)
							* we use the fact that 
							* (I-X*M*X'*W)i*V*(I-X*M*X'*W)j' = 
							*
							* if i!=j
							* - Vi*Wi*Xi*M*Xj'   - Xi*M*Xj'*Wj*Vj     + Xi*(M*X'*W*V*W*X*M)*Xj'
							* we call PPj = Vj*Wj*Xj*M
							*/				
							middle_PThetaP = -PPi*Xj'-Xi*PPj' + Xi*MXWTWXM*Xj'
							temp_calc = gs'*C_Ftest* Pi_relevant*middle_PThetaP*Pj_relevant'*C_Ftest'*gt*gt'*C_Ftest*Pi_relevant*middle_PThetaP*Pj_relevant'*C_Ftest'*gs + gs'*C_Ftest*Pi_relevant*middle_PThetaP*Pj_relevant'*C_Ftest'*gs*gt'*C_Ftest*Pi_relevant*middle_PThetaP*Pj_relevant'*C_Ftest'*gt

						}
						else {
							
							temp_calc = gs'*C_Ftest*Big_P_relevant[starti..endi,1..p]'*MXWTWXM*Big_P_relevant[startj..endj,1..p]*C_Ftest'*gt*gt'*C_Ftest*Big_P_relevant[starti..endi,1..p]'*MXWTWXM*Big_P_relevant[startj..endj,1..p]*C_Ftest'*gs + gs'*C_Ftest*Big_P_relevant[starti..endi,1..p]'*MXWTWXM*Big_P_relevant[startj..endj,1..p]*C_Ftest'*gs*gt'*C_Ftest*Big_P_relevant[starti..endi,1..p]'*MXWTWXM*Big_P_relevant[startj..endj,1..p]*C_Ftest'*gt
					
							
						}
						
						
						
					}
					// We use the symmetry here, since that temp(i,j) =temp(j,i)
					if (i==j){
						temp_calc2 = temp_calc2 + temp_calc[1,1]	
					}
					else {
						temp_calc2 = temp_calc2 + 2*temp_calc[1,1]
					}
				
					// Update startj and endj
					startj = endj+1
					endj  =  startj + p - 1
				} 
			}
			
			// update SumSum temp_calc
			// We use the symmetry here, since that temp_calc2(i,j) = temp_calc2(j,i)
			if (s==t) {
				Sum_temp_calc2 = Sum_temp_calc2 + temp_calc2
			}
			else {
				Sum_temp_calc2 = Sum_temp_calc2 + 2*temp_calc2
			}
			
		}
	}
	
	return(strofreal(Sum_temp_calc2))
			
}

/* Save and clear */
mata mosave test_sandwich_ftests(), replace
mata drop test_sandwich_ftests()

end
