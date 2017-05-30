clear all
mata: mata clear

version 14
mata

real matrix function matsqrt(real matrix A)
{
	real matrix sq, X
	real colvector L
	symeigensystem(A, X, L)
	sq = X * diag(sqrt(L)) * X'
	return(sq)
}

end

sysuse auto
qui reg price mpg weight length turn displacement
mat v = e(V)

mata: vsq = matsqrt(st_matrix("v")); st_matrix("vsq",vsq)

// check against Jeff Pitblado's routine
matsqrt v
matlist sq_v
matlist vsq
mat diff = sq_v - vsq
matlist diff

