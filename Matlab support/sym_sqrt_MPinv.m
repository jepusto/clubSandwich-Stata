function B = sym_sqrt_MPinv(A)

    [E,L] = eig(A);
    MM = L.^(-1/2);
    MM(L<=0)=0;
    B= real(E*MM*E');