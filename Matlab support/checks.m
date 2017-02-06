
%A_i = D_i'*sqrtm(pinv(B_i))*D_i;
Check_1 = A_i_1;
Check_sqrtpinvBi =  inv(D_i')*Check_1*inv(D_i);
Check_pinvBi = Check_sqrtpinvBi^2;
Check_Bi = pinv(Check_pinvBi);
Check_Bimiddle = inv(D_i')*Check_Bi*inv(D_i);

Check_XMX = -Check_Bimiddle+theta_i;
Check_theta = inv(I_H_i)*Check_Bimiddle*inv(I_H_i');




Check_1a = A_i_1a;
Check_sqrtpinvBia =  inv(D_i')*Check_1a*inv(D_i);
Check_pinvBia = Check_sqrtpinvBia^2;
Check_Bia = pinv(Check_pinvBia);
Check_Bimiddlea = inv(D_i')*Check_Bia*inv(D_i);

sum(sum((Check_1-Check_1a).^2))
sum(sum((Check_sqrtpinvBi-Check_sqrtpinvBia).^2))
sum(sum((Check_pinvBi-Check_pinvBia).^2))
sum(sum((Check_Bi-Check_Bia).^2))
sum(sum((Check_Bimiddle-Check_Bimiddlea).^2))