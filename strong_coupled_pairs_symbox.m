%trisection to symbox
%1st SCRP transformation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%trisection to box
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f0=1250.4;
f_even=880;%assign appropriate f_even
k_14=Mr(2+1,4+1);%assign coupling coefficient
k_12=Mr(2+1,3+1);%assign coupling coefficient
k_24=Mr(3+1,4+1);%assign coupling coefficient
f_1=f(2);f_2=f(3);f_4=f(4);%assign frequency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
omega_even=(f_even/f0-f0/f_even);
k_13=sqrt(abs(omega_even*k_14));
k_34=sign(omega_even*k_14)*k_13;
f_3=f_even;
omega_1=(f_1/f0-f0/f_1)+k_13^2/omega_even;
omega_4=(f_4/f0-f0/f_4)+k_34^2/omega_even;
f_1=f0*(omega_1/2+sqrt((omega_1/2)^2+1));
f_4=f0*(omega_4/2+sqrt((omega_4/2)^2+1));
%box to symbox
omega_2=(f_2/f0-f0/f_2);
omega_3=(f_3/f0-f0/f_3);
omega_N=(omega_2+omega_3)/2;
f_2=f0*(omega_N/2+sqrt((omega_N/2)^2+1));
f_3=f_2;%f2=f3 and f1=f4 do not change
M_23=(omega_3-omega_2)/2;
M_12=(k_12-k_13)/sqrt(2);
M_13=(k_12+k_13)/sqrt(2);
M_24=(k_24-k_34)/sqrt(2);
M_34=(k_24+k_34)/sqrt(2);
% M_12 M_24 M_23 positive
M_12=abs(M_12);M_24=abs(M_24);M_23=abs(M_23);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2nd SCRP transformation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%trisection to box
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f0=1250.4;
f_even=930;%assign appropriate f_even
k_14=Mr(6+1,8+1);%assign coupling coefficient
k_12=Mr(6+1,7+1);%assign coupling coefficient
k_24=Mr(7+1,8+1);%assign coupling coefficient
f_1=f(6);f_2=f(7);f_4=f(8);%assign frequency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
omega_even=(f_even/f0-f0/f_even);
k_13=sqrt(abs(omega_even*k_14));
k_34=sign(omega_even*k_14)*k_13;
f_3=f_even;
omega_1=(f_1/f0-f0/f_1)+k_13^2/omega_even;
omega_4=(f_4/f0-f0/f_4)+k_34^2/omega_even;
f_1=f0*(omega_1/2+sqrt((omega_1/2)^2+1));
f_4=f0*(omega_4/2+sqrt((omega_4/2)^2+1));
%box to symbox
omega_2=(f_2/f0-f0/f_2);
omega_3=(f_3/f0-f0/f_3);
omega_N=(omega_2+omega_3)/2;
f_2=f0*(omega_N/2+sqrt((omega_N/2)^2+1));
f_3=f_2;%f2=f3 and f1=f4 do not change
M_23=(omega_3-omega_2)/2;
M_12=(k_12-k_13)/sqrt(2);
M_13=(k_12+k_13)/sqrt(2);
M_24=(k_24-k_34)/sqrt(2);
M_34=(k_24+k_34)/sqrt(2);
% M_12 M_24 M_23 positive
M_12=abs(M_12);M_24=abs(M_24);M_23=abs(M_23);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f0=1250.4;
% f_even=950.5;%assign appropriate f_even
% k_14=-0.0195;%assign coupling coefficient
% k_12=0.0159;%assign coupling coefficient
% k_24=0.0159;%assign coupling coefficient
% f_1=1252.76;f_2=1238.26;f_4=1252.76;%assign frequency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% omega_even=(f_even/f0-f0/f_even);
% k_13=sqrt(abs(omega_even*k_14));
% k_34=sign(omega_even*k_14)*k_13;
% f_3=f_even;
% omega_1=(f_1/f0-f0/f_1)+k_13^2/omega_even;
% omega_4=(f_4/f0-f0/f_4)+k_34^2/omega_even;
% f_1=f0*(omega_1/2+sqrt((omega_1/2)^2+1));
% f_4=f0*(omega_4/2+sqrt((omega_4/2)^2+1));
% %box to symbox
% omega_2=(f_2/f0-f0/f_2);
% omega_3=(f_3/f0-f0/f_3);
% omega_N=(omega_2+omega_3)/2;
% f_2=f0*(omega_N/2+sqrt((omega_N/2)^2+1));
% f_3=f_2;%f2=f3 and f1=f4 do not change
% M_23=(omega_3-omega_2)/2;
% M_12=(k_12-k_13)/sqrt(2);
% M_13=(k_12+k_13)/sqrt(2);
% M_24=(k_24-k_34)/sqrt(2);
% M_34=(k_24+k_34)/sqrt(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%