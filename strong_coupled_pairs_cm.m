clear k;
%since roots of Fs are on imag axis Fs*=Fs N=9
F22S=-Fs;
%convert to Y matrix
yd=zeros(1,10);
for k=1:10
    if mod(k,2)==0
       yd(k)=i*imag(Es(k)+Fs(k));
    else
       yd(k)=real(Es(k)+Fs(k));
    end
end
k=0;
y11n=zeros(1,9);
y22n=zeros(1,9);
for k=1:9
    if mod(k,2)==0
       y11n(k)=i*imag(Es(k+1)-Fs(k+1));
       y22n(k)=i*imag(Es(k+1)+Fs(k+1));
    else
       y11n(k)=real(Es(k+1)-Fs(k+1));
       y22n(k)=real(Es(k+1)+Fs(k+1));
    end
end
k=0;
y21n=-Ps/E;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yd=yd/2;y11n=y11n/2;y22n=y22n/2;y21n=y21n/2;
%y11=y11n/yd;y12=y21=y21n/yd;y22=y22n/yd
[r11k,pk,k11]=residue(y11n,yd);
[r21k,pk,k21]=residue(y21n,yd);
lamdak=pk/i;
r22k=r11k;
%N+2 transversal matrix M
M=zeros(11,11);
%M(1,1)=k11=0;M(11,11)=k22=0;M111=M111=k21=0;
for k=2:10
    M(k,k)=-lamdak(k-1);
    M(k,11)=sqrt(r22k(k-1));
    M(11,k)=M(k,11);
    M(1,k)=r21k(k-1)/sqrt(r22k(k-1));
    M(k,1)=M(1,k);
    
end
k=0;M=real(M);
%rotation to arrow form %%%%%%%%%%%%%%%%%%%
Mr=M;
%rotation matrix R
R=eye(11,11);
%annihilate element m,n 0-->S 1-9 11-->L
%inner(m,n)-->M(m+1,n+1)
% m=0;n=9;
for m=0:7
    n=9;
    while n-m>=2
%pivotinner(n-1,n)-->pivotM(n,n+1)=pivotM(i,j)
%calculate theta:Mkl/Mmn-->k=m+1,l=n+1,m=m+1,n=n;
          theta_r=-atan(Mr(m+1,n+1)/Mr(m+1,n));
          cr=cos(theta_r);
          sr=sin(theta_r);
          R(n,n)=cr;
          R(n+1,n+1)=cr;
          R(n+1,n)=sr;
          R(n,n+1)=-sr;
          Mr=R*Mr*R';
%           for k=1:11
%               for j=1:11
%                   if abs(Mr(j,k))<0.0001
%                      Mr(j,k)=0;
%                   end
%               end
%           end
%           clear k;
%           clear j;
          R=eye(11,11);
          n=n-1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rotation 1 to trisection %%%%%%%%%%%%%%%%%%%
m=8;n=9;w=TZ(1);%(2,4)
for k=1:6
    theta_r=atan(Mr(m+1,n+1)/(w+Mr(n+1,n+1)));
    cr=cos(theta_r);
    sr=sin(theta_r);
    R(n,n)=cr;
    R(n+1,n+1)=cr;
    R(n+1,n)=sr;
    R(n,n+1)=-sr;
    Mr=R*Mr*R';
    R=eye(11,11);
    m=m-1;
    n=n-1;
end
k=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rotation 2 to trisection %%%%%%%%%%%%%%%%%%%
m=8;n=9;w=TZ(2);%(6,8)
for k=1:2
    theta_r=atan(Mr(m+1,n+1)/(w+Mr(n+1,n+1)));
    cr=cos(theta_r);
    sr=sin(theta_r);
    R(n,n)=cr;
    R(n+1,n+1)=cr;
    R(n+1,n)=sr;
    R(n,n+1)=-sr;
    Mr=R*Mr*R';
    R=eye(11,11);
    m=m-1;
    n=n-1;
end
k=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set main coupling be positive
for k=1:10
    if Mr(k,k+1)<0
       Mr(:,k+1)=-Mr(:,k+1);
       Mr(k+1,:)=-Mr(k+1,:);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%verify coulping matrix M%%%%%%%%%%%%%%%%%%
Gm=zeros(11,11);
Gm(1,1)=1;
Gm(11,11)=1;
Cm=eye(11,11);
Cm(1,1)=0;
Cm(11,11)=0;
syms x;
Zm=Gm+i*x.*Cm+i.*Mr;% M
Zm=inv(Zm);
s21=2*Zm(11,1);
s11=-1+2*Zm(1,1);
% mags11=10*log(abs(s11));
% mags21=10*log(abs(s21));
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
figure(4);
% fplot(mags11,[-5,5],'b');hold on;
% fplot(mags21,[-5,5],'r');
plot(sample,10*log(abs(subs(s21,x,sample))));
hold on;
plot(sample,10*log(abs(subs(s11,x,sample))));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mr=Mr/(f0/BW);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=zeros(1,9);
for k=1:9
       f(k)=0.5*(2-Mr(k+1,k+1))*f0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%self coupling
