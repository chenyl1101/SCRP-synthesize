%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
freq=linspace(1000,1400,1000);
f1=1212.5;%MHz
f2=1284.9;
f0=sqrt(f1*f2);
BW=f2-f1;
omega=(f0/BW).*(freq./f0-f0./freq);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TZ=[1203.22 1206.86];%num of TZ=2
TZ=(f0/BW).*(TZ./f0-f0./TZ);
Pw=poly(TZ);
Ps=poly(i*TZ);
Fw=zeros(1,9);%num of order=9
Fw(1)=1;
%derivation of Cw function
ref=[-1 -0.8 -0.6 -0.3 -0.1 0.1 0.3 0.6 0.8 1];
for num=1:6
    %solve the linear equations
    G=[ref(1)^8 ref(1)^7 ref(1)^6 ref(1)^5 ref(1)^4 ref(1)^3 ref(1)^2 ref(1) 1 polyval(Pw,ref(1)),
       ref(2)^8 ref(2)^7 ref(2)^6 ref(2)^5 ref(2)^4 ref(2)^3 ref(2)^2 ref(2) 1 -polyval(Pw,ref(2)),   
       ref(3)^8 ref(3)^7 ref(3)^6 ref(3)^5 ref(3)^4 ref(3)^3 ref(3)^2 ref(3) 1 polyval(Pw,ref(3)),
       ref(4)^8 ref(4)^7 ref(4)^6 ref(4)^5 ref(4)^4 ref(4)^3 ref(4)^2 ref(4) 1 -polyval(Pw,ref(4)),
       ref(5)^8 ref(5)^7 ref(5)^6 ref(5)^5 ref(5)^4 ref(5)^3 ref(5)^2 ref(5) 1 polyval(Pw,ref(5)),
       ref(6)^8 ref(6)^7 ref(6)^6 ref(6)^5 ref(6)^4 ref(6)^3 ref(6)^2 ref(6) 1 -polyval(Pw,ref(6)),
       ref(7)^8 ref(7)^7 ref(7)^6 ref(7)^5 ref(7)^4 ref(7)^3 ref(7)^2 ref(7) 1 polyval(Pw,ref(7)),
       ref(8)^8 ref(8)^7 ref(8)^6 ref(8)^5 ref(8)^4 ref(8)^3 ref(8)^2 ref(8) 1 -polyval(Pw,ref(8)),
       ref(9)^8 ref(9)^7 ref(9)^6 ref(9)^5 ref(9)^4 ref(9)^3 ref(9)^2 ref(9) 1 polyval(Pw,ref(9)),
       ref(10)^8 ref(10)^7 ref(10)^6 ref(10)^5 ref(10)^4 ref(10)^3 ref(10)^2 ref(10) 1 -polyval(Pw,ref(10))];
    V=-[ref(1)^9 ref(2)^9 ref(3)^9 ref(4)^9 ref(5)^9 ref(6)^9 ref(7)^9 ref(8)^9 ref(9)^9 ref(10)^9];
    X=inv(G)*V';
    for k=1:9
        Fw(k+1)=X(k);
    end
    k=0;
   
    [dFw,dPw]=polyder(Fw,Pw);
    new=roots(dFw);
    %choose inner passband extreme point
    for k=1:length(new)
        if abs(real(new(k)))<1 && abs(imag(new(k)))<0.0001
           new(k)=new(k);
        else 
           new(k)=-2;
        end
    end
    k=0;
    %order those point
    new=real(new);%
    new=sort(new);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this part is related to num of TZs
    for k=1:length(new)-1
        if new(k+1)==-2
           ref(k)=ref(k);
        else ref(k)=new(k+1);
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
end
clear k;
% Cw after interation
Cw=poly2sym(Fw)/(poly2sym(Pw)*X(10));
% figure(1);
% fplot(Cw,[-1,1],'r');hold on;
% line([-1,1],[0,0],'linestyle','--');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E=abs(polyval(Pw,[-1])/polyval(Fw,[-1]))/sqrt(10^1.8-1);%RL=-18dB
Fs=poly(roots(Fw)*i);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rEs=roots([zeros(1,7) Pw]./E-i.*Fw);
clear k;
for k=1:length(rEs)
    if imag(rEs(k))<0
       rEs(k)=conj(rEs(k));
    else
       rEs(k)=rEs(k);
    end
end
Ew=poly(rEs);
Es=poly(i*rEs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fs_conj=[-conj(Fs(1)),conj(Fs(2)),-conj(Fs(3)),conj(Fs(4)),-conj(Fs(5)),conj(Fs(6))];
% Ps_conj=[conj(Ps(1)),-conj(Ps(2)),conj(Ps(3)),-conj(Ps(4)),conj(Ps(5))];
% sum=[zeros(1,2) conv(Ps,Ps_conj)./E^2]+conv(Fs,Fs_conj);
% rsum=roots(sum);
% rEs=zeros(1,5);
% j=1;
% for k=1:length(rsum)
%     if real(rsum(k))<0
%        rEs(j)=rsum(k);
%        j=j+1;
% %     else
% %        rEs(j)=rEs(j);
%     end
% end
% clear j;clear k;
% Es=sqrt(-sum(1))*poly(rEs);%conservation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N-nfz is odd
figure(2);
sample=linspace(-2,2,500);
plot(sample,10.*log(abs(polyval(Fs,i.*sample)./polyval(Es,i.*sample))));
hold on;
plot(sample,10.*log(abs(polyval(Ps,i.*sample)./polyval(Es,i.*sample))./E));
figure(3);
plot(freq.*1e6,10.*log(abs(polyval(Fs,i.*omega)./polyval(Es,i.*omega))));
hold on;
plot(freq.*1e6,10.*log(abs(polyval(Ps,i.*omega)./polyval(Es,i.*omega))./E));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%