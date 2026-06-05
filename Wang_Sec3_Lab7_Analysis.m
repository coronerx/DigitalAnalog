k=18.506;
mSprEff=0.0289;
g=9.80665;
fin=load('partA1.txt');
tA1=fin(:,1);
xA1=fin(:,2);
vA1=fin(:,3);
accA1=fin(:,4);
m_d=0.1500;
m_plate=0.0104;
m_a=0.3;
m_FA=0.5;

fin=load('partA1-1.txt');
tA11=fin(:,1);
xA11=fin(:,2);
vA11=fin(:,3);
accA11=fin(:,4);
fin=load('partA1-3.txt');
tA12=fin(:,1);
xA12=fin(:,2);
vA12=fin(:,3);
accA12=fin(:,4);

fin=load('partA2.txt');
tA2=fin(:,1);
xA2=fin(:,2);
vA2=fin(:,3);
accA2=fin(:,4);
fin=load('partA2-1.txt');
tA21=fin(:,1);
xA21=fin(:,2);
vA21=fin(:,3);
accA21=fin(:,4);
fin=load('partA2-2.txt');
tA22=fin(:,1);
xA22=fin(:,2);
vA22=fin(:,3);
accA22=fin(:,4);

fin=load('partA3.txt');
tA3=fin(:,1);
xA3=fin(:,2);
vA3=fin(:,3);
accA3=fin(:,4);

fin=load('partA3-1.txt');
tA31=fin(:,1);
xA31=fin(:,2);
vA31=fin(:,3);
accA31=fin(:,4);

fin=load('partA3-2.txt');
tA32=fin(:,1);
xA32=fin(:,2);
vA32=fin(:,3);
accA32=fin(:,4);

fin=load('partB1-1.txt');
tB11=fin(:,1);
xB11=fin(:,2);
vB11=fin(:,3);
accB11=fin(:,4);
omegaB11=fin(:,6);
alphaB11=fin(:,7);
fin=load('partB1-2.txt');
tB12=fin(:,1);
xB12=fin(:,2);
vB12=fin(:,3);
accB12=fin(:,4);
omegaB12=fin(:,6);
alphaB12=fin(:,7);

fin=load('partB1-3.txt');
tB13=fin(:,1);
xB13=fin(:,2);
vB13=fin(:,3);
accB13=fin(:,4);
omegaB13=fin(:,6);
alphaB13=fin(:,7);

fin=load('partB1-4.txt');
tB14=fin(:,1);
xB14=fin(:,2);
vB14=fin(:,3);
accB14=fin(:,4);
omegaB14=fin(:,6);
alphaB14=fin(:,7);

fin=load('partB1-5.txt');
tB15=fin(:,1);
xB15=fin(:,2);
vB15=fin(:,3);
accB15=fin(:,4);
omegaB15=fin(:,6);
alphaB15=fin(:,7);
%%
% part A
vA=[vA1;vA11;vA12];
peaks=[];
tA=[tA1;tA11;tA12];
flag=true;
for i=2:length(vA)-1
    if flag&&(vA(i)<0 && vA(i+1)>0 || vA(i)>0 && vA(i+1)<0)
        peaks=[peaks,tA(i)];
        flag=false;
    else
        flag=true;
    end
end
Ta=[];
for i=1:2:length(peaks)-2
    if peaks(i+2)-peaks(i)>0
        Ta=[Ta,peaks(i+2)-peaks(i)];
    end
end

omega_A1=2*pi/mean(Ta);
omega_0A1=sqrt(k/(.35+mSprEff+m_d));

vA=[vA2,vA21,vA22];
peaks=[];
for i=2:length(vA)-1
    if flag&&(vA(i)<0 && vA(i+1)>0 || vA(i)>0 && vA(i+1)<0)
        peaks=[peaks,tA(i)];
        flag=false;
    else
        flag=true;
    end
end
Ta=[];
for i=1:2:length(peaks)-2
    if peaks(i+2)-peaks(i)>0
        Ta=[Ta,peaks(i+2)-peaks(i)];
    end
end

omega_A2=2*pi/mean(Ta);
omega_0A2=sqrt(k/(.35+mSprEff));

vA=[vA3,vA31,vA32];
peaks=[];
for i=2:length(vA)-1
    if flag&&(vA(i)<0 && vA(i+1)>0 || vA(i)>0 && vA(i+1)<0)
        peaks=[peaks,tA(i)];
        flag=false;
    else
        flag=true;
    end
end
Ta=[];
for i=1:2:length(peaks)-2
    if peaks(i+2)-peaks(i)>0
        Ta=[Ta,peaks(i+2)-peaks(i)];
    end
end

omega_A3=2*pi/mean(Ta);
omega_0A3=sqrt(k/(.35+mSprEff+.15));

%%
% part C
t_C=[tB11,tB12,tB13,tB14,tB15];
omega_0=omega_A3;
m=m_d+m_a;
x=[xB11,xB12,xB13,xB14,xB15];
x=x-0.871;
omega=abs([omegaB11,omegaB12,omegaB13,omegaB14,omegaB15]);
a=[accB11,accB12,accB13,accB14,accB15];
t0=pi/2./omega; 
F=(m_d+m_a).*[accB11,accB12,accB13,accB14,accB15];
v_C=[vB11,vB12,vB13,vB14,vB15];
omega7=omega;
F_0=m*max(a);

R=abs((m*a-F_0.*sin(omega.*t_C)+k*x)./[vB11,vB12,vB13,vB14,vB15]);
gamma=mean(R./2/(m_d+m_a));
A=(F_0./m./((omega_0^2-omega.^2).^2+4*gamma.^2.*omega.^2).^(1/2));
phi=atan(2.*gamma.*omega./(omega.^2-omega_0^2));
tmp=(abs(A-A/sqrt(2)));
omegaHalf=zeros(5,1);

for i=1:5
    [minVal,minTmp]=min(tmp(:,i),[],1);
    omegaHalf(i)=omega(minTmp,i);
end

color=['b','r','y','m','g'];
%%

figure;
plot(omega,A);
xlabel('omega rad/s');
ylabel('A/m');
title('A/m vs omega rad/s for 5 voltages applied to motor');
legend('voltages/V: 5.54','4.89','5.75','3.47','7.01');
for i=1:5
    text(omega(150,i),A(150,i),strcat('γ=',num2str(gamma(i))));
end
%%
v_max=omega.*A;
vOrig=v_max;
figure;
plot(log(omega),v_max);
xlabel('log omega rad/s');
ylabel('v_{max} m/s');
title('v_{max} m/s vs log(omega) rad/s for 5 voltages applied to motor');
legend('voltages/V: 5.54','4.89','5.75','3.47','7.01');
for i=1:5
    text(log(omega(150,i)),v_max(150,i),strcat('γ=',num2str(gamma(i))));
end
hold on;
tmp=max(v_max);
for i=1:5
    yline(0.7*tmp(i),color(i),'LineWidth',1);
    text(-4,0.7*tmp(i),strcat('gamma=',num2str(gamma(i))));
end
legend('voltages/V: 5.54','4.89','5.75','3.47','7.01');
%%
figure;
plot(omega,phi);

xlabel('omega rad/s');
ylabel('phase shift/rad');
title('phase shift/rad vs omega rad/s for 5 voltages applied to motor');
legend('voltages/V: 5.54','4.89','5.75','3.47','7.01');
for i=1:5
    text(omega(150,i),phi(150,i),strcat('γ=',num2str(gamma(i))));
end

%%
Vpk=0.7*max(v_C);
omegaP=zeros(size(Vpk));
omegaN=zeros(size(Vpk));
for i=1:5
    [~, Ind] = min(abs(v_C(:,i) - Vpk(:,i)));
    omegaP(i)=omega(Ind,i);
end
% disp(Vpk);
Vpk=0.7*min(v_C);
for i=1:5
    [~, Ind] = min(abs(v_C(:,i) - Vpk(:,i)));
    omegaN(i)=omega(Ind,i);
end
d2omega=(omegaP-omegaN);
Q=omega_0/2./d2omega;
%%

errF_0=m*std(a)/(length(a));
omegaErr=std(omega)/sqrt(length(omega));
errD2omega=2*omegaErr;
errQ=sqrt((omegaErr/2./d2omega).^2+(omega_0./d2omega.*errD2omega).^2);
omegaPNErr=sqrt(2)*std(omega);
uncK=3.1144e-04;
uncAcc=std(a)/sqrt(length(a));
uncV=std(v_C)/sqrt(length(v_C));
uncX=std(x)/sqrt(length(x));

Rerr=sqrt((mean(m./v_C).*uncAcc).^2+ ...
    (mean(sin(omega.*t_C)./v_C.*errF_0).^2+ ...
    (mean(F_0.*t_C.*cos(omega.*t_C)./v_C).*omegaErr).^2)+ ...
    (mean(x./v_C*uncK)).^2+ ...
    ((mean(k./v_C).*uncX).^2)+ ...
    (mean((m*a-F_0.*sin(omega.*t_C)+k*x)./v_C.^2.*uncV).^2)...
    );
gammaErr=Rerr/2/m;
disp(Rerr);disp(gammaErr);
%{
fprintf('%.4f ',(mean(m./v_C).*uncAcc).^2);
fprintf('\n');
    fprintf('%.4f ',mean(sin(omega.*t_C)./v_C.*errF_0).^2); 
    fprintf('\n');
    fprintf('%.4f ',(mean(F_0.*t_C.*cos(omega.*t_C)./v_C).*omegaErr).^2); 
    fprintf('\n');
    fprintf('%.4f ',(mean(x./v_C*uncK)).^2);
    fprintf('\n');
    fprintf('%.4f ',(mean(k./v_C).*uncX).^2);
    fprintf('\n');
%}
    fprintf('%.4f ',mean((m*a-F_0.*sin(omega.*t_C)+k*x)./v_C.^2).^2);fprintf('\n');
%%
omega=omega7;
v_max=vOrig;
figure;
plot(log(omega),v_max);
hold on;
omega=[];
for i=1:5
    omega=[omega;linspace(0,omega_0,600)];
end
omega=omega';
A=(F_0./m./((omega_0^2-omega.^2).^2+4*gamma.^2.*omega.^2).^(1/2));
v_max=omega.*A;

plot(log(omega),v_max,'r:','LineWidth',2);
xlabel('log omega rad/s');
ylabel('v_{max} m/s');
title('v_{max} m/s vs log(omega) rad/s for 11 voltages applied to motor');
xlabel('log(omega) rad/s');ylabel('v_{max} m/s');

hold off;
%%
figure;
plot(abs(omega7),phi);

phiTheory=atan(2.*gamma.*omega./(omega.^2-omega_0^2));
xlabel('omega rad/s');
ylabel('phase shift/rad');
title('phase shift/rad vs omega rad/s');
hold on;
for i=1:5
    plot(omega(:,i),phiTheory(:,i),':','LineWidth',3);
end
title('theoretical and actual phi/rad vs omega rad/s');
xlabel('omega rad/s');ylabel('phi/rad');

hold off;
%%
voltage1=[4.13,4.59,5.10,5.38,5.46,5.60,5.83,6.37,6.85,7.37,7.98];
voltage2=[4.13,4.58,5.11,5.38,5.46,5.61,5.83,6.39,6.85,7.37,7.98];
%%
omega=abs(omega7)/omega_0;
f=1./sqrt((1-omega.^2).^2+(gamma.*omega).^2);
%figure;
plot(omega,f,'b:','LineWidth',2);
m=0.51;
F_0=m*max(a);

R=abs((m*a-F_0.*sin(omega.*t_C)+k*x)./v_C);

hold on;

omega=[];
for i=1:5
    omega=[omega;linspace(0,1,200)];
end
omega=omega';
omegaT=omega;
f=1./sqrt((1-omega.^2).^2+(gamma.*omega).^2);

plot(omega,f,'r:','LineWidth',2);
title('normalized theoretical and actual f vs. omega');xlabel('omega rad/s');ylabel('f');
%%
%figure;
omega=abs(omega7)/omega_0;
phiN=atan(gamma/omega_0.*omega./(omega.^2-1));
plot(omega,phiN,'b:','LineWidth',2);

hold on;
omega=[];
for i=1:5
    omega=[omega;linspace(0,1,200)];
end
omega=omega';
phiN=atan(gamma/omega_0.*omega./(omega.^2-1));
plot(omega,phiN,'r:','LineWidth',2);
xlabel('normalized omega');ylabel('normalized phi');title('normalized theoretical and actual phi vs. omega');

%%
omega=abs(omega7)/omega_0;
zeta=gamma/2/omega_0;
g=2*zeta./sqrt((omega-1./omega).^2+4*zeta.^2);
figure;
plot(log(omega),g,'b:','LineWidth',1.5);
hold on;
omega=omegaT;
g=2*zeta./sqrt((omega-1./omega).^2+4*zeta.^2);
plot(log(omega),g,'r:','LineWidth',2);
xlabel('log normalized omega');ylabel('normalized v_{max}');
title('normalized theoretical and actual v_{max} vs. omega in logarithmic view');
