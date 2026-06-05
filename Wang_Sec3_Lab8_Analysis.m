k=18.506;
mSprEff=0.0289;
g=9.80665;
fin=load('partf1.txt');
tA1=fin(:,1);
xA1=fin(:,2);
vA1=fin(:,3);
accA1=fin(:,4);
m_d=0.0104;

m_FA=0.5;

fin=load('partf2.txt');
tA2=fin(:,1);
xA2=fin(:,2);
vA2=fin(:,3);
accA2=fin(:,4);
fin=load('partf3.txt');
tA3=fin(:,1);
xA3=fin(:,2);
vA3=fin(:,3);
accA3=fin(:,4);
t=[];
x=[];
a=[];
omega=[];
v=[];
for i=1:11
    fin=load(strcat('partf1',int2str(i),'.txt'));
    t=[t,fin(:,1)];
    x=[x,fin(:,2)];
    v=[v,fin(:,3)];
    a=[a,fin(:,4)];
    omega=[omega,fin(:,6)];
end


x=x-0.821;


voltage1=[4.13,4.59,5.10,5.38,5.46,5.60,5.83,6.37,6.85,7.37,7.98];
voltage2=[4.13,4.58,5.11,5.38,5.46,5.61,5.83,6.39,6.85,7.37,7.98];
omega8=omega;
%%
% part A
vA=vA1;
peaks=[];
tA=tA1;
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
omega_0A1=sqrt(k/(m_FA+0.01+mSprEff));

vA=vA2;
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
omega_0A2=sqrt(k/(m_FA+mSprEff));

vA=vA3;
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
omega_0A3=sqrt(k/(m_FA+mSprEff+m_d));
%%
% part C
omega_0=omega_A3;
m=m_d+m_FA+mSprEff
t0=pi/2./omega; 
v_C=round(v,5);
t_C=t;


F_0=m*max(a);
omega8=omega;
gamma=[];
R=abs((m*a-F_0.*sin(omega.*t_C)+k*x)./v_C);
for i=1:size(R,2)
    tmp=R(:,i);
    RnoOutlier=tmp(find(tmp~=Inf));
    %disp(mean(RnoOutlier));
    gamma=[gamma,mean(RnoOutlier./2/m)];
end

omega=abs(omega);
A=(F_0./m./((omega_0^2-omega.^2).^2+4*gamma.^2.*omega.^2).^(1/2));
phi=atan(2.*gamma.*omega./(omega.^2-omega_0^2));
tmp=(abs(A-A/sqrt(2)));
omegaHalf=zeros(11,1);

for i=1:11
    [minVal,minTmp]=min(tmp(:,i),[],1);
    omegaHalf(i)=omega(minTmp,i);
end

color=['b','r','y','m','g'];
fprintf("%.4f ",gamma);fprintf("\n");

%%
figure;
plot(omega,A);
xlabel('omega rad/s');
ylabel('A/m');
title('A/m vs omega rad/s for 11 voltages applied to motor');
legend(strcat('voltages/V:',num2str(voltage1)));
for i=1:11
    text(abs(omega(50+i,i)),A(50+i,i),strcat('γ=',num2str(gamma(i))));
end

%%

v_max=omega.*A;
[sortedMat,indices]=sort(omega);
figure;
plot(log(omega),v_max);
xlabel('log omega rad/s');
ylabel('v_{max} m/s');
title('v_{max} m/s vs log(omega) rad/s for 5/11 voltages applied to motor');
for i=1:11
    [val,ind]=max(v_max(:,i));
    text(log(omega(ind,i)),v_max(ind,i),strcat('γ=',num2str(gamma(i))));
end
hold on;
tmp=max(v_max);
%{
for i=1:11
    yline(0.7*tmp(i),'LineWidth',1);
    text(-4+i/5,0.7*tmp(i),strcat('gamma=',num2str(gamma(i))));
end
legend('voltages/V: 5.54','4.89','5.75','3.47','7.01');
%}
%%
figure;
plot(omega,phi);

xlabel('omega rad/s');
ylabel('phase shift/rad');
title('phase shift/rad vs omega rad/s for 5 voltages applied to motor');
for i=1:11
    [val,ind]=min(phi(:,i));
    text(omega(ind-1+mod(ind,2),i),phi(ind-1+mod(ind,2),i),strcat('γ=',num2str(gamma(i))));
end

%%
Vpk=0.7*max(v_C);
omegaP=zeros(size(Vpk));
omegaN=zeros(size(Vpk));
Inds=zeros(30,11);
length0=[];
for i=1:11
    tmp=abs(v_C(:,i) - Vpk(:,i));
    Ind=find(tmp<0.01);
    for j=1:length(Ind)
        Inds(j,i)=Ind(j);
    end
    length0=[length0,length(Ind)];
end
% disp(Vpk);
Vpk=0.7*min(v_C);
Inds1=zeros(30,11);
length1=[];
for i=1:11
    tmp=abs(v_C(:,i) - Vpk(:,i));
    Ind1=find(tmp<0.01);
    for j=1:length(Ind1)
        Inds1(j,i)=Ind1(j);
    end
    length1=[length1,length(Ind1)];
end
d2omega=[];
for k=1:11
    tmp=1000;
    ti=100;tj=100;
    for i=1:length0(k)
        
        for j=1:length1(k)
            if abs(Inds(i,k)-Inds1(j,k))<tmp && Inds(i,k)~=Inds1(j,k)
                tmp=abs(Inds(i,k)-Inds1(j,k));
                ti=Inds(i,k);tj=Inds1(j,k);
            end
        end        
    end
    %fprintf('%d %d\n',ti,tj);
    d2omega=[d2omega,abs(omega(ti,k)-omega(tj,k))];
end
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
Rerr=[];

for i=1:11
    T1=0;T2=0;T3=0;T4=0;T5=0;T6=0;
    for j=1:length(v_C)
        t1=m./v_C(j,i)*uncAcc(i);
        t2=(sin(omega(j,i)*t_C(j,i))/v_C(j,i).*errF_0(i));
        t3=((F_0(i).*t_C(j,i).*cos(omega(j,i).*t_C(j,i))./v_C(j,i)).*omegaErr(i));
        t4=((x(j,i)./v_C(j,i)*uncK));
        t5=(((k./v_C(j,i)).*uncX(i)));
        t6=((m*a(j,i)-F_0(i).*sin(omega(j,i).*t_C(j,i))+k*x(j,i))./v_C(j,i).^2.*uncV(i));
        if [t1, t2, t3, t4, t5, t6]<Inf
            T1=T1+t1;T2=T2+t2;T3=T3+t3;T4=t4+T4;T5=T5+t5;T6=T6+t6;
        end        
        
        
    end
    T1=T1/200;T2=T2/200;T3=T3/200;T4=T4/200;T5=T5/200;T6=T6/200;
    Rerr=[Rerr,sqrt(T1^2+T2^2+T3^2+T4^2+T5^2+T6^2)];
end
%{
Rerr=sqrt((mean(m./v_C).*uncAcc).^2+ ...
    (mean(sin(omega.*t_C)./v_C.*errF_0).^2+ ...
    (mean(F_0.*t_C.*cos(omega.*t_C)./v_C).*omegaErr).^2)+ ...
    (mean(x./v_C*uncK)).^2+ ...
    ((mean(k./v_C).*uncX).^2)+ ...
    (mean((m*a-F_0.*sin(omega.*t_C)+k*x)./v_C.^2.*uncV).^2)...
    );
%}
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
%%
t=tA;
figure;
plot(log(omega),v_max);
hold on;
omega=[];
for i=1:11
    omega=[omega;linspace(0,omega_0,200)];
end
omega=omega';
A=(F_0./m./((omega_0^2-omega.^2).^2+4*gamma.^2.*omega.^2).^(1/2));
v_max=omega.*A;
phiTheory=atan(2.*gamma.*omega./(omega.^2-omega_0^2));

plot(log(omega),v_max,'LineWidth',1);
xlabel('log omega rad/s');
ylabel('v_{max} m/s');
title('v_{max} m/s vs log(omega) rad/s for 11 voltages applied to motor');
hold on;
for i=1:11
    plot(log(omega(:,i)),omega(:,i).*A(:,i),':','LineWidth',3);
end

title('theoretical v_{max} m/s vs log(omega) rad/s');
xlabel('log(omega) rad/s');ylabel('v_{max} m/s');

hold off;

figure;
plot(abs(omega8),phi);

xlabel('omega rad/s');
ylabel('phase shift/rad');
title('phase shift/rad vs omega rad/s');
hold on;
for i=1:11
    plot(omega(:,i),phiTheory(:,i),':','LineWidth',3);
end
title('theoretical phi/rad vs omega rad/s');
xlabel('omega rad/s');ylabel('phi/rad');

hold off;
%%
t=[];
x=[];
a=[];
omega=[];
v=[];
for i=1:11
    fin=load(strcat('partg',int2str(i),'.txt'));
    t=[t,fin(:,1)];
    x=[x,fin(:,2)];
    v=[v,fin(:,3)];
    a=[a,fin(:,4)];
    omega=[omega,fin(:,6)];
end
x=x-0.823;
figure;
plot(t,x);
xlabel('t/s');ylabel('x/m');title('x/m vs t/s for forced oscillation with only mass');
figure;
plot(t,omega);
vA=v;
peaks=[];
tA=t;
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
omegaG=2*pi/mean(Ta);
%%
omega=abs(omega8)/omega_0;
f=1./sqrt((1-omega.^2).^2+(gamma.*omega).^2);
%figure;
plot(omega,f,'g:','LineWidth',1);
m=0.51;
F_0=m*max(a);

R=abs((m*a-F_0.*sin(omega.*t_C)+k*x)./v_C);

hold on;

omega=[];
for i=1:11
    omega=[omega;linspace(0,1,200)];
end
omega=omega';
omegaT=omega;
f=1./sqrt((1-omega.^2).^2+(gamma.*omega).^2);

plot(omega,f,'black:','LineWidth',2);
title('theoretical and actual f vs. omega');xlabel('omega rad/s');ylabel('f');
%%
figure;
omega=abs(omega8)/omega_0;
phiN=atan(gamma/omega_0.*omega./(omega.^2-1));
plot(omega,phiN,'g:','LineWidth',2);

hold on;
omega=omegaT;
phiN=atan(gamma/omega_0.*omega./(omega.^2-1));
plot(omega,phiN,'black:','LineWidth',2);
xlabel('normalized omega');ylabel('normalized phi');title('normalized theoretical and actual phi vs. omega');

%%
omega=abs(omega8)/omega_0;
zeta=gamma/2/omega_0;
g=2*zeta./sqrt((omega-1./omega).^2+4*zeta.^2);
%figure;
plot(log(omega),g,'g:','LineWidth',1.5);
hold on;
omega=[];
for i=1:11
    omega=[omega;linspace(0,1,200)];
end
omega=omega';
g=2*zeta./sqrt((omega-1./omega).^2+4*zeta.^2);
plot(log(omega),g,'black:','LineWidth',2);
xlabel('log normalized omega');ylabel('normalized v_{max}');
title('normalized theoretical and actual v_{max} vs. omega in logarithmic view');