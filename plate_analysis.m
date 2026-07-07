%% Lab 8: Forced Oscillations (spring-mass system with small paper plate)
% Physics 326, Section 3
%
% PURPOSE OF THIS SCRIPT
% -----------------------
% Same analysis pipeline as Lab 7, but repeated for the "plate" damping
% configuration (a much lighter 10.40 g paper plate instead of the large
% damping disc), and extended from 5 to 11 different motor voltages.
% It also includes an extra "mass only, no plate/disc" driven trial
% (partg*.txt) used for the x(t) sinusoidal-growth plot (Figure 14).
%
% Structure mirrors Lab7_Analysis.m:
%   Part A  - free oscillation, used to get measured omega_0
%   Part C  - driven oscillation: fit damping R/gamma, amplitude A,
%             phase phi, v_max, and Q for each of the 11 voltages
%   Part G  - driven oscillation with no damping element, to show
%             amplitude growth over time (near resonance)
% See the Lab 7 script's comments for detailed derivations of the
% physics (Eq. 1-3); comments here focus on what differs for Lab 8
% (11 voltages instead of 5, vectorized loops over columns, etc.)

% ---- Physical constants and system parameters ----
k=18.506;          % spring constant (N/m), same value as Lab 7 (same spring)
mSprEff=0.0289;    % effective mass of the spring (kg)
g=9.80665;         % gravitational acceleration (m/s^2), unused directly

% ---- Part A (f1/f2/f3): free oscillation trials for the plate setup ----
fin=load('partf1.txt');
tA1=fin(:,1);
xA1=fin(:,2);
vA1=fin(:,3);
accA1=fin(:,4);
m_d=0.0104;        % mass of the small paper plate (kg) -- much lighter than Lab 7's disc

m_FA=0.5;          % mass mA (500 g) added to the weight holder for this setup

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

% ---- Driven oscillation trials, 11 motor voltages (partf1_1 .. partf1_11) ----
% Loaded in a loop and horizontally concatenated so each column
% corresponds to one voltage trial (unlike Lab 7 which listed each
% trial's variables individually since there were only 5).
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
    omega=[omega,fin(:,6)];        % column 6 = motor's driving angular velocity
end


x=x-0.821;        % subtract sensor offset so x=0 is the equilibrium position


% Motor voltage readings from two separate multimeter/trial readings
% (used for figure legends below; the tiny discrepancies between
% voltage1 and voltage2 reflect measurement noise/rounding).
voltage1=[4.13,4.59,5.10,5.38,5.46,5.60,5.83,6.37,6.85,7.37,7.98];
voltage2=[4.13,4.58,5.11,5.38,5.46,5.61,5.83,6.39,6.85,7.37,7.98];
omega8=omega;      % keep an unmodified copy of the measured omega matrix for later re-use
%%
% ================== PART A: measured vs. theoretical omega_0 ==================
% Same zero-crossing / half-period method as in Lab 7: find times where
% velocity changes sign (turning points of x), keep every OTHER
% crossing (skipping the half-period crossing), then average the time
% between crossings 2 apart to get the period T, and omega = 2*pi/T.

% --- f1: mA (500g) + plate ---
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
omega_0A1=sqrt(k/(m_FA+0.01+mSprEff));   % theoretical omega_0 (mA + ~10g plate + spring)

% --- f2: mA (500g) only, no plate ---
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

% --- f3: mA + plate + equivalent extra mass (m_d) ----
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
% ================== PART C: driven (forced) oscillation analysis ==================
omega_0=omega_A3;                 % use experimentally measured omega_0 (A3 case)
m=m_d+m_FA+mSprEff                % total effective oscillating mass (note: no semicolon -> prints to console)
t0=pi/2./omega;
v_C=round(v,5);                   % round velocity to 5 decimals to reduce floating-point noise
                                   % when later comparing values with find()/tolerances
t_C=t;


F_0=m*max(a);                     % driving force amplitude estimate (per column/voltage)
omega8=omega;
gamma=[];
% Solve Eq. 1 for R at every sample, then average per voltage column,
% excluding any Inf values (which occur when v_C is exactly 0, causing
% division by zero -- these are treated as outliers and dropped).
R=abs((m*a-F_0.*sin(omega.*t_C)+k*x)./v_C);
for i=1:size(R,2)
    tmp=R(:,i);
    RnoOutlier=tmp(find(tmp~=Inf));
    %disp(mean(RnoOutlier));
    gamma=[gamma,mean(RnoOutlier./2/m)];   % gamma = R/2m, averaged per voltage
end

omega=abs(omega);
% Driven damped-oscillator predictions (same formulas as Lab 7):
%   A(omega)   = (F0/m) / sqrt((omega_0^2-omega^2)^2 + 4*gamma^2*omega^2)
%   phi(omega) = atan( 2*gamma*omega / (omega^2 - omega_0^2) )
A=(F_0./m./((omega_0^2-omega.^2).^2+4*gamma.^2.*omega.^2).^(1/2));
phi=atan(2.*gamma.*omega./(omega.^2-omega_0^2));
tmp=(abs(A-A/sqrt(2)));
omegaHalf=zeros(11,1);

for i=1:11
    [minVal,minTmp]=min(tmp(:,i),[],1);
    omegaHalf(i)=omega(minTmp,i);
end

color=['b','r','y','m','g'];
fprintf("%.4f ",gamma);fprintf("\n");   % print fitted gamma for all 11 voltages

%%
% ---- Figure 8: Amplitude/mass vs. driving frequency, 11 voltages ----
figure;
plot(omega,A);
xlabel('omega rad/s');
ylabel('A/m');
title('A/m vs omega rad/s for 11 voltages applied to motor');
legend(strcat('voltages/V:',num2str(voltage1)));
for i=1:11
    text(abs(omega(50+i,i)),A(50+i,i),strcat('\gamma=',num2str(gamma(i))));
end

%%
% ---- Figure 9/10: v_max vs. log(omega) (Fig 10 = zoomed-in view) ----
v_max=omega.*A;
[sortedMat,indices]=sort(omega);   % sorted omega kept for potential future use (not directly plotted)
figure;
plot(log(omega),v_max);
xlabel('log omega rad/s');
ylabel('v_{max} m/s');
title('v_{max} m/s vs log(omega) rad/s for 5/11 voltages applied to motor');
for i=1:11
    [val,ind]=max(v_max(:,i));
    text(log(omega(ind,i)),v_max(ind,i),strcat('\gamma=',num2str(gamma(i))));
end
hold on;
tmp=max(v_max);
%{
  Reference-line annotation code from Lab 7 (yline at 0.7*peak),
  disabled here since with 11 overlapping curves the lines/labels
  became too cluttered to read.
for i=1:11
    yline(0.7*tmp(i),'LineWidth',1);
    text(-4+i/5,0.7*tmp(i),strcat('gamma=',num2str(gamma(i))));
end
legend('voltages/V: 5.54','4.89','5.75','3.47','7.01');
%}
%%
% ---- Figure 11: phase shift vs. driving frequency ----
figure;
plot(omega,phi);

xlabel('omega rad/s');
ylabel('phase shift/rad');
title('phase shift/rad vs omega rad/s for 5 voltages applied to motor');
for i=1:11
    [val,ind]=min(phi(:,i));
    % label near the most negative phi value for each voltage; the
    % ind-1+mod(ind,2) offset nudges the label to a nearby index so it
    % doesn't sit exactly on top of the curve's minimum point
    text(omega(ind-1+mod(ind,2),i),phi(ind-1+mod(ind,2),i),strcat('\gamma=',num2str(gamma(i))));
end

%%
% ================== Quality factor Q from half-power bandwidth (11 voltages) ==================
% More robust version of the Lab 7 bandwidth search: instead of just
% taking the single closest sample to 0.7*v_peak, this collects ALL
% sample indices within a small tolerance (0.01) of the 0.7*peak level
% on both the rising (Vpk=max) and falling (Vpk=min) sides, then picks
% the pair of indices (one from each side) that are CLOSEST together in
% time/index -- i.e. adjacent half-power crossings around a single
% resonance-like peak, rather than crossings from unrelated regions of
% the noisy signal.
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
    tmp=1000;     % large initial "closest distance" sentinel
    ti=100;tj=100;
    for i=1:length0(k)
        
        for j=1:length1(k)
            % find the closest pair of (rising-side, falling-side)
            % crossing indices, excluding a candidate matching itself
            if abs(Inds(i,k)-Inds1(j,k))<tmp && Inds(i,k)~=Inds1(j,k)
                tmp=abs(Inds(i,k)-Inds1(j,k));
                ti=Inds(i,k);tj=Inds1(j,k);
            end
        end        
    end
    %fprintf('%d %d\n',ti,tj);
    d2omega=[d2omega,abs(omega(ti,k)-omega(tj,k))];   % half-power bandwidth for voltage k
end
Q=omega_0/2./d2omega;             % Q = omega_0 / (2*Delta_omega), Eq. 2
%%

% ================== Uncertainty propagation (per voltage) ==================
errF_0=m*std(a)/(length(a));
omegaErr=std(omega)/sqrt(length(omega));
errD2omega=2*omegaErr;
errQ=sqrt((omegaErr/2./d2omega).^2+(omega_0./d2omega.*errD2omega).^2);
omegaPNErr=sqrt(2)*std(omega);
uncK=3.1144e-04;                  % uncertainty in spring constant k, from prior lab's fit
uncAcc=std(a)/sqrt(length(a));
uncV=std(v_C)/sqrt(length(v_C));
uncX=std(x)/sqrt(length(x));
Rerr=[];

% Unlike Lab 7 (which used mean() across the whole matrix in one line),
% Lab 8 loops explicitly per voltage (i) and per sample (j), summing the
% 6 individual error-propagation terms from Eq. 3 (T1..T6, corresponding
% to dR/dF0, dR/domega, dR/da, dR/dk, dR/dx, dR/dv respectively) while
% skipping any samples where a term blows up to Inf (e.g. v_C(j,i)=0).
% This is mathematically equivalent to Lab 7's vectorized mean(), just
% written out explicitly and with outlier (Inf) rejection added.
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
    % Divide by 200 (samples per trial) to get the mean contribution of
    % each term, matching Lab 7's mean()-based approach.
    T1=T1/200;T2=T2/200;T3=T3/200;T4=T4/200;T5=T5/200;T6=T6/200;
    Rerr=[Rerr,sqrt(T1^2+T2^2+T3^2+T4^2+T5^2+T6^2)];
end
%{
  Vectorized equivalent (as used in Lab 7), left here for reference/
  comparison against the explicit-loop version above.
Rerr=sqrt((mean(m./v_C).*uncAcc).^2+ ...
    (mean(sin(omega.*t_C)./v_C.*errF_0).^2+ ...
    (mean(F_0.*t_C.*cos(omega.*t_C)./v_C).*omegaErr).^2)+ ...
    (mean(x./v_C*uncK)).^2+ ...
    ((mean(k./v_C).*uncX).^2)+ ...
    (mean((m*a-F_0.*sin(omega.*t_C)+k*x)./v_C.^2.*uncV).^2)...
    );
%}
gammaErr=Rerr/2/m;                % propagate R's uncertainty into gamma = R/2m
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
% ================== Figure 12: v_max, actual vs theoretical overlay ==================
t=tA;
figure;
plot(log(omega),v_max);
hold on;
% Dense theoretical omega grid (0 to omega_0) for smooth overlay curves
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

% ================== Figure 13: phase shift, actual vs theoretical overlay ==================
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
% ================== Part G: driven oscillation with mass only (no damper) ==================
% Loads 11 trials (partg1..partg11.txt) of forced oscillation with just
% the hanging mass and no disc/plate attached, so that damping is
% minimal and the resonance build-up (growing amplitude near omega_0)
% is directly visible in x(t) -- see Figure 14 in the writeup.
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
x=x-0.823;         % subtract sensor offset for this configuration
figure;
plot(t,x);
xlabel('t/s');ylabel('x/m');title('x/m vs t/s for forced oscillation with only mass');
figure;
plot(t,omega);
% Estimate the (barely damped) oscillation frequency from velocity
% zero-crossings, same method as Part A above.
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
% ================== Normalized response curves (f, phi, g vs omega/omega_0) ==================
% Same normalization scheme as Lab 7: collapse all 11 voltage trials
% onto a single dimensionless curve using omega_n = omega/omega_0.
omega=abs(omega8)/omega_0;
f=1./sqrt((1-omega.^2).^2+(gamma.*omega).^2);
%figure;
plot(omega,f,'g:','LineWidth',1);
m=0.51;                            % re-estimated mass, matching Lab 7's convention here
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
% ---- Normalized phase shift vs normalized omega ----
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
% ---- Normalized v_max (g function) vs log(normalized omega) ----
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
