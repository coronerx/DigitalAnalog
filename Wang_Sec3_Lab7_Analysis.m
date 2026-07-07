%% Lab 7: Forced Oscillations (spring-mass system with damping disc)
% Physics 326, Section 3
%
% PURPOSE OF THIS SCRIPT
% -----------------------
% Loads position/velocity/acceleration data recorded by the ultrasonic
% position sensor (via Logger Pro) for a mass hanging from a spring,
% under three "free oscillation" configurations (Part A) and five
% "driven oscillation" configurations at different motor voltages
% (Part B/C). From these it:
%   1. Estimates the natural (undamped) angular frequency omega_0 for
%      each free-oscillation configuration, and compares it to the
%      theoretical omega_0 = sqrt(k/m).
%   2. For the driven case, fits the damping resistance R (and hence the
%      damping rate gamma = R/2m) from the equation of motion
%           m*a = F0*sin(omega*t) - k*x - R*v          ... (Eq. 1)
%      solved for R:
%           R = -(m*a - F0*sin(omega*t) + k*x) / v
%   3. Computes the steady-state amplitude A(omega), the velocity
%      extrema v_max(omega), and the phase lag phi(omega) predicted by
%      the driven damped harmonic oscillator model, and compares them to
%      the measured values.
%   4. Estimates the quality factor Q = omega_0 / (2*Delta_omega) from
%      the half-power (0.7 x peak velocity) bandwidth.
%   5. Propagates measurement uncertainty through R, gamma, and Q.
%
% DATA FILES
% ----------
% partA*.txt   - free oscillation (no motor), columns: [t, x, v, a]
% partB1-*.txt - driven oscillation at 5 different motor voltages,
%                columns: [t, x, v, a, ?, omega, alpha] (col 6 = omega,
%                col 7 = angular acceleration, from Logger Pro's motion
%                fit)
%
% NOTE ON CODE QUALITY: This script was written for a single one-off lab
% analysis, so it re-uses variable names (e.g. "vA", "peaks", "Ta") and
% relies on execution order/side effects rather than functions. Comments
% below explain intent, but the control flow itself has not been
% refactored so it stays faithful to the version that produced the
% figures/numbers in the writeup.

% ---- Physical constants and system parameters ----
k=18.506;          % spring constant (N/m), measured in a previous lab (Hooke's law fit)
mSprEff=0.0289;    % effective mass of the spring itself (kg), from previous lab
g=9.80665;         % gravitational acceleration (m/s^2), unused directly but kept for reference

% ---- Part A1: free oscillation, mass + damping disc ----
fin=load('partA1.txt');
tA1=fin(:,1);      % time (s)
xA1=fin(:,2);      % position (m), sensor reading
vA1=fin(:,3);      % velocity (m/s), Logger Pro derivative of x
accA1=fin(:,4);    % acceleration (m/s^2), Logger Pro second derivative of x
m_d=0.1500;        % mass of the damping disc (kg)
m_plate=0.0104;    % mass of the small paper plate used later in Lab 8 (kg)
m_a=0.3;           % mass mA added to the weight holder (kg)
m_FA=0.5;          % total hanging mass used in Lab 8 setup (kg), for reference

% Additional repeated trials for Part A1 (mass + disc), concatenated below
% for more reliable period estimation
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

% ---- Part A2: free oscillation, mass mA only (350 g), no disc ----
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

% ---- Part A3: free oscillation, mass mA + mass equivalent to disc (150 g) ----
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

% ---- Part B: driven oscillation, mass + disc, 5 motor voltages ----
% Columns 6 and 7 (omega, alpha) come from Logger Pro's angular-motion
% fit of the driving arm, i.e. the motor's angular velocity/acceleration.
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
% ================== PART A: measured vs. theoretical omega_0 ==================
% Strategy: find times where the velocity crosses zero (sign change),
% which occur at the turning points (extrema) of x(t). Every OTHER
% zero-crossing corresponds to a full period, so we take the time
% difference between crossings 2 apart (peaks(i+2)-peaks(i)) to estimate
% the period T, then average over many periods and convert to angular
% frequency omega = 2*pi/T.

% --- A1: mass + disc ---
vA=[vA1;vA11;vA12];             % concatenate repeated trials
peaks=[];                        % times of velocity zero-crossings (turning points of x)
tA=[tA1;tA11;tA12];
flag=true;                       % used to only record every OTHER crossing (full period, not half)
for i=2:length(vA)-1
    if flag&&(vA(i)<0 && vA(i+1)>0 || vA(i)>0 && vA(i+1)<0)
        peaks=[peaks,tA(i)];     % record crossing time
        flag=false;              % skip the next crossing (half period)
    else
        flag=true;               % re-arm for the crossing after next
    end
end
Ta=[];                            % measured periods
for i=1:2:length(peaks)-2
    if peaks(i+2)-peaks(i)>0
        Ta=[Ta,peaks(i+2)-peaks(i)];
    end
end

omega_A1=2*pi/mean(Ta);                          % measured angular frequency for A1
omega_0A1=sqrt(k/(.35+mSprEff+m_d));             % theoretical omega_0 = sqrt(k/m_total)

% --- A2: mass mA only (no disc) ---
% NOTE: vA here is built with commas (row concat) unlike A1's semicolons
% (column concat) above -- this looks like a copy/paste inconsistency in
% the original script, but is left as-is since it reproduces the
% reported results.
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
omega_0A2=sqrt(k/(.35+mSprEff));    % theoretical omega_0, no disc mass added

% --- A3: mass mA + mass equivalent to disc (150 g) ---
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
omega_0A3=sqrt(k/(.35+mSprEff+.15));   % theoretical omega_0 with 150 g equivalent mass

%%
% ================== PART C: driven (forced) oscillation analysis ==================
% Combine the 5 voltage trials into single matrices, one column per
% voltage, so all subsequent operations (fits, plots) run on all 5 at
% once via vectorized MATLAB operations.
t_C=[tB11,tB12,tB13,tB14,tB15];
omega_0=omega_A3;                 % use the experimentally measured omega_0 (A3 case) rather
                                   % than the theoretical value, per writeup discussion
m=m_d+m_a;                        % total oscillating mass (disc + mA)
x=[xB11,xB12,xB13,xB14,xB15];
x=x-0.871;                        % subtract sensor offset so x=0 is the equilibrium position
omega=abs([omegaB11,omegaB12,omegaB13,omegaB14,omegaB15]);  % driving angular frequency (motor)
a=[accB11,accB12,accB13,accB14,accB15];
t0=pi/2./omega;                   % quarter-period time scale (not used further below)
F=(m_d+m_a).*[accB11,accB12,accB13,accB14,accB15];  % raw inertial force m*a (unused downstream)
v_C=[vB11,vB12,vB13,vB14,vB15];
omega7=omega;                     % keep an unmodified copy of omega for later re-use/plots
F_0=m*max(a);                     % driving force amplitude, estimated as m * (max acceleration)

% Solve Eq. 1 for the damping resistance R at every sample point:
%   R = |m*a - F0*sin(omega*t) + k*x| / v
R=abs((m*a-F_0.*sin(omega.*t_C)+k*x)./[vB11,vB12,vB13,vB14,vB15]);
gamma=mean(R./2/(m_d+m_a));       % damping rate gamma = R/(2m), averaged over all samples/voltages

% Steady-state driven-oscillator amplitude and phase lag predicted by
% the standard damped-driven harmonic oscillator solution:
%   A(omega)   = (F0/m) / sqrt((omega_0^2-omega^2)^2 + 4*gamma^2*omega^2)
%   phi(omega) = atan( 2*gamma*omega / (omega^2 - omega_0^2) )
A=(F_0./m./((omega_0^2-omega.^2).^2+4*gamma.^2.*omega.^2).^(1/2));
phi=atan(2.*gamma.*omega./(omega.^2-omega_0^2));

% Locate the half-power point (amplitude = A_peak/sqrt(2)) per voltage,
% used as a secondary/alternative bandwidth estimate (not the one
% ultimately reported for Q; see Vpk-based method further below).
tmp=(abs(A-A/sqrt(2)));
omegaHalf=zeros(5,1);

for i=1:5
    [minVal,minTmp]=min(tmp(:,i),[],1);
    omegaHalf(i)=omega(minTmp,i);
end

color=['b','r','y','m','g'];   % plot colors, one per voltage trial
%%

% ---- Figure 2: Amplitude/mass vs. driving frequency ----
figure;
plot(omega,A);
xlabel('omega rad/s');
ylabel('A/m');
title('A/m vs omega rad/s for 5 voltages applied to motor');
legend('voltages/V: 5.54','4.89','5.75','3.47','7.01');
for i=1:5
    % annotate each curve with its fitted damping rate gamma
    text(omega(150,i),A(150,i),strcat('\gamma=',num2str(gamma(i))));
end
%%
% ---- Figure 3: v_max vs. log(omega), with 0.7*peak reference lines ----
% v_max = omega * A is the velocity amplitude at resonance-like extrema;
% used later to estimate the half-power bandwidth for Q.
v_max=omega.*A;
vOrig=v_max;                      % keep a copy since 'v_max'/'omega' get reused/overwritten below
figure;
plot(log(omega),v_max);
xlabel('log omega rad/s');
ylabel('v_{max} m/s');
title('v_{max} m/s vs log(omega) rad/s for 5 voltages applied to motor');
legend('voltages/V: 5.54','4.89','5.75','3.47','7.01');
for i=1:5
    text(log(omega(150,i)),v_max(150,i),strcat('\gamma=',num2str(gamma(i))));
end
hold on;
tmp=max(v_max);
for i=1:5
    % horizontal reference line at 70% of peak v_max -- crossing points
    % of this line with the curve define the half-power bandwidth 2*d_omega
    yline(0.7*tmp(i),color(i),'LineWidth',1);
    text(-4,0.7*tmp(i),strcat('gamma=',num2str(gamma(i))));
end
legend('voltages/V: 5.54','4.89','5.75','3.47','7.01');
%%
% ---- Figure 4: phase shift vs. driving frequency ----
figure;
plot(omega,phi);

xlabel('omega rad/s');
ylabel('phase shift/rad');
title('phase shift/rad vs omega rad/s for 5 voltages applied to motor');
legend('voltages/V: 5.54','4.89','5.75','3.47','7.01');
for i=1:5
    text(omega(150,i),phi(150,i),strcat('\gamma=',num2str(gamma(i))));
end

%%
% ================== Quality factor Q from half-power bandwidth ==================
% Find the driving frequency omegaP where the measured velocity first
% reaches 70% of its maximum value (rising side of resonance) and
% omegaN where it reaches 70% of its minimum (falling side). The
% difference d2omega = omegaP - omegaN approximates the full width
% between the two half-power points, i.e. 2*Delta_omega in Eq. 2.
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
d2omega=(omegaP-omegaN);          % 2*Delta_omega, the half-power bandwidth
Q=omega_0/2./d2omega;             % quality factor, Eq. 2: Q = omega_0 / (2*Delta_omega)
%%

% ================== Uncertainty propagation ==================
% Standard-error-of-the-mean style uncertainties for the key measured
% quantities, then propagated through R (Eq. 3), gamma = R/2m, and Q.
errF_0=m*std(a)/(length(a));
omegaErr=std(omega)/sqrt(length(omega));
errD2omega=2*omegaErr;
errQ=sqrt((omegaErr/2./d2omega).^2+(omega_0./d2omega.*errD2omega).^2);
omegaPNErr=sqrt(2)*std(omega);
uncK=3.1144e-04;                  % uncertainty in spring constant k, carried over from prior lab's fit
uncAcc=std(a)/sqrt(length(a));
uncV=std(v_C)/sqrt(length(v_C));
uncX=std(x)/sqrt(length(x));

% Error propagation for R via Eq. 3:
%   dR/dF0    = sin(omega*t)/v
%   dR/domega = F0*t*cos(omega*t)/v
%   dR/da     = m/v
%   dR/dk     = x/v
%   dR/dv     = (m*a - F0*sin(omega*t) + k*x) / v^2
% (the dR/dx term, x/v * uncX-analog, is combined into the k/v term below)
Rerr=sqrt((mean(m./v_C).*uncAcc).^2+ ...
    (mean(sin(omega.*t_C)./v_C.*errF_0).^2+ ...
    (mean(F_0.*t_C.*cos(omega.*t_C)./v_C).*omegaErr).^2)+ ...
    (mean(x./v_C*uncK)).^2+ ...
    ((mean(k./v_C).*uncX).^2)+ ...
    (mean((m*a-F_0.*sin(omega.*t_C)+k*x)./v_C.^2.*uncV).^2)...
    );
gammaErr=Rerr/2/m;                % propagate R's uncertainty into gamma = R/2m
disp(Rerr);disp(gammaErr);
%{
  Diagnostic printouts of each individual term in the Rerr sum, left
  commented out; useful for identifying which term dominates the total
  uncertainty (per the writeup, the dR/dv term dominates because v
  appears squared in the denominator and v varies a lot near zero
  crossings).
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
% ================== Figure 6: v_max with theoretical overlay ==================
omega=omega7;                     % restore the actual measured omega values
v_max=vOrig;                      % restore actual measured v_max
figure;
plot(log(omega),v_max);
hold on;
% Build a dense, evenly-spaced omega grid from 0 to omega_0 for each
% voltage trial, to plot a smooth theoretical curve alongside the
% (sparser, noisier) actual data.
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
% ================== Figure 7: phase shift, actual vs theoretical ==================
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
% Motor voltage readings recorded from two multimeters/trials (used in
% Lab 8's legend rather than here; kept for reference/consistency).
voltage1=[4.13,4.59,5.10,5.38,5.46,5.60,5.83,6.37,6.85,7.37,7.98];
voltage2=[4.13,4.58,5.11,5.38,5.46,5.61,5.83,6.39,6.85,7.37,7.98];
%%
% ================== Normalized response curves (f, phi, g vs omega/omega_0) ==================
% Normalize the driving frequency by omega_0 and plot the dimensionless
% resonance curve f = 1/sqrt((1-omega_n^2)^2 + (gamma*omega_n)^2), where
% omega_n = omega/omega_0. This collapses all 5 voltage trials onto a
% common curve determined only by their (differing) gamma values.
omega=abs(omega7)/omega_0;
f=1./sqrt((1-omega.^2).^2+(gamma.*omega).^2);
%figure;
plot(omega,f,'b:','LineWidth',2);
m=0.51;                            % re-estimate total mass (slightly different from m_d+m_a above)
F_0=m*max(a);

R=abs((m*a-F_0.*sin(omega.*t_C)+k*x)./v_C);

hold on;

omega=[];
for i=1:5
    omega=[omega;linspace(0,1,200)];
end
omega=omega';
omegaT=omega;                      % keep the smooth theoretical omega grid for later plots
f=1./sqrt((1-omega.^2).^2+(gamma.*omega).^2);

plot(omega,f,'r:','LineWidth',2);
title('normalized theoretical and actual f vs. omega');xlabel('omega rad/s');ylabel('f');
%%
% ---- Normalized phase shift vs normalized omega ----
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
% ---- Normalized v_max (g function) vs log(normalized omega) ----
% zeta is the dimensionless damping ratio (zeta = gamma / (2*omega_0));
% g(omega_n) is the normalized velocity-response curve used to compare
% the disc and plate configurations on the same axes in Lab 8.
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
