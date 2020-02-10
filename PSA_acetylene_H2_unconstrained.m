%% OPTIMIZATION OF SKARSTROM PRESSURE SWING PROCESSES
%% Case 1: Removal of (A) acetylene from (B) H2 using activated carbon
%% (ethylene production)
clear all; clc; clf;

%% Thermodynamic variables
global R T PL cp cv kappa;
R=8.31446; % gas constant, J/(mol K)
T=298.15; % temperature, K
PL=101325; % low pressure, Pa
cp=14.32; % heat capacity at constant pressure, kJ/(kg K), *mass* basis
cv=10.16; % heat capacity at constant volume, kJ/(kg K), *mass* basis
kappa=cp/cv; % heat capacity ratio, dimensionless

%% Adsorber column
global L Lmargin V D area S E theta psi alower aupper Pmargin;
Lmargin=0.05; % robustness margin, m
D=0.1; % column diameter, m
area=pi*D^2/4; % column cross-sectional area, m^2
S=20*1000*6895; % maximum allowable stress of material used, Pa
E=0.70; % welded joint efficiency, dimensionless
theta=3.5e-3; % column wall thickness, m
psi=3e-3; % corrosion allowance, m
alower=2; % lower practical design limit of aspect ratio, dimensionless
aupper=10; % upper practical design limit of aspect ratio, dimensionless
Pmargin=0.9; % margin PH/MAWP (max allowable working P), dimensionless

%% Adsorbent bed
global rhoads eps nu Nads;
rhoads=2100; % adsorbent density, kg/m^3
eps=0.62; % adsorbent bed voidage, dimensionless
nu=(1-eps)/eps; % dimensionless packing ratio, dimensionless
Nads=1000; %adsorbent lifetime, number of regenerations (cycles)

%% Flow variables
global ufeed vfeed Qfeed phi;
ufeed=0.45; % superficial feed velocity, m/s
vfeed=ufeed/eps; % interstiatial feed velocity, m/s
Qfeed=area*ufeed; % superficial feed volumetric flow rate, m^3/s
phi=1.05; % purge ratio, dimensionless

%% Feed variables
global yAfeed;
yAfeed=0.005; % molar fraction of adsorbing species A, dimensionless

%% Adsorption data
global betaA betaB;
HA=0.0356*rhoads; % linear adsorption constant, dimensionless
betaA=1+nu*HA; % adsorption factor for A, dimensionless
betaB=1; % adsorption factor for B, dimensionless (B doesn't adsorb)

%% Compressor data
global etacomp tlife maxratio Pmaxdiff;
etacomp=0.8; % compressor efficiency, dimensionless
tlife=60^2*24*365*10; % compressor lifetime, dimensionless
maxratio=4; % maximum desirable compression ratio, dimensionless
Pmaxdiff=10e5; % maximum allowable pressure differential, Pa

%% Prices & costs
global KB Kelec Kads Ka Kb Kc K1 K2 K3;
KB=0.003; % retail price of B, USD/mol
Kelec=3.6e-8; % electricity tariff, USD/J
Kads=1.20; % adsorbent unit cost, USD/kg
% compressor cost correlating parameters, dimensionless
Ka=[5.8e5,2.6e5]; Kb=[2e4,2.7e3]; Kc=[0.6,0.75];
% column cost correlating parameters, dimensionless
K1=3.5565; K2=0.3776; K3=0.0905;

%% Optimisation type
% set true for PH/PL and L/D to be operationally constrained
% by safety / practical design limits
constrained = false;

%% Optimisation problem formulation and solving via fmincon
fun = @(x) objfun(x);
nonlcon = @(x) constraints(x);
x0 = [10,58,10];
A = [];
b = [];
Aeq = [];
beq = [];
if constrained
    lb = [1,...
        0,...
        max(alower,Lmargin/D)];
    ub = [1/PL*...
            min(Pmaxdiff+PL,...
            Pmargin*min(...
                2*S*E*theta/(D+1.2*theta),...
                4*S*E*theta/(D-0.8*theta))...
            ),...
        tlife,...
        aupper];
else
    Lmargin=0; % Redundant in the absence of
               % safety / practical design limits
    lb = [1,0,Lmargin/D];
	ub = [Inf,tlife,Inf];
end
options = optimoptions('fmincon','Algorithm','sqp','Display','iter',...
    'StepTolerance',1e-12,'ConstraintTolerance',1e-9);
[x,fval] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)

%% Evaluate the remaining variables at optimal point
Pswing = x(1); tfeed = x(2); aspect = x(3);

L=aspect.*D; % column length, m
V=L.*area; % column volume, m^3
PH=Pswing.*PL; % high pressure, Pa

PB=PL.*Qfeed.*tfeed./(R.*T)...
    .*(Pswing.*(1-yAfeed)...
    -phi.*(1-yAfeed.*(Pswing).^(1-betaB./betaA))); % production per cycle, mol

Ncomp=max(1,log(Pswing)./log(maxratio)); % number of compression stages

Wcomp=Ncomp/etacomp.*kappa/(kappa-1).*R.*T...
    .*(Pswing.^((kappa-1)./(kappa.*Ncomp))-1); % work done by compressor, J/mol

Pcomp=Wcomp.*(PH./(R.*T).*Qfeed); % power required by compressor, W

tpres=(1-1./Pswing).*V./Qfeed; % time for pressurisation, s

tcycle=2.*(tfeed+tpres); % total cycle time, s

Revenue=PB.*KB.*tlife./tcycle; % production revenue of B, USD

Celec=Pcomp.*(tlife./2).*Kelec; % cost of electricy, USD

% cost of compressor, USD
Ccomp=min(Ka+Kb.*(Pcomp./1000).^Kc);

% cost of column, USD
Ccolbase=10.^(K1+K2.*log10(V)+K3.*(log10(V)).^2);
Fcol=(PH.*D./(2.*S.*E-1.2.*PH)+psi)/(3.*D./1000+psi);
Ccol=Fcol.*Ccolbase;

% cost of adsorbent, USD
Cads=Kads.*rhoads.*(1-eps).*V.*tlife./(Nads.*tcycle);

cost=-Revenue+Celec+Ccomp+Ccol+Cads; % overall cost (should be negative!)

%% Display results
format bank
results.AspectRatio=x(3);
results.PressureSwing=x(1);
results.ProductionRate=PB;
results.FeedVelocity=vfeed;
results.CompressorWork=Wcomp;
results.CompressorPower=Pcomp;
results.PresTime=tpres;
results.FeedTime=x(2);
results.CycleTime=tcycle;
results.RevenueB=Revenue;
results.ElectricityCost=Celec;
results.CompressorCapital=Ccomp;
results.ColumnCapital=Ccol;
results.AdsorbentCost=Cads;
results.Profit=-cost;
results
format

%% Plot PH/PL vs tfeed global surface
figure(1); clf; hold on;
clearvars xcoord ycoord zcoord fcoord;
xcoord = linspace(1,300); % Pswing
ycoord = linspace(1,100); % tfeed
zcoord = x(3); % L/D
for i=1:length(xcoord)
    for j=1:length(ycoord)
        fcoord(j,i) = objfun([xcoord(i),ycoord(j),zcoord]);
    end
end
%surf(xcoord,ycoord,-fcoord,'EdgeColor','none');
levels=[0:0.2:1.2].*1e7;
contour(xcoord,ycoord,-fcoord,levels,'k--','ShowText','on');
contour(xcoord,ycoord,-fcoord,[0,0],'k-','ShowText','on');
xlabel('P_H/P_L'), ylabel('t_{feed} / s'), zlabel('Lifetime profit/USD');
plot([x(1) x(1)],[x(2) x(2)],'k*');

adsfrontcon = xcoord.^(-betaB./betaA).*betaA.*(L-Lmargin)./vfeed;
desfrontcon = (phi.^(betaA./betaB)).*ones(length(ycoord));
moleboundcon = (1./yAfeed).^(1./(1-betaB./betaA)).*ones(length(ycoord));
plot(xcoord,adsfrontcon,'r-',desfrontcon,ycoord,'b-',...
    moleboundcon,ycoord,'c-');
legend('Lifetime profit/USD','Profit boundary','Optimal point');
title(strcat('Acetylene/H_2, Lifetime profit (L/D=',num2str(x(3)),')'));

hold off;

%% Plot PH/PL vs L/D global surface
figure(2); clf; hold on;
clearvars xcoord ycoord zcoord fcoord;
xcoord = linspace(1,300); % Pswing
ycoord = x(2); % tfeed
zcoord = linspace(1,1e3); % L/D
for i=1:length(xcoord) % Pswing
    for j=1:length(zcoord) % L/D
        fcoord(j,i) = objfun([xcoord(i),ycoord,zcoord(j)]);
    end
end
%surf(xcoord,zcoord,-fcoord,'EdgeColor','none');
levels=[-0.075,-0.05,-0.025,0,0.1,0.2,0.6,1].*1e7;
contour(xcoord,zcoord,-fcoord,levels,'k--','ShowText','on');
contour(xcoord,zcoord,-fcoord,[0,0],'k-','ShowText','on');
xlabel('P_H/P_L'), ylabel('L/D'), zlabel('Lifetime profit/USD');
plot([x(1) x(1)],[x(3) x(3)],'k*');

adsfrontcon = (xcoord.^(betaB./betaA).*vfeed.*x(2)./betaA+Lmargin)./D;
desfrontcon = (phi.^(betaA/betaB)).*ones(length(zcoord));
moleboundcon = (1./yAfeed).^(1./(1-betaB./betaA)).*ones(length(zcoord));
plot(xcoord,adsfrontcon,'r-',desfrontcon,zcoord,'b-',...
    moleboundcon,zcoord,'c-');
legend('Lifetime profit/USD','Profit boundary','Optimal point');
title(strcat('Acetylene/H_2, Lifetime profit (t_{feed} / s=',num2str(x(2)),')'));

hold off;

%% Plot tfeed vs L/D global surface
figure(3); clf; hold on;
clearvars xcoord ycoord zcoord fcoord;
xcoord = x(1); % Pswing
ycoord = linspace(1,100); % tfeed
zcoord = linspace(1,50); % L/D
for i=1:length(ycoord)
    for j=1:length(zcoord)
        fcoord(j,i) = objfun([xcoord,ycoord(i),zcoord(j)]);
    end
end
%surf(ycoord,zcoord,-fcoord,'EdgeColor','none');
levels=[0,0.5,0.8,0.9,1,1.05,1.075].*1e7;
contour(ycoord,zcoord,-fcoord,levels,'k--','ShowText','on');
contour(ycoord,zcoord,-fcoord,[0,0],'k-','ShowText','on');
xlabel('t_{feed} / s'), ylabel('L/D'), zlabel('Lifetime profit/USD');
plot([x(2) x(2)],[x(3) x(3)],'k*');

adsfrontcon = (vfeed.*ycoord./betaA.*(x(1)).^(betaB./betaA)+Lmargin)./D;
plot(ycoord,adsfrontcon,'r-');
legend('Lifetime profit/USD','Profit boundary','Optimal point');
title(strcat('Acetylene/H_2, Lifetime profit (P_H/P_L=',num2str(x(1)),')'));

hold off;