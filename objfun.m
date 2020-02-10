%% OPTIMIZATION OF SKARSTROM PRESSURE SWING PROCESSES
%% Objective function

function cost = objfun(x)
global R T PL cp cv kappa;
global L Lmargin V D area S E theta psi alower aupper Pmargin;
global rhoads eps nu Nads;
global ufeed vfeed Qfeed phi;
global yAfeed;
global betaA betaB;
global etacomp tlife maxratio Pmaxdiff;
global KB Kelec Kads Ka Kb Kc K1 K2 K3;

Pswing = x(1); tfeed = x(2); aspect = x(3);

L=aspect.*D; % column length, m
V=L.*area; % column volume, m^3
PH=Pswing.*PL; % high pressure, Pa

PB=PL.*Qfeed.*tfeed./(R.*T)...
    .*(Pswing.*(1-yAfeed)...
    -phi.*(1-yAfeed.*(Pswing).^(1-betaB./betaA))); % production per cycle, mol

Ncomp=max(1,log(Pswing)./log(maxratio));

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
end