%% OPTIMIZATION OF SKARSTROM PRESSURE SWING PROCESSES
%% Constraints (nonlinear) for fmincon

function [c,ceq] = constraints(x)
global R T PL cp cv kappa;
global L Lmargin V D area S E theta psi alower aupper Pmargin;
global rhoads eps nu;
global ufeed vfeed Qfeed phi;
global yAfeed;
global betaA betaB;
global etacomp tlife maxratio Pmaxdiff;
global KB Kelec Kads Ka Kb Kc K1 K2 K3;

Pswing = x(1); tfeed = x(2); aspect = x(3);

L=aspect*D;
V=aspect*area*D;
PH=Pswing*PL;

c = [tfeed*vfeed/(betaA*(L-Lmargin))-Pswing^(-betaB/betaA);...
    Pswing^(-betaB/betaA)-1/phi;...
    yAfeed*Pswing^(1-betaB/betaA)-1;...
    -yAfeed*Pswing^(1-betaB/betaA)];
ceq = [];
end