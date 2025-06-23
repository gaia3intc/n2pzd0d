more off
%close all
clear all
format short g

% diary on

%------------------------------------------------------------------------
% MATLAB / GNU OCTAVE -- PROGRAM "MarinEcosysModellingCourse_NPZD0D.m" 
%------------------------------------------------------------------------
% This script solves the basic NPZD model in 0D (single water volume):
%
% dP/dt =  Fp - Ep - Gp - Mp; %Phytoplankton
% dZ/dt =  Fz - Ez - Gz - Mz; %Zooplankton
% dN/dt = -Fp + (epsPhy*Ep + epsZoo*Ez) + (omePhy*Mp + omeZoo*Mz) + md*D; %Nitrates
% dD/dt = (1-epsPhy)*Ep + (1-epsZoo)*Ez + (1-omePhy)*Mp + (1-omeZoo)*Mz - md*D; %Detritus
%
%------------------------------------------------------------------------
% AUTHORS 
% 
% Sergio M. Vallina (24 April 2022) - basic NPZD used in summer school
%
% Katsumi Matsumoto, U. Minnesota (April 2024) 
% (a) added Q10 temperature dependence to kinetic reaction rates
% (b) added vars: PO4, DIC, DOM (C, N, P) and D (P, C) w/o full C chemistry
% (c) linked phytoplankton C:N:P stoichiometry by the power law formula
% (d) added second P (phytoplankton); (possibly another Z) 
% (e) implemented Z grazing that depends on population and food quality
% (f) carry out experiments under different temperature regimes
%       - standing biomass, DOM and POM stoichiometry, nutrient cycling? 
%------------------------------------------------------------------------


% GLOBAL VARIABLES:------------------------------------------------------
global mup knp kpp betap mp Isat 
global gzmax kgz betaz mz cz 
global epsPhy omePhy epsZoo omeZoo md mdom
global ntot0 ptot0 ctot0
global parz0 tempC0 q10 npower
global P2C0 sP2C_PO4 sP2C_NO3 sP2C_temp sP2C_I
global N2C0 sN2C_PO4 sN2C_NO3 sN2C_temp sN2C_I
global PO40 NO30 tempK0 par0  
global H ZooCP_ref ZooCN_ref 
global enable_pf sigma_pf
global c2p_uptake_t c2n_uptake_t n2p_uptake_t nplimit_phyto nplimit_zoo
global gppn_t gppp_t gppc_t gspn_t gspp_t gspc_t gtpn_t gtpp_t gtpc_t 
global epn_t epp_t epc_t ezn_t ezp_t ezc_t ddn_t ddp_t ddc_t
global mpn_t mpp_t mpc_t mzn_t mzp_t mzc_t ddomn_t ddomp_t ddomc_t
global gspn_adj_t gspp_adj_t gspc_adj_t pf_t


% TEMPORAL RESOLUTION:---------------------------------------------------
% The time step (dt) value does not matter for the "ode45" solver, which
% uses its own internal adaptative dt.
dt = 1;     %time step [days] for output resolution. 
t0 = dt;
days = 360;
years = 3;
tmax = days*years;
tspan = t0:dt:tmax;


%MODEL INITIALIZATION AND PARAMETERS:-----------------------------------
N2PZD0D_params
N2PZD0D_init


%FILENAME AND DIARY:-----------------------------------------------------
filename = [date '_Nn' num2str(Nn0) '_pf' num2str(enable_pf) '_sigma' num2str(sigma_pf)'']
eval(['diary ' filename '.diary'])


%EXTERNAL FORCINGS:------------------------------------------------------
T = days;
tau = days/4;
time = dt:dt:tmax;
ttime = time-tau;
wrad = (2*pi/T);

parmin = 5; %winter [W*m-2]
parmax = 140; %summer [W*m-2]

A0 = parmin;
A = (parmax-parmin)/2; %amplitud.
parsin = A0 + A*( sin(wrad*ttime) + 1); %Sinusoidal PAR
parz0 = parsin; %PAR at the surface [W*m-2]

tempCmin = 12; %winter [deg C]; tresp=1 when temp 10C
tempCmax = 24; %summer [deg C]; set to 10 to get Vallina solution

temp0 = tempCmin;
tempC = (tempCmax-tempCmin)/2;                 %amplitude
tempC0 = temp0 + tempC*( sin(wrad*ttime) + 1); %sinusoidal temperature


%SOLVE EQUATIONS USING RUNGE-KUTTA ODE45:---------------------------------
ode45options = odeset('AbsTol',1e-12,'RelTol',1e-6);
%ode45options = odeset('AbsTol',1e-12,'RelTol',1e-6,'NonNegative',[1 3]);
[tode,Vode] = ode45(@N2PZD0D_eqs,tspan,V0,ode45options);

% Vode = Vode'; %transpose matrix [rows-StateVariables, columns-time]


%STATE VARIABLES:--------------------------------------------------------

yr2start = length(tspan)/years+1;
tplot = tode(yr2start:end); %Exclude the first year

NO3 = Vode(yr2start:end,1);
PO4 = Vode(yr2start:end,7);
DIC = Vode(yr2start:end,13);

P1n = Vode(yr2start:end,2);
P1p = Vode(yr2start:end,8);
P1c = Vode(yr2start:end,14);

P2n = Vode(yr2start:end,3);
P2p = Vode(yr2start:end,9);
P2c = Vode(yr2start:end,15);

Zn = Vode(yr2start:end,4);
Zp = Vode(yr2start:end,10);
Zc = Vode(yr2start:end,16);

Dn = Vode(yr2start:end,5);
Dp = Vode(yr2start:end,11);
Dc = Vode(yr2start:end,17);

DON = Vode(yr2start:end,6);
DOP = Vode(yr2start:end,12);
DOC = Vode(yr2start:end,18);

[V0 Vode(end,:)']

ntotf = sum(Vode(end,1:6));   %final total N mass [mmolN*m-3]
ptotf = sum(Vode(end,7:12));  %final total P mass [mmolN*m-3]
ctotf = sum(Vode(end,13:18)); %final total C mass [mmolN*m-3]

['start, end total C:P = ' num2str(ctot0/ptot0) ', ' num2str(ctotf/ptotf)]
['start, end total C:N = ',num2str(ctot0/ntot0) ', ' num2str(ctotf/ntotf)]
['start, end total N:P = ',num2str(ntot0/ptot0) ', ' num2str(ntotf/ptotf)]


%SAVE and PLOT:---------------------------------------------------------

N2PZD0D_plots


return



% figure(1), print -dpng sinfig1_Nn3_pf1.png
% figure(2), print -dpng sinfig2_Nn3_pf1.png
% figure(3), print -dpng sinfig3_Nn3_pf1.png
% figure(4), print -dpng sinfig4_Nn3_pf1.png
% figure(5), print -dpng sinfig5_Nn3_pf1.png

% figure(1), print -dpng sinfig1_Nn3_pf2_sigma1.png
% figure(2), print -dpng sinfig2_Nn3_pf2_sigma1.png
% figure(3), print -dpng sinfig3_Nn3_pf2_sigma1.png
% figure(4), print -dpng sinfig4_Nn3_pf2_sigma1.png
% figure(5), print -dpng sinfig5_Nn3_pf2_sigma1.png
% 
% 
% clear *0 sN2C* sP2C* *min *max Red* *ref *totf ode45* *pf umol*
% clear mup knp kpp betap mp Isat epsPhy omePhy fignum
% clear gzmax kgz betaz mz cz epsZoo omeZoo md mdom
% clear A H T ans days dt npower tau tempC wrad years yr2start
% 
% diary off
% %eval(['save ./Results/' filename ' *'])
% eval(['save ' filename '.mat *'])
% 
% 
% return