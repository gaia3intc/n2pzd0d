% KM - part of the N2PZD0D model's parameter setting file

mup   = [3.0 0.5];   %Phytoplankton max growth specific rate [d-1] (KM: at 10C)
knp   = [4.0 0.2];   %Phytoplankton NO3 half-saturation uptake constant [mmolN*m-3]
kpp   = [0.4 0.2];   %Phytoplankton PO4 half-saturation uptake constant [mmolp*m-3]
betap = [0.8 0.8];   %Phytoplankton assimilation efficiency [n.d.]
mp    = [0.05 0.05]; %Phytoplankton mortality specific rate [d-1]

Isat  = 120;  %Phytoplankton saturating irradiance [W*m-2]

gzmax = 1.0;  %Zooplankton maximum ingestion specific rate [d-1] (KM: at 10C)
kgz   = 0.75; %Zooplankton half-saturation ingestion constant [mmolN*m-3]
betaz = 0.40; %Zooplankton assimilation efficiency [n.d.]
mz    = 0.05; %Zooplankton mortality specific rate [d-1]
cz    = 0.10; %Zooplankton closure mortality specific rate [d-1]

epsPhy = [1/4 1/4]; %Phyplankton exudation fraction going to DOM [n.d]
epsZoo = 1/4; %Zooplankton exudation fraction going to DOM [n.d]

omePhy = [1/4 1/4]; %Phyplankton mortality fraction going to DOM [n.d]
omeZoo = 1/4; %Zooplankton mortality fraction going to DOM [n.d]

md     = 0.2; %Detritus degradation specific rate [d-1]
mdom   = 0.2; %DOM degradation specific rate [d-1]
q10    = 2.0; %Q10=2.0 (kinetic rates increase x2 for every 10C warming)
npower = 2.0; %Grazing type: 1 for Holling type II; 2 for Holling Type III

umolkg2mmolm3 = (1e-3)/(1e-6)/1025; 

RedCP = 106/1;  %Redfield C:P ratio
RedCN = 106/16; %Redfield C:N ratio
RedNP = 16/1;   %Redfield N:P ratio

H = 0.08;        %TT (2021) GRL 0.08; H=0 100% homeostatic; H=1 you are what you eat
ZooCP_ref = RedCP;
ZooCN_ref = RedCN;

%CNP power law constants (P1=eukaryotes, P2=cyanobacteria; from TM meta-analysis)
P2C0      = [11.6 6.3];    %permil: P/C=11.6 permil==>C:P=86.2; P/C=6.3 permil==>C:P=158.73
sP2C_PO4  = [0.58 0.28];
sP2C_NO3  = [0 0];
sP2C_temp = [0 -8.0];
sP2C_I    = [0 0];

N2C0      = [151.0 151.0]; %permil: N/C=151.0 permil==>C:N=6.62
sN2C_PO4  = [0 0];
sN2C_NO3  = [0.22 0.22];
sN2C_temp = [0 0];
sN2C_I    = [-0.05 -0.05];

PO40 = 0.57*umolkg2mmolm3; % reference PO4
NO30 = 5.70*umolkg2mmolm3; % reference NO3
tempK0 = 291;              % temp in Kelvin (18C)
par0 = 70;                 % reference light level (W/m2)

enable_pf = 2    % 0 = equal preference, fixed
                 % 1 = only prey density (squared) dependent (rho)
                 % 2 = only food quality-dependent (fq)
                 % 3 = both rho and fq dependent

sigma_pf = 1.0     % default = 1.0; option for when enable_pf=2 or 3
                 % selectiveness (smaller value peaks higher w/ narrower distribution)
                 % large sigma_pf --> enable_pf=1
                 % v. small sigma_pf --> preferred prey and thus zooplankton wiped out


return