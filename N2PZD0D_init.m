% KM - part of the N2PZD0D model's initialization file

%NPZD in Nitrogen units
%Nn0   = 3.0;  %Nutrient DIN [mmolN*m-3] KM: important for limitation and CNP
Nn0   = 0.5;  %Nutrient DIN [mmolN*m-3] KM: important for limitation and CNP
P1n0  = 0.75; %Phytoplankton type 1 [mmolN*m-3]
P2n0  = 0.75; %Phytoplankton type 2
Zn0   = 0.25; %Zooplankton [mmolN*m-3]
Dn0   = 0.00; %Detrits PON [mmolN*m-3]
DOMn0 = 30/RedCN; %DON [mmolN*m-3];not Redfield, use Letscher's DOM stoichiometry to set

%NPZD in Phosphorus units - start out with Redfield relative to N
Np0   = Nn0/RedNP; %PO4 [mmolP*m-3]
P1p0  = P1n0/RedNP;
P2p0  = P2n0/RedNP;
Zp0   = Zn0/RedNP;
Dp0   = Dn0/RedNP; %Detritus P=POP
DOMp0 = 30/RedCP; %not Redfield, use Letscher's DOM stoichiometry to set

%NPZD in Carbon units
Nc0   = 2000; % DIC; mean surface ocean DIC~2000 umol/kg [mmolC/m3]
P1c0  = P1n0*RedCN;
P2c0  = P2n0*RedCN;
Zc0   = Zn0*RedCN;
Dc0   = Dn0*RedCN; %Detritus C=POC
DOMc0 = 30; % warm surface ocean labile/semilabile DOC umol/kg [mmolC/m3]

%........................................................................
% 4 (NPZD N only)--> 15 (P,C,DOM)--> 18(P2)
V0 = [Nn0;P1n0;P2n0;Zn0;Dn0;DOMn0;
      Np0;P1p0;P2p0;Zp0;Dp0;DOMp0;
      Nc0;P1c0;P2c0;Zc0;Dc0;DOMc0]; 

%........................................................................
ntot0 = sum(V0(1:6));   %initial total N mass [mmolN*m-3]
ptot0 = sum(V0(7:12));  %initial total P mass [mmolN*m-3]
ctot0 = sum(V0(13:18)); %initial total C mass [mmolN*m-3]

%........................................................................
% The following arrays that end with "_t" are calculated in the 
% function N2PZD0D_eqs during ode45 integration. ode45 uses adaptive
% time stepping. The returned values of ode45 were interpolated at 
% the span times. Typically, ode45 integrates multiple times
% within a single "jday." Thus, the values in the arrays are
% not an average of the "jday" but the last appearance. At the 
% same time, some "jday" maybe skipped altogether. This actually happened 
% in one of the runs. The value for C2P_uptake_t, for example, 
% remained zero as initialized for the skipped "jday." 

c2p_uptake_t  = zeros(length(tspan),2); % phytoplankton uptake at each dt
c2n_uptake_t  = zeros(length(tspan),2);
n2p_uptake_t  = zeros(length(tspan),2); 
nplimit_phyto = zeros(length(tspan),2); 
nplimit_zoo   = zeros(length(tspan),2);

gppn_t = zeros(length(tspan),2); % GPPn each dt
gppp_t = zeros(length(tspan),2);
gppc_t = zeros(length(tspan),2);

gspn_t = zeros(length(tspan),2); % GSPn each dt
gspp_t = zeros(length(tspan),2);
gspc_t = zeros(length(tspan),2);

gspn_adj_t = zeros(length(tspan),2); % GSPn each dt
gspp_adj_t = zeros(length(tspan),2);
gspc_adj_t = zeros(length(tspan),2);

gtpn_t = zeros(length(tspan),1); % GSPn each dt
gtpp_t = zeros(length(tspan),1);
gtpc_t = zeros(length(tspan),1);

epn_t = zeros(length(tspan),2); % GSPn each dt
epp_t = zeros(length(tspan),2);
epc_t = zeros(length(tspan),2);

ezn_t = zeros(length(tspan),2); % GSPn each dt
ezp_t = zeros(length(tspan),2);
ezc_t = zeros(length(tspan),2);

mpn_t = zeros(length(tspan),2); % GSPn each dt
mpp_t = zeros(length(tspan),2);
mpc_t = zeros(length(tspan),2);

mzn_t = zeros(length(tspan),1); % GSPn each dt
mzp_t = zeros(length(tspan),1);
mzc_t = zeros(length(tspan),1);

ddn_t = zeros(length(tspan),1); % GSPn each dt
ddp_t = zeros(length(tspan),1);
ddc_t = zeros(length(tspan),1);

ddomn_t = zeros(length(tspan),1); % GSPn each dt
ddomp_t = zeros(length(tspan),1);
ddomc_t = zeros(length(tspan),1);

pf_t = zeros(length(tspan),2); % pf each dt


return