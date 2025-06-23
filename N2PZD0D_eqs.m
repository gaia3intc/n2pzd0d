function [dXdt]=N2PZD0D_eqs(iTime,V0)

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


%STATE VARIABLES =======================================================
Nn    = V0(1);  %NO3 [mmolN*m-3]
Pn(1) = V0(2);  %Phyplankton type 1 [mmolN*m-3]
Pn(2) = V0(3);  %Phyplankton type 1 [mmolN*m-3]
Zn    = V0(4);  %Zooplankton [mmolN*m-3]
Dn    = V0(5);  %Detritus [mmolN*m-3]
DOMn  = V0(6);  %DON

Np    = V0(7);  %PO4 [mmolP*m-3]
Pp(1) = V0(8);  %P
Pp(2) = V0(9);  %P
Zp    = V0(10); %Z
Dp    = V0(11); %D
DOMp  = V0(12); %DOP

Nc    = V0(13); %DIC [mmolC*m-3]
Pc(1) = V0(14); %Phyplankton Carbon content
Pc(2) = V0(15); %Phyplankton Carbon content
Zc    = V0(16); %Zooplankton Carbon content
Dc    = V0(17); %Detritus Carbon content
DOMc  = V0(18); %DOC

jday = ceil(iTime);
% if mod(jday,100)==0, jday, end

nvars = length(V0);
dXdt = zeros(nvars,1);


%CHECK MASS CONSERVATION ===============================================
ntoti = sum(V0(1:6));
ptoti = sum(V0(7:12));
ctoti = sum(V0(12:18));

ndist = ntoti - ntot0;
pdist = ptoti - ptot0;
cdist = ctoti - ctot0;
ndistmax = 1d-6;

if abs(ndist) > ndistmax
    masscheck = [iTime,ntot0,ntoti,ndist]
    disp('Error!!! N mass is NOT conserved!')
    if abs(pdist) > ndistmax
        masscheck = [iTime,ptot0,ptoti,ndist]
        disp('Error!!! P mass is NOT conserved!')
        if abs(cdist) > ndistmax
            masscheck = [iTime,ctot0,ctoti,ndist]
            disp('Error!!! C mass is NOT conserved!')
            pause
        end
        pause
    end
    pause
end


%CHECK FOR NEGATIVE CONCENTRATIONS =====================================
Ineg = find(V0(:) < 0);
if length(Ineg > 0)
    disp('Error!!! 1. there are NEGATIVE concentrations!')
    iTime
    V0
    pause
end


%MODEL TERMS ===========================================================
jpar = parz0(jday);
jtempC = tempC0(jday);
jtempK = jtempC + 273.15; % Kelvin
tresp = q10^((jtempC-10)/10)-q10^((jtempC-32)/3); %Temp response after Blackford 2004
tresp_max = 3.7;  %calculated separately
tresp_min = 0.5;
tresp_range = tresp_max - tresp_min;

%Phytoplankton Equations.................................................
Qpar = (jpar / Isat)*exp(1 - (jpar/Isat)); %Phy light limitation [n.d.] values between 0 and 1.
Qdin = Nn ./ (knp + Nn); %Phytoplankton NO3 limitation [n.d.] values between 0 and 1
Qdip = Np ./ (kpp + Np); %Phytoplankton PO4 limitation [n.d.] values between 0 and 1

%Phytoplankton uptake C:N:P ratio according to TT's meta-analysis; permil
P2Cpermil = P2C0 .* (Np/PO40).^sP2C_PO4.*(Nn/NO30).^sP2C_NO3.*(jtempK/tempK0).^sP2C_temp.*(jpar/par0).^sP2C_I;
N2Cpermil = N2C0 .* (Np/PO40).^sN2C_PO4.*(Nn/NO30).^sN2C_NO3.*(jtempK/tempK0).^sN2C_temp.*(jpar/par0).^sN2C_I;

C2P = 1./(P2Cpermil/1000);
C2N = 1./(N2Cpermil/1000);

C2P_uptake = (C2P<26.6)*26.6 + (C2P>546.7)*546.7 + (C2P>=26.6 & C2P<=546.7).*C2P;
C2N_uptake = (C2N<2)*2 + (C2N>30)*30 + (C2N>=2 & C2N<=30).*C2N;
N2P_uptake = C2P_uptake ./ C2N_uptake;

% % adaptive time stepping in ode45 skipped some jdays
% if ((jday>718) & (jday<720))
%     ['718<jday<720']
%     pause
% end

c2p_uptake_t(jday,:) = C2P_uptake;
c2n_uptake_t(jday,:) = C2N_uptake;
n2p_uptake_t(jday,:) = N2P_uptake;

%Phytoplankton cellular quota C:N:P ratio
C2P_phyto = Pc./Pp;
C2N_phyto = Pc./Pn;
N2P_phyto = Pn./Pp;

%Convert to N units to determine limiting nutrient
for i=1:2
    if (Nn*Qdin(i))<=(Np*Qdip(i)*N2P_uptake(i)) %N limitation
        Qnut(i) = Qdin(i);
        nplimit_phyto(jday,i)=1; %=1 N limitation
        %['N limitation. Qdin: ',num2str(Qdin), ', Qdip: ', num2str(Qdip)]
    else                         %P limitation
        Qnut(i) = Qdip(i);
        nplimit_phyto(jday,i)=2; %=2 P limitation
        %[' P limitation. Qdin: ',num2str(Qdin), ', Qdip: ', num2str(Qdip)]
    end
end


GPPn = (mup*tresp)*Qpar.*Qnut.*Pn;       %Phy gross primary production
EPn = (1- min(1,betap)).*GPPn;           %Phy exudation that increases w/ temperature...5-20%; betap is assimilated
MPn = (mp*tresp).*Pn;                    %Phy natural mortality [mmolN*m-3*d-1]

GPPp = GPPn./N2P_uptake;                 %Linked to P by N:P uptake stoichiometry
EPp = (1- min(1,betap)).*GPPp; 
MPp = MPn./N2P_phyto;                    %Linked to P by phytoplankton N:P              

GPPc = GPPn.*C2N_uptake;                 %Linked to C by C:N uptake stoichiometry
EPc = (1- min(1,betap)).*GPPc; 
MPc = MPn.*C2N_phyto;                    %Linked to P by phytoplankton C:P                            

gppn_t(jday,:) = GPPn;
gppp_t(jday,:) = GPPp;
gppc_t(jday,:) = GPPc;

epn_t(jday,:) = EPn;
epp_t(jday,:) = EPp;
epc_t(jday,:) = EPc;

mpn_t(jday,:) = MPn;
mpp_t(jday,:) = MPp;
mpc_t(jday,:) = MPc;


%Zooplankton Equations...................................................

%Zooplankton cellular quota C:N:P ratio
C2N_zoo = Zc/Zn;
N2P_zoo = Zn/Zp;

%Qphyn = Pn.^npower ./ (kgz^npower + Pn.^npower);    %Zoo grazing limitation [n.d.] for 1P
%Qphyn = Pn / (kgz + sum(Pn)); %Zoo grazing limitation [n.d.] for P>1; npower=1, Hollying Type 2 (small difference)

Qphyn = Pn.^npower / (kgz^npower + sum(Pn.^npower)); %Zoo grazing limitation [n.d.] for P>1
%Qphyn = Pn.^npower / (kgz^npower + sum(Pn)^npower); %Zoo grazing limitation [n.d.] for P>1 WRONG

%Various grazing preferences
%1) Nominal preference (equal=no preference)
pf_nominal = 0.5*ones(size(Qphyn)); 

%2a) Density (linear) dependent grazing (Fasham; denominator is total food)
% e.g., pf_nominal(1)*Pn(1) / (pf_nominal(1)*Pn(1) + pf_nominal(2)*Pn(2))
pf_rho1 = pf_nominal.*Pn/sum(pf_nominal.*Pn);

%2b) Density (squared) dependent grazing preference (as in NPPZ0D model)
% e.g., Pn(1)^2 ./  (Pn(1)^2 + Pn(2)^2); assumes pf_nominal=[0.5 0.5]
% density-squared dependence comes from independent encounter and interaction probabilities
pf_rho2 = (pf_nominal.*Pn).^2 ./ sum( (pf_nominal.*Pn).^2 );

%3) Food quality dependent grazing preference
C2P_optimal = C2N_zoo * N2P_zoo; % optimal CP is zooplankton C:P, assuming gross growth efficiency 
x = (C2P_phyto - C2P_optimal) / C2P_optimal;
fq = normpdf(x,0.0,sigma_pf);  

if enable_pf == 0
    pf_weighted = pf_nominal;  % sum to 1
elseif enable_pf == 1 
%    pf_weighted = pf_rho1;
    pf_weighted = pf_rho2;
elseif enable_pf == 2
    pf_weighted = fq;
elseif enable_pf == 3
    % simple multiplication might not make sense...FQ preference stronger when food is abundant?
    pf_weighted = fq .* pf_rho2; 
end

pf_weighted = pf_weighted./(sum(pf_weighted)); %Normalize so sums to 1

GSPn  = (gzmax*tresp)*pf_weighted.*Qphyn*Zn;  %Zoo grazing (gross secondary production) [mmolN*m-3*d-1]
MZn   = (mz*tresp)*Zn;                        %Zoo natural mortality [mmolN*m-3*d-1]
GTPn  = (cz*tresp)*(Zn*Zn);                   %Higher order predatation on Zoo (gross tertiary prod; closure) [mmolN*m-3*d-1]

GSPp = GSPn./N2P_phyto;                        %GSPp is related to GSPn by phytolankton N:P ratio
MZp = MZn/N2P_zoo;                            %Linked to P by zooplankton N:P ratio
GTPp = GTPn/N2P_zoo;                          %Linked to P by zooplankton N:P ratio

GSPc = GSPn.*C2N_phyto;                        %GSPc is related to GSPn by cellular C:N ratio
MZc = MZn*C2N_zoo;                            %Linked to P by zooplankton N:P ratio
GTPc = GTPn*C2N_zoo;                          %Linked to P by zooplankton N:P ratio

% Zooplankton homeostatic adjustment of GSP fluxes (NOTE: not a direct adjustment of Z biomass)
% First, adjust either N or P whichever is more unbalanced from ideal zoo C:N and C:P
P2C_zoo_ideal = (1/ZooCP_ref)^(1-H) .* (1./C2P_phyto).^H;
N2C_zoo_ideal = (1/ZooCN_ref)^(1-H) .* (1./C2N_phyto).^H;

Nrelease_byZ = max(0, GSPn - GSPc.*N2C_zoo_ideal);
Prelease_byZ = max(0, GSPp - GSPc.*P2C_zoo_ideal);

Nlimit_vsP_1 = GSPn - (N2C_zoo_ideal./P2C_zoo_ideal).*GSPp; %N excess if >0
Plimit_vsN_1 = GSPp - (P2C_zoo_ideal./N2C_zoo_ideal).*GSPn; %P excess if >0

Nadj_flux_1 = Nrelease_byZ.*double(Nlimit_vsP_1>0); %Either N release
Padj_flux_1 = Prelease_byZ.*double(Plimit_vsN_1>0); % or P release, not both

GSPn_1 = GSPn - Nadj_flux_1; %Net change after first adjustment
GSPp_1 = GSPp - Padj_flux_1;

% Second, adjust C flux relative to N or P whichever is more limiting
Cflux_vsNP = max(GSPc-GSPn_1./N2C_zoo_ideal, GSPc-GSPp_1./P2C_zoo_ideal);
Cadj_flux_1 = Cflux_vsNP.*double(Cflux_vsNP>0);

% Third, go back to N and P, and adjust N:P flux
Nlimit_vsP_2 = GSPn_1 - (N2C_zoo_ideal./P2C_zoo_ideal).*GSPp_1; %N excess: P limitation if >0
Plimit_vsN_2 = GSPp_1 - (P2C_zoo_ideal./N2C_zoo_ideal).*GSPn_1; %P excess: N limitation if >0

Nadj_flux_2 = Nlimit_vsP_2.*double(Nlimit_vsP_2>0);
Padj_flux_2 = Plimit_vsN_2.*double(Plimit_vsN_2>0);

nplimit_zoo(jday,:) = 12*(Nlimit_vsP_2>0) + 11*(Plimit_vsN_2>0); %11=N limit; 12=P limit

GSPn_adj = Nadj_flux_1 + Nadj_flux_2; %Total homeostatic adjustment in flux (not biomass)
GSPp_adj = Padj_flux_1 + Padj_flux_2;
GSPc_adj = Cadj_flux_1;

GSPn_2 = GSPn - GSPn_adj; %New GSP after homeostatic adjustment; 
GSPp_2 = GSPp - GSPp_adj;
GSPc_2 = GSPc - GSPc_adj;

EZn = (1- min(1,(betaz)))*GSPn_2; %Z exudation calculated on the basis of adjusted GSP
EZp = (1- min(1,(betaz)))*GSPp_2;
EZc = (1- min(1,(betaz)))*GSPc_2;

gspn_t(jday,:) = GSPn;
gspp_t(jday,:) = GSPp;
gspc_t(jday,:) = GSPc;

gspn_adj_t(jday,:) = GSPn_adj;
gspp_adj_t(jday,:) = GSPp_adj;
gspc_adj_t(jday,:) = GSPc_adj;

gtpn_t(jday) = GTPn;
gtpp_t(jday) = GTPp;
gtpc_t(jday) = GTPc;

ezn_t(jday,:) = EZn;
ezp_t(jday,:) = EZp;
ezc_t(jday,:) = EZc;

mzn_t(jday) = MZn;
mzp_t(jday) = MZp;
mzc_t(jday) = MZc;

pf_t(jday,:) = pf_weighted;


%Detritus Equations......................................................
DDn = (md*tresp)*Dn; %PON degradation rate;
DDp = (md*tresp)*Dp; %POP degradation rate;
DDc = (md*tresp)*Dc; %POC degradation rate;

ddn_t(jday) = DDn;
ddp_t(jday) = DDp;
ddc_t(jday) = DDc;

%DOM Equations...........................................................
DDOMn = (mdom*tresp)*DOMn; %DON degradation rate;
DDOMp = (mdom*tresp)*DOMp; %DOP degradation rate;
DDOMc = (mdom*tresp)*DOMc; %DOC degradation rate;

ddomn_t(jday) = DDOMn;
ddomp_t(jday) = DDOMp;
ddomc_t(jday) = DDOMc;


%ECOSYSTEM MODEL EQUATIONS ==============================================
% Nitrogen NPZD + DOM [mmolN*m-3*d-1]
dXdt(1)  = DDOMn - sum(GPPn); 
dXdt(2)  = GPPn(1) - EPn(1) - GSPn(1) - MPn(1); 
dXdt(3)  = GPPn(2) - EPn(2) - GSPn(2) - MPn(2); 
dXdt(4)  = sum(GSPn) - sum(EZn) - GTPn - MZn - sum(GSPn_adj); 
dXdt(5)  = GTPn + sum((1-epsPhy).*EPn) + (1-epsZoo)*sum(EZn+GSPn_adj) + sum((1-omePhy).*MPn) + (1-omeZoo)*MZn - DDn; 
dXdt(6)  = sum(epsPhy.*EPn) + epsZoo*sum(EZn+GSPn_adj) + sum(omePhy.*MPn) + omeZoo*MZn + DDn - DDOMn;
% Phosphorus NPZD + DOM
dXdt(7)  = DDOMp - sum(GPPp); 
dXdt(8)  = GPPp(1) - EPp(1) - GSPp(1) - MPp(1); 
dXdt(9)  = GPPp(2) - EPp(2) - GSPp(2) - MPp(2); 
dXdt(10) = sum(GSPp) - sum(EZp) - GTPp - MZp - sum(GSPp_adj); 
dXdt(11) = GTPp + sum((1-epsPhy).*EPp) + (1-epsZoo)*sum(EZp+GSPp_adj) + sum((1-omePhy).*MPp) + (1-omeZoo)*MZp - DDp;
dXdt(12) = sum(epsPhy.*EPp) + epsZoo*sum(EZp+GSPp_adj) + sum(omePhy.*MPp) + omeZoo*MZp + DDp - DDOMp;
% Carbon NPZD + DOM
dXdt(13)  = DDOMc - sum(GPPc); 
dXdt(14) = GPPc(1) - EPc(1) - GSPc(1) - MPc(1); 
dXdt(15) = GPPc(2) - EPc(2) - GSPc(2) - MPc(2); 
dXdt(16) = sum(GSPc) - sum(EZc) - GTPc - MZc - sum(GSPc_adj); 
dXdt(17) = GTPc + sum((1-epsPhy).*EPc) + (1-epsZoo)*sum(EZc+GSPc_adj) + sum((1-omePhy).*MPc) + (1-omeZoo)*MZc - DDc;
dXdt(18) = sum(epsPhy.*EPc) + epsZoo*sum(EZc+GSPc_adj) + sum(omePhy.*MPc) + omeZoo*MZc + DDc - DDOMc;


%CHECK MASS CONSERVATION ================================================
sumODEsN = sum(dXdt(1:6));
sumODEsP = sum(dXdt(7:12));
sumODEsC = sum(dXdt(13:18));

xdistmax = 1d-6;
xdistN = abs(sumODEsN - 0d0);
xdistP = abs(sumODEsP - 0d0);
xdistC = abs(sumODEsC - 0d0);

if (xdistN > xdistmax)
    sumODEsN
    disp('Error!!!! Nitrogen Mass is NOT conserved!!!')
    pause
elseif (xdistP > xdistmax)
    sumODEsP
    disp('Error!!!! Phosphorus Mass is NOT conserved!!!')
    pause
elseif (xdistC > xdistmax)
    sumODEsC
    disp('Error!!!! Carbon Mass is NOT conserved!!!')
    pause
end


return