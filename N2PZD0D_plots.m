fignum = 1;

figure(fignum), clf
subplot(2,2,1)
plot(tplot,P1n,'-g.'), hold on
plot(tplot,P2n,'-m.')
plot(tplot,Zn,'-k.')
plot(tplot,NO3,'-r.')
plot(tplot,Dn,'-b.')
plot(tplot,DON,'-c.')
title('Stocks in terms of nitrogen'),xlabel('Time'), ylabel('Nitrogen (mmolN*m-3)')
%legend('Phyto1','Phyto2','Zoo','NO3','Detritus','DON')
grid on

subplot(2,2,2)
plot(tplot,P1p,'-g.'), hold on
plot(tplot,P2p,'-m.')
plot(tplot,Zp,'-k.')
plot(tplot,PO4,'-r.')
plot(tplot,Dp,'-b.')
plot(tplot,DOP,'-c.')
title('...in terms of phosphorus'),xlabel('Time'), ylabel('Phosphorus (mmolP*m-3)')
legend('Phyto1','Phyto2','Zoo','PO4','Detritus','DOP')
grid on

subplot(2,2,3)
plot(tplot,P1c,'-g.'), hold on
plot(tplot,P2c,'-m.')
plot(tplot,Zc,'-k.')
plot(tplot,DIC/100,'-r.')
plot(tplot,Dc,'-b.')
plot(tplot,DOC,'-c.')
title('...in terms of carbon'),xlabel('Time'), ylabel('Carbon (mmolC*m-3)')
%legend('Phyto1','Phyto2','Zoo','DIC/100','Detritus','DOC')
grid on

subplot(2,2,4)
plot(tplot,pf_t(yr2start:end,:))
title('Grazing preference'),xlabel('Time'), ylabel('N.D.')
legend('Grazing on P1','Grazing on P2')

figure(fignum+1), clf
%subplot(2,2,1), plot(tplot,DIC/PO4)
subplot(2,2,1)
plot(tplot,P1c./(P1c+P2c))
legend('fraction eukary'), ylabel('Fraction in C')
subplot(2,2,2)
plot(tplot,P1c./P1p,'-g.'), hold on
plot(tplot,P2c./P2p,'-m.')
plot(tplot,Zc./Zp,'-k.')
plot(tplot,c2p_uptake_t(yr2start:end,:))
plot([tplot(1) tplot(end)],[ZooCP_ref ZooCP_ref],'k-')
ylabel('C:P')
legend('P1 biomass','P2 biomass','Z biomass','P1 uptake','P2 uptake','Z ref')
subplot(2,2,3)
plot(tplot,DOC./DOP,'-r.'), hold on
plot(tplot,Dc./Dp,'-b.')
plot(tplot,c2p_uptake_t(yr2start:end,:))
plot([tplot(1) tplot(end)],[RedCP RedCP],'k-')
ylabel('C:P')
legend('DOM','Detritus','P1 uptake','P2 uptake','Redfield')
subplot(2,2,4)
plot(tplot,sum(gspn_adj_t(yr2start:end,:),2),'-r.'), hold on
plot(tplot,sum(gspp_adj_t(yr2start:end,:),2),'-g.')
plot(tplot,sum(gspc_adj_t(yr2start:end,:),2),'-b.')
ylabel('Zoo homeostatic flux (mmolX*m-3/day)')
legend('N','P','C')

figure(fignum+2), clf
subplot(2,2,1)
plot(tplot,tempC0(yr2start:end)), ylabel('Temperature (C)')
subplot(2,2,2)
plot(tplot,P1c./P1n,'-g.'), hold on
plot(tplot,P2c./P2n,'-m.')
plot(tplot,Zc./Zn,'-k.')
plot(tplot,c2n_uptake_t(yr2start:end,:))
plot([tplot(1) tplot(end)],[ZooCN_ref ZooCN_ref],'k-')
ylabel('C:N')
legend('P1 biomass','P2 biomass','Z biomass','P1 uptake','P2 uptake','Z ref')
subplot(2,2,3)
plot(tplot,DOC./DON,'-r.'), hold on
plot(tplot,Dc./Dn,'-b.')
plot(tplot,c2n_uptake_t(yr2start:end,:))
plot([tplot(1) tplot(end)],[RedCN RedCN],'k-')
ylabel('C:N')
legend('DOM','Detritus','P1 uptake','P2 uptake','Redfield')
subplot(2,2,4)
plot(tplot,parz0(yr2start:end)), ylabel('PAR (W/m2)')

figure(fignum+3), clf
subplot(2,2,1)
plot(tplot,NO3/PO4), hold on
ylabel('N:P'), legend('NO3/PO4')
subplot(2,2,2)
plot(tplot,P1n./P1p,'-g.'), hold on
plot(tplot,P2n./P2p,'-m.')
plot(tplot,Zn./Zp,'-k.')
plot(tplot,n2p_uptake_t(yr2start:end,:))
plot([tplot(1) tplot(end)],[ZooCP_ref/ZooCN_ref ZooCP_ref/ZooCN_ref],'k-')
ylabel('N:P')
legend('P1 biomass','P2 biomass','Z biomass','P1 uptake','P2 uptake','Z ref')
subplot(2,2,3)
plot(tplot,DON./DOP,'-r.'), hold on
plot(tplot,Dn./Dp,'-b.')
plot(tplot,n2p_uptake_t(yr2start:end,:))
plot([tplot(1) tplot(end)],[RedCP/RedCN RedCP/RedCN],'k-')
ylabel('N:P')
legend('DOM','Detritus','P1 uptake','P2 uptake','Redfield')
subplot(2,2,4)
plot(tplot,nplimit_phyto(yr2start:end,1),'g-'), hold on
plot(tplot,nplimit_phyto(yr2start:end,2),'b.')
plot(tplot,nplimit_zoo(yr2start:end,1),'r-')
plot(tplot,nplimit_zoo(yr2start:end,2),'m.')
ylabel('N vs. P limitation')
title('Phyto limit: 1=N, 2=P; Zoo limit: 11=N, 12=P')
%if phyto is N limited; zoo becomes P limited...flips...make sense?

figure(fignum+4), clf
subplot(2,2,1), plot(tplot,gppn_t(yr2start:end,1),'-b.'), hold on
plot(tplot,epn_t(yr2start:end,1),'-g.')
plot(tplot,gspn_t(yr2start:end,1))
plot(tplot,mpn_t(yr2start:end,1))
axis([300 1100 0 1.3])
ylabel('P1 Flux terms (mmolN*m-3/day)')
legend('GPPn1','EPn1','GSPn1','MPn1')

subplot(2,2,2), plot(tplot,gppn_t(yr2start:end,2),'-b.'), hold on
plot(tplot,epn_t(yr2start:end,2),'-g.')
plot(tplot,gspn_t(yr2start:end,2))
plot(tplot,mpn_t(yr2start:end,2))
axis([300 1100 0 1.3])
ylabel('P2 Flux terms (mmolN*m-3/day)')
legend('GPPn2','EPn2','GSPn2','MPn2')

subplot(2,2,3), plot(tplot,sum(gspn_t(yr2start:end,:),2),'-b.'), hold on
plot(tplot,sum(ezn_t(yr2start:end,:),2),'-g.')
plot(tplot,gtpn_t(yr2start:end,1))
plot(tplot,mzn_t(yr2start:end,1))
plot(tplot,sum(gspn_adj_t(yr2start:end,:),2),'-r.')
axis([300 1100 0 1.3])
ylabel('Zoo Flux terms (mmolN*m-3/day)')
legend('Sum(GSPn)','Sum(EZn)','GTPn','MZn','Sum(Homeo N)')

subplot(2,2,4)
plot(tplot,gppn_t(yr2start:end,:),'-'), hold on
plot(tplot,gspn_t(yr2start:end,:),'-.')
axis([300 1100 0 1.3])
ylabel('N flux (mmolN*m-3/day)'),legend('GPPn1','GPPn2','GSPn1','GSPn2')

return
