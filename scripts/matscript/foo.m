clear all
%
%fname1=['sst.dat']
fname1=['taux.dat'];
%fname1=['prcp.dat'];
%fname2=['ssta.dat'];
%fname2=['tauxa.dat'];
%fname2=['prcpa.dat'];
%fname2=['thrmcla.dat'];
fid1=fopen(fname1,'r');
ind = 3
nyear = 20;
nseason = 4;
nx = [84 84 80 322]; 
ny = [30 32 24 20];
nt = nyear*nseason

ii = 1;
for j = 1:nyear
  for i = 1:nseason
  seasonalmeridmean(:,ii)=mean(fread(fid1,[nx(ind),ny(ind)],'real*4'),2);
  ii = ii+1;
  end
end
 
ii = 1;
for i = 1:nseason
  seasonalclim = sum(seasonalmeridmean(:,i:nseason:nt-nseason+ii),2)/nyear; 
  for j = 0:nyear-1
    interannualvar(:,ii+j*nseason) = seasonalmeridmean(:,ii+j*nseason) - seasonalclim;
  end
  ii = ii+1;
  ss(:,i) = seasonalclim;
end
  yearlyclim = sum(ss,2)/nseason;
for i = 1:nseason
  seasonala(:,i) = ss(:,i) - yearlyclim;
end
figure
contourf(120:(280-120)/nx(ind):280-160/nx(ind),1980:20/nt:2000-20/nt,interannualvar',20,'linestyle','none')
colorbar
title('Interannual Seasonal taux Annomaly')
ylabel('month')
xlabel('longitude')
%print(gcf,'InterannualSSTa','-djpeg','-r200')
%print(gcf,'Interannualtauxa','-djpeg','-r200')

figure
contourf(120:(280-120)/nx(ind):280-160/nx(ind),1:12/nseason:12,seasonala',20,'linestyle','none')
colorbar
title('Seasonal taux Anomaly');
ylabel('month')
xlabel('longitude')
%print(gcf,'SeasonSSTa','-djpeg','-r200')
%print(gcf,'Seasonaltauxa','-djpeg','-r200')
