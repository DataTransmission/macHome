clear all
addpath('~/matscript')
%
ind = 2 
fname1={'sst.dat','taux.dat','prcp.dat'};
fid1=fopen(fname1{ind},'r');
nyear = 20;
nseason = 4;
nmonth = 12;
nx = [84 87 64];
ny = [30 32 24];
cmin1 = [-3 -0.04 -6];
cmax1 = [ 3  0.04  6];
cmin2 = [-1.5 -0.02 -1.5];
cmax2 = [ 1.5  0.02  1.5];
nt = nyear*nseason;
nymingrid = ny(ind) - round(ny(ind)*2/3);
nymaxgrid = ny(ind) - round(ny(ind)*1/3);
nxmin = 120;
nxmax = 280;
ntmin = 1980;
ntmax = 2000;
ii = 1;


for j = 1:nyear
  for i = 1:nseason
  seasonal = fread(fid1,[nx(ind),ny(ind)],'real*4');
  seasonalmeridmean(:,ii)=mean(seasonal(:,nymingrid:nymaxgrid),2);
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
contourf(nxmin:(nxmax-nxmin)/(nx(ind)-1):nxmax,ntmin:nyear/nt:ntmax-nyear/nt,interannualvar',30,'linestyle','none')
colorbar
ntitle={'SST','TAUX','PRCP'};
nsource={'[CAC Jan 1980-2000]','[NOAA NCEP-NCAR CDAS-1 1980-2000]','[NOAA NCEP CPC Merged-Analysis CMAP]'};
title(['Interannual Seasonal ' ntitle{ind} ' Annomaly '  nsource{ind}]);
ylabel('year')
xlabel('longitude')
colormap(myblue2red)
caxis([cmin1(ind) cmax1(ind)])
filename={'InterannualSSTa','Interannualtauxa','Interannualprcpa'};
print(gcf,filename{ind},'-djpeg','-r200')

figure
contourf(nxmin:(nxmax-nxmin)/(nx(ind)-1):nxmax,1:nmonth/nseason:nmonth,seasonala',30,'linestyle','none')
colorbar
title(['Seasonal ' ntitle{ind} ' Anomaly ' nsource{ind}]);
ylabel('month')
xlabel('longitude')
colormap(myblue2red)
caxis([cmin2(ind) cmax2(ind)])
filename={'SeasonalSSTa','Seasonaltauxa','Seasonalprcpa'};
print(gcf,filename{ind},'-djpeg','-r200')
