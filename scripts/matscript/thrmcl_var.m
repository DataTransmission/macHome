clear all
fname1=['thrmcl.dat'];
fid1=fopen(fname1,'r');
ind = 1
nyear = 20;
nseason = 4;
nmonth = 12;
nx = [320]; 
ny = [20]; %thrm is actually nz
nt = nyear*nseason;
cmin1 = -90;
cmax1 =  90;
cmin2 = -15;
cmax2 =  15;
isotherm = 20;
nxmin = 120;
nxmax = 280;
ntmin = 1980;
ntmax = 2000; 
nz = [5 15 25 35 46 57 70 82 96 112 129 148 171 197 229 268 317 381 465 579];
ii = 1;
for j = 1:nyear
  for i = 1:nseason
  [lon,dep]=find(fread(fid1,[nx(ind),ny(ind)],'real*4') - isotherm < 1); 
  %row is lon col is depth, this will compare all temp with 20 degree, 
  %see if it's close enough to 20 for the defined depth of thermocline
    for k = 1:nx(ind)
        seasonalthrmcldepth(k,ii) = nz(dep(find(lon==k,1))); 
  %first depth reaching in each k lon
    end
  ii= ii+1;
  end
end
 
ii = 1;
for i = 1:nseason
  seasonalclim = sum(seasonalthrmcldepth(:,i:nseason:nt-nseason+ii),2)/nyear; 
  for j = 0:nyear-1
    interannualvar(:,ii+j*nseason) = seasonalthrmcldepth(:,ii+j*nseason) - seasonalclim;
  end
  ii = ii+1;
  ss(:,i) = seasonalclim;
end
  yearlyclim = sum(ss,2)/nseason;
for i = 1:nseason
  seasonala(:,i) = ss(:,i) - yearlyclim;
end
figure
contourf(nxmin:(nxmax-nxmin)/(nx(ind)-1):nxmax,ntmin:nyear/nt:ntmax-nyear/nt,interannualvar',20,'linestyle','none')
colorbar
ntitle={'thermocline depth (m)'}
nsource={'[CARTON-GIESE SODA Version 1.2]'}
title(['Interannual Seasonal ' ntitle{ind} ' Anomaly ' nsource{ind}])
ylabel('year')
xlabel('longitude')
colormap(blue2red)
caxis([cmin1(ind) cmax1(ind)])
filename={'Interannualthrmcla'};
print(gcf,filename{ind},'-djpeg','-r200')

figure
contourf(nxmin:(nxmax-nxmin)/(nx(ind)-1):nxmax,1:nmonth/nseason:nmonth,seasonala',20,'linestyle','none')
colorbar
title(['Seasonal ' ntitle{ind} ' Anomaly ' nsource{ind}]);
ylabel('month')
xlabel('longitude')
colormap(blue2red)
caxis([cmin2(ind) cmax2(ind)])
filename={'Seasonalthrmcla'};
print(gcf,filename{ind},'-djpeg','-r200')
