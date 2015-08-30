
clear

nmonth = 11;
indx = 1
nvar = 4;
nx = 180;
ny = 91;
longdim  = 1;
latdim   = 2;
monthdim = 3;
dlat = 2*pi/180;
lat1 = -90:180/90:90;
lat = -0.5*pi:dlat:0.5*pi;
dx1 = 2*110.574*10^3; %m
dx = cos(lat)*dx1;
varsum = 0;
area = dx*dx1;


for ind = 1:nvar
  ind
  fname1={'solr.dat','lwfx.dat','lhfx.dat','shfx.dat'};
  fid1=fopen(fname1{ind},'r');
  for i = 1:nmonth
      var(:,:,i,ind) = fread(fid1,[nx,ny],'real*4');
%contourf(var(:,:,i,1));pause
  end
  for i = 1:nx
    for j = 1:ny
      for k = 1:nmonth
        if var(i,j,k,ind)==-999
          var(i,j,k,ind)=NaN;
        end
      end
    end
  end
  
  for j = 1:nx
    for i = 1:nmonth
      var1(j,:,i,ind) = var(j,:,i,ind);       %W/m/m
    end
  end 
  varsum = var1(:,:,:,ind) + varsum;
end

totalheatflux_xavg1 = squeeze(mean(nanmean(varsum(:,:,:),longdim),monthdim));
%totalheatflux_xavg2 = totalheatflux_xavg1; 
totalheatflux_xavg2 = squeeze(nansum(mean(varsum(:,:,:),monthdim),longdim)); 
% use monthly mean to kick out the missing data at each point on earth
% we will have many missing data at higher latitudes 
ii=1;
NHeatTransport(1) = totalheatflux_xavg2(1);
for i=1:length(totalheatflux_xavg2)-1
  NHeatTransport(i+1) = nansum(totalheatflux_xavg2(ii) + NHeatTransport(i));
  ii=ii+1;
end
%plot the seasonal cycle of heat flux pick Jan, Aug
ntitle={'Jan Surface Heat Flux Climatology','August Surface Heat Flux Climatology',...
        'Annual Surface Heat Flux Climatology','Northward Ocean Heat Flux' ''};
nsource={' [OBERHUBER: Max Planck 2x2 Global] '}
nfile ={'JanHF.jpg','AugHF.jpg','AnnualHF.jpg','NHF.jpg'};
figure
plot(lat1,nanmean(varsum(:,:,1),longdim)); %without x and var dep, 
refline(0,0)
xlabel('latitude')
ylabel('W/m/m')
xlim([-90 90])
title([ntitle{indx} nsource{indx}])
print(nfile{1},'-djpeg','-r150')


figure
plot(lat1,nanmean(varsum(:,:,8),longdim)); %without x and var dep
refline(0,0)
xlabel('latitude')
ylabel('W/m/m')
xlim([-90 90])
title([ntitle{2} nsource{indx}])
print(nfile{2},'-djpeg','-r150')


figure
plot(lat1,totalheatflux_xavg1); %without x and var dep
%don't use nanmean in monthly mean, b/c many month have missing value at certain latitudes
refline(0,0)
xlabel('latitude')
ylabel('W/m/m')
xlim([-90 90])
title([ntitle{3} nsource{indx}])
print(nfile{3},'-djpeg','-r150')
%plot the annual mean heat flux 

figure
plot(lat1,NHeatTransport.*area); %without x and var dep
%plot(NHeatTransport); %without x and var dep
refline(0,0)
xlabel('latitude')
ylabel('W')
xlim([-90 90])
title([ntitle{4} nsource{indx}],'position',[5 2.56433e+15 1.00011])
print(nfile{4},'-djpeg','-r150')
%text(-100, -96, nsource, 'clipping', 'off');
