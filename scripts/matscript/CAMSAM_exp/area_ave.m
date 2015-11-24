%ncfile = 'MJO_5120x2560x32_4km_10s_QOBS_EQX_nopert_1280_0000135000_U.nc';
ncfile_CAM = '/projects/rsmas/kirtman/gchen/cesm/run/cam_only_aquaplanet/cam_only_aquaplanet.cam.rh0.0001-01-01-10800.nc' 
lev = ncread(ncfile_CAM,'lev')
lat = ncread(ncfile_CAM,'lat')
lon = ncread(ncfile_CAM,'lon')
varname = 'U' % V, W
ncfile = [varname '_tav.nc'];
dimname = {'x','y','z','time','p'}
for id = 1:numel(dimname) 
   eval(sprintf('%s = ncread(ncfile,dimname{id});',dimname{id})); 
end
nx = numel(x); % nx = 5120; x_crm = 20480000m; earth circumference = 39690000m
ny = numel(y); % ny = 2560; y_crm = 10240000m;
nz = numel(z);
nt = 1; % number of time samples 2592000/10800 = 240 = 30 days
nbox_x = 32; % nbox_x of 4km boxes adds up to nbox_x*4 = 128km approximately 1.16 degrees
nbox_y = 32; 
var = ncread(ncfile,varname); % time-summed var
% var(time, z, y, x)
var = squeeze(var(1,:,:,:))/nt; % average of the time-summed var
%
% CRM xy average (xy bar)
% Purpose : treat turbulence to be locally deviated from local box
% climatology
for iz = 1:nz
   ii = 1; 
   for ix = 1:nbox_x:nx
      jj = 1;
      for iy = 1:nbox_y:ny
         varxyb(ii,jj,iz) = mean(mean(var(iz,iy:iy+nbox-1,ix:ix+nbox-1),2),3); 
         jj = jj + 1;
      end
      ii = ii + 1;
   end
end
%
% CRM y average (y bar)
% Purpose : treat turbulence to be locally deviated from latitudinal 
% climatology (zonally symmetric)
varyb_tmp = squeeze(mean(var,3)); % remove zonal dependence
for iz = 1:nz
   jj = 1;
   for iy = 1:nbox:ny
      varyb(jj,iz) = mean(varyb_tmp(iz,iy:iy+nbox-1),2);
      jj = jj + 1;
   end
end
%
% TKE = xymean(u'^2 + v'^2 + w'^2)/2
% TPERT : Perturbation temperature (eddies in PBL)
% PBLH 
%

