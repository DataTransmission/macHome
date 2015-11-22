ncfile = 'MJO_5120x2560x32_4km_10s_QOBS_EQX_nopert_1280_0000135000_U.nc';
mode = 'NC_NOWRITE';
varname = {'x','y','z','time','p','U'}
nx = numel(x);
ny = numel(y);
nz = numel(z);
ndx = 27; % ndx of 4km boxes adds up to ndx*4 = 108km approximately one degree
%function (ncfile,mode,varname,edge)
% read in CRM data and coarse-grain to SCAM resolution
% x, y, z, time, p, U (V,W)
% U(time, z, y, x)
nv = numel(varname);

for iv = 1:nv 
% get var value
   eval(sprintf('%s = ncread(ncfile,varname{iv});',varname{iv})); 
end
% find the spaced-index corresponding to the coars grid of 1 degree
%
% grid spacing for each level of CRM (dxj,dyj)_k
for ix = 1:ndx-1:nx
   for iy =
      Ub() = U();
   end
end
% dx
% dy
% dz
