format long
%---------------------------------------
%matlab_script to calc domain parameters
%---------------------------------------
zfilename='grd';
nz=30;
norder=3; % fitting polynomial to create the vertical spacing from 'grd' data
dolinearAbove = 1;
maxdz = 80;

dx=0.08 ;% km
dy=0.08 ;% km
domain_x=120 ;% SCAM box = 312 km (64x128 spectral grid resolution), otherwise typical climate predition usese 100km
domain_y=120 ;% 
subdomain_x=1.2;%
subdomain_y=1.2;%
[nsubdomain_x,nsubdomain_y,nx_gl,ny_gl] = domainparm(domain_x,domain_y,subdomain_x,subdomain_y,dx,dy)
%rem(domain_x,dx)==0 %& rem(domain_y,dy)==0
%rem(domain_x,subdomain_x)==0 & rem(domain_y,subdomain_y)==0


% Calculate vertical grid 'grd' by fitting data from 'grd.original'
[zfit dz] = grdprofile(zfilename,nz,norder,dolinearAbove,maxdz);
if (dolinearAbove == 0)
   maxdz = max(dz);
end
nz_gl = length(zfit)
