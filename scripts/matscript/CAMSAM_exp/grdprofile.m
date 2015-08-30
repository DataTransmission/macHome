function [zfit dz] = grdprofile(zfilename, nz, norder, dolinearAbove, maxdz)

% obtain the vertical grid from 'grd'
ztmp = readgrd(zfilename); % 0 is add to the surface layer
z = [0,ztmp'];
%z = [0,ztmp(ztmp<=izlinearAbove)'];
% nz : desired number of vetical level

% make a linearly spaced coordinate
x=[1:length(z)];

% least-square fit a norder standard polynomial through the paired points (x,z), z=span(1,x,x^2,...)
coeff=polyfit(x,z,norder);

% evaluate the desired points of z(x), where x=linspace(1,length(z),nz) 
xfit = linspace(1,length(z),nz);
zfit=polyval(coeff,xfit)';

% vertical spacing
dz = zfit(2:end)-zfit(1:end-1);

% find the index that reaches maxdz and set all level interval above that equaling maxdz
if (dolinearAbove)
   i = min(find(dz>=maxdz))-1; % minimum index that has dz value exceeding maxdz

   while (zfit(i-1)<=max(z))
      zfit(i) = zfit(i-1) + maxdz;
      i=i+1;
   end
   dz = zfit(2:end)-zfit(1:end-1);
end

