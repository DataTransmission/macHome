
clear all
addpath('~/matscript')
ind = 1
fname1={'bogusvort'};
fid1=fopen(fname1{ind});
nx = [48];
ny = [40];
ndays = 5
cmin = 
for t = 1:ndays
    bogus(:,:,t) = fread(fid1,[nx(ind),ny(ind)],'real*4')
    contourf(bogus(:,:,t),20,'linstyle','none')
    caxis([cmin cmax])
    print(gcf,['bogusvort' num2str(t)],'-djpeg','-r200')
    pause
end
