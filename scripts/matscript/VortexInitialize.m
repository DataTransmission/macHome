%p-coor has p = [1000, 850, 700, 500, 300, 200, 100, 50]

%bg stands for background, tc stands for tropical cyclone

%Bckgrd zbg is in p-coor
%All TC field variables are calculated in bckgrd p-coor, by plugging in bckgrd zbg 
%we have zbg data in p-coor initially, 
%TC ztc is from TC ptc (directly from bckgrd zbg) still all in p-coor

%calc the r(lam,phi,z) in spherical coor, instead of polar coor,
%give a 3D profile of grid zbg(j,i,k) sounding
%calc the axismm variable (p(z(p)),t(z(p)),u(z(p)),v(z(p))) distr for a vortex in p-coor
%calc the p(lam,phi,z(p)), t(lam,phi,z(p)), u(lam,phi,z(p)), v(lam,phi,z(p)) 
%vortex distribution in spherical p-coor

%i,j,k represent the grids
%i,j,k is the grid point number, not the exact coordinate position
%i.e phi(i) means phi on ith grid, p(j,i,k) means pressure on grid (j,i,k) 



format long

%============READIN BCKGRD VARIABLES=======================================
%==========================================================================
%nc_dump('gT42NoOMP.cam.0301-06.nc')
addpath('/Users/ginochen/Documents/MATLAB/Research/mexcdf/mexnc')
addpath('/Users/ginochen/Documents/MATLAB/Research/mexcdf/snctools')
addpath('/Users/ginochen/Documents/MATLAB/Research')
addpath('/Users/ginochen/Documents/MATLAB/GC/proj')
%============Lo Res========================================================
%psobs = nc_varget('gT42NoOMP.cam.0301-03.nc','PS'); %ps(1,64,128)
%zbg1 = nc_varget('gT42NoOMP.cam.0301-03.nc','Z3'); %z3(z,y,x)=z3(8,64,128)
%ubg = nc_varget('gT42NoOMP.cam.0301-03.nc','U');
%vbg = nc_varget('gT42NoOMP.cam.0301-03.nc','V');
%============Hi Res========================================================
%psobs = nc_varget('gcontrol.cam.0301-06.nc','PS')./100; %ps(1,64,128)
   %this is the bckgrd ps, not the TC ps
zbg1 = nc_varget('gcontrol.cam.0301-06.nc','Z3'); %z3(z,y,x)=z3(8,64,128)
ubg = nc_varget('gcontrol.cam.0301-06.nc','U');
vbg = nc_varget('gcontrol.cam.0301-06.nc','V');
%==========================================================================



%============PARAMETERS====================================================
%==========================================================================
a=6.37122*10^6,                 %a=radius of earth, 
omg=7.292115*10^(-5),           %omg=rot speed of earth, 
rd=287.04,                      %rd=dry air gas ctc, 
g=9.80616,                      %g=gravity, 
zt=15000,                       %zt=tropopaus ht,
e=10^-25,                       %e=small ctc to avoid division by zero, 
eps=2*10^-13                    %eps=convergence limit for fixed point iterations
%----determined param----------------
q0=21*10^-3,                    %q0=sfc spec hum, 
qt=10^-8*10^-3,                 %qt=upper atm spec hum, 
zq1=3000,                       %zq1 zq2=ctc for spec hum profile, 
zq2=8000, 
t0=302.15                       %t0=sfc temp sst, 
tv0=t0*(1 + 0.608*q0),          %tv0=sfc bckgrd tv(virtual temp), 
gam=0.007,                      %gam=tv lapse rate, 
tvt=tv0-gam*zt,                 %tvt=upper atm t set to a ctc,
p0=1000,                        %p0=bckgrd sfc p, 
dp=11.15,                       %dp=p0 - ps, sfc p diff btwn bckgrd sfc p0 and ps at ini vortex center, 
rp=282*10^3,                    %rp & zp=ctc for p fit, 
zp=7000,
%------specify genesis domain---------
phic=20;                        %phic(lamc)=center lat(lon) of initial vortex, 
lamc=300; 
fc=2*omg*sind(phic);            %fc=coriolis param, 
jj=65:90;                       %genesis zone lat grid point, not exact longitude
ii=200:220;                     %genesis zone lon grid point
kk=1:8;
nn=1:10;                        %iterate times for ztc
%==========================================================================




%==============GRID POINT==================================================
%==========================================================================
pmod = [1000, 850, 700, 500, 300, 200, 100, 50]
lam1(1) = 0; %take index 58,72
% for i = 1:127
for i = 1:255 %i grids for hi res
%   lam1(i+1) = lam1(i) + 2.8125; %low res
    lam1(i+1) = lam1(i) + 1.40625;%high res
end
%============Low Res=======================================================
% phi1 = [...
%     -87.8638, -85.0965, -82.3129, -79.5256, -76.7369,...
%     -73.9475, -71.1578, -68.3678, -65.5776, -62.7874,...
%     -59.997,  -57.2066, -54.4162, -51.6257, -48.8352,...
%     -46.0447, -43.2542, -40.4636, -37.6731, -34.8825,...
%     -32.0919, -29.3014, -26.5108, -23.7202, -20.9296,... 
%     -18.139,  -15.3484, -12.5578, -9.76715, -6.97653,...
%     -4.18592, -1.39531, 1.39531,  4.18592,  6.97653,...
%     9.76715,  12.5578,  15.3484,  18.139,   20.9296,...
%     23.7202,  26.5108,  29.3014,  32.0919,  34.8825,...
%     37.6731,  40.4636,  43.2542,  46.0447,  48.8352,...
%     51.6257,  54.4162,  57.2066,  59.997,   62.7874,...
%     65.5776,  68.3678,  71.1578,  73.9475,  76.7369,...
%     79.5256,  82.3129,  85.0965,  87.8638];
%============Hi Res========================================================
phi1 = [... %128
    -88.9277353522959, -87.5387052130273, -86.1414721015279, -84.7423855907142,  -83.3425960440704,... 
    -81.9424662991732, -80.5421464346171, -79.1417096486217, -77.7411958655139,  -76.3406287023715,...
    -74.9400230196494, -73.5393886337675, -72.1387322891624, -70.7380587725176,  -69.3373715749609,...
    -67.9366733025785, -66.5359659401756, -65.1352510260352, -63.7345297708429,  -62.3338031405324,...
    -60.9330719152074, -59.5323367318266, -58.1315981156439, -56.7308565037137,  -55.3301122627028,...
    -53.9293657025561, -52.5286170870997, -51.1278666423533, -49.7271145631097,  -48.3263610181882,...
    -46.9256061546646, -45.5248501013023, -44.1240929713558, -42.723334864877,   -41.3225758706231,...
    -39.9218160676465, -38.5210555266244, -37.1202943109789, -35.719532477824,   -34.3187700787707,...
    -32.918007160614,  -31.5172437659226, -30.1164799335463, -28.7157156990552,  -27.3149510951204,...
    -25.9141861518467, -24.5134208970629, -23.1126553565776, -21.7118895544042,  -20.3111235129604,...
    -18.9103572532454, -17.5095907949986, -16.1088241568413, -14.7080573564048,  -13.3072904104462,...
    -11.9065233349538, -10.5057561452436, -9.10498885604852, -7.70422148160049,  -6.3034540357076,...
    -4.90268653182654, -3.5019189831313,  -2.10115140257898, -0.700383802973324, 0.700383802973324,...
    2.10115140257898,  3.5019189831313,   4.90268653182654,  6.3034540357076,    7.70422148160049,...
    9.10498885604852,  10.5057561452436,  11.9065233349538,  13.3072904104462,   14.7080573564048,...
    16.1088241568413,  17.5095907949986,  18.9103572532454,  20.3111235129604,   21.7118895544042,...
    23.1126553565776,  24.5134208970629,  25.9141861518467,  27.3149510951204,   28.7157156990552,...
    30.1164799335463,  31.5172437659226,  32.918007160614,   34.3187700787707,   35.719532477824,...
    37.1202943109789,  38.5210555266244,  39.9218160676465,  41.3225758706231,   42.723334864877,...
    44.1240929713558,  45.5248501013023,  46.9256061546646,  48.3263610181882,   49.7271145631097,...
    51.1278666423533,  52.5286170870997,  53.9293657025561,  55.3301122627028,   56.7308565037137,...
    58.1315981156439,  59.5323367318266,  60.9330719152074,  62.3338031405324,   63.7345297708429,...
    65.1352510260352,  66.5359659401756,  67.9366733025785,  69.3373715749609,   70.7380587725176,...
    72.1387322891624,  73.5393886337675,  74.9400230196494,  76.3406287023715,   77.7411958655139,...
    79.1417096486217,  80.5421464346171,  81.9424662991732,  83.3425960440704,   84.7423855907142,...
    86.1414721015279,  87.5387052130273,  88.9277353522959]
%==========================================================================



%============Genesis Zone==================================================
%==========================================================================
%============Lo Res========================================================
% for k=1:8
%     for j=33:40
%         for i=114:128
%         zbg(i-113,j-32,k)=zbg1(k,j,i);
%         end
%     end
% end
%============Hi Res========================================================
for k=kk
    for j=jj
        for i=ii
        zbg(j-jj(1)+1,i-ii(1)+1,k)=zbg1(k,j,i); 
        %restrict z to the area of TC initialize zone
        end
    end
end
%============Lo Res========================================================
%lam(1:15) = lam1(114:128);
%phi(1:8) = phi1(33:40);
%============Hi Res========================================================
lam(1:length(ii)) = lam1(ii); %restrict lam, phi to the area of TC initialize zone
phi(1:length(jj)) = phi1(jj);
%==========================================================================



%================CALCULATE VAR=============================================
%==========================================================================
for n=1
   for i=1:length(ii) %1:15
      for j=1:length(jj) %1:8
         r  = a*acos(sind(phic)*sind(phi(j)) + cosd(phic)*cosd(phi(j))*cosd(lam(i)-lamc));
         d1 = sind(phic)*cosd(phi(j)) - cosd(phic)*sind(phi(j))*cosd(lam(i)-lamc);
         d2 = cosd(phic)*sind(lam(i)-lamc);
         d  = max(e, sqrt(d1^2 + d2^2));
         ps = p0 - dp*exp(-(r/rp)^1.5); %ps is the sfc-p for initial TC
         pt = p0*(tvt/tv0)^(g/(rd*gam)); %pt=tropop zt p
         for k=kk
            if zbg(j,i,k) < zt %check whether it's ztc < zt instead
              %********Pressure****************
              ptc(j,i,k) = (p0 - dp*exp(-(r/rp)^1.5)*exp(-(zbg(j,i,k)/zp)^2))*((tv0 - ...
                           gam*zbg(j,i,k))/tv0)^(g/(rd*gam)); 
              %p(phi,lam,z) distribution for vortex ztc use, since zbg
              %cannot directly transfer to ztc, need to transfer to p
              %
              %********GeopotentialHeight****** ps >= p >= pt (below tropopause and above sfc p)
              ztc(j,i,k) = tv0/gam*(1 - (ptc(j,i,k)/ps)^(rd*gam/g)); 
              for n=nn
                f = pmod(k) - ptc(j,i,k);
                dpdz = 2*dp*zbg(j,i,k)/zp^2*exp(-(r/rp)^1.5)*exp(-(zbg(j,i,k)/zp)^2)*...
                      ((tv0 - gam*zbg(j,i,k))/tv0)^(g/(rd*gam)) - g/(rd*tv0)*...
                      (p0 - dp*exp(-(r/rp)^1.5)*exp(-(zbg(j,i,k)/zp)^2))*...
                      ((tv0 - gam*zbg(j,i,k))/tv0)^(g/(rd*gam) - 1);
                %derivative of ptc
                ztc(j,i,k) = ztc(j,i,k) + f./dpdz; 
                %this is the ztc interatively calc at each model level,
                %each level should interate 10 times, 
                %until ztc(10) approx equal ztc(9), variation 
                %use ztc to calc qb, vt and t instead of using zbg
              end
              %********SpecificHumidity********
              qb = q0*exp((-ztc(j,i,k)/zq1))*exp(-(ztc(j,i,k)/zq2)^2); %q(z) dist for vortex
              %vt(phi,lam,z) dist. for vortex
              vt(j,i,k) = -fc*r/2 + sqrt(fc^2*r^2/4 - (1.5*(r/rp)^1.5*(tv0-gam*ztc(j,i,k))*rd) / ...
                          (1+2*rd*(tv0-gam*ztc(j,i,k))*ztc(j,i,k)/(g*zp^2) - p0/dp*exp((r/rp)^1.5)*...
                          exp((ztc(j,i,k)/zp)^2))); 
              %********Temperature*************
              t(j,i,k) = (tv0-gam*ztc(j,i,k))/(1+0.608*qb)*(1+(2*rd*(tv0-gam*ztc(j,i,k))*ztc(j,i,k)) / ...
                         (g*zp^2*(1-p0/dp*exp((r/rp)^1.5)*exp((ztc(j,i,k)/zp)^2))))^-1; 
              %t(phi,lam,z) dist for vortex and bckgrd, t approach bckgrd temp as r becomes large  
            else
              %********Pressure****************
              ptc(j,i,k) = pt*exp(-((g*(zbg(j,i,k)-zt))/(rd*tvt)));
              %********Geopotential Height****** pt > p (above tropopause)
              ztc(j,i,k) = zt + rd*tvt/g*log(pt/ptc(j,i,k));
              %********Temperature*************
              t(j,i,k) = tvt;
              %********SpecificHumidity********             
              qb = qt;
              %********Tangential Velocity*****
              vt(j,i,k) = 0;
            end
            u(j,i,k) = vt(j,i,k)*d1/d;
            v(j,i,k) = vt(j,i,k)*d2/d;
         end
      end
   end
end



%============SAVE VAR======================================================
%==========================================================================
fid = fopen('ztc_T85.dta','w');
for k=kk
    fprintf(fid,'%8.2f',ztc(:,:,k)); %save ztc into ztc_T85.dta file
    fprintf(fid,'\n');
end
fclose(fid);
%==========================================================================



%============PLOT VAR======================================================
%==========================================================================
%============Lo Res========================================================
% plot(1:15,vt(:,4,1)')
% contour(ptc(:,:,1)')
% quiver(squeeze(ubg(1,:,:)),squeeze(vbg(1,:,:))); hold on;
% quiver(114:128,33:40,u(:,:,1)',v(:,:,1)')
% xlim([0 128]);ylim([0 64]);
%============Hi Res========================================================
%--------3D UV Vector plot for TC-------------------------------
figure
[x,y,z]=meshgrid(lam1(ii),phi1(jj),kk);
w=zeros([length(jj),length(ii),length(kk)]);
scale=2;
quiver3(x,y,z,u,v,w,scale);hold on;
text(x(1),y(1),z(1),'u = 10 m/s');
uscl(1,1,1)=10;
uscl(2:length(jj),2:length(ii),2:length(kk))=zeros;
quiver3(x,y,z,uscl,w,w,scale);
daspect([1,1,1]); %adjust to correct aspect ratio
%--------2D UV Vector plot for TC on Global field---------------
figure
k=1
quiver(0:360/(length(lam1)-1):360,-90:180/(length(phi1)-1):90,squeeze(ubg(k,:,:)),squeeze(vbg(k,:,:))); hold on;
quiver(x(:,:,k),y(:,:,k),u(:,:,k),v(:,:,k))
xlim([0 360]);ylim([-90 90]);
%--------2D T plot for TC on Global field---------------
figure
for k=1:6
if k<=2
    subplot(3,2,k);
    contourf(x(:,:,k),y(:,:,k),t(:,:,k),1000,'linestyle','none');colorbar;caxis([299,max(max(t(:,:,k)))]);
else 
    subplot(3,2,k);
    contourf(x(:,:,k),y(:,:,k),t(:,:,k),20,'linestyle','none');colorbar;caxis([min(min(t(:,:,k))),max(max(t(:,:,k)))]);
end
end
%==========================================================================
