clear
format long
%==============parameter================================
lat   = 30*pi/180;
angfreq   = 1/(24*3600);
f0    = 2*angfreq*sin(lat);
g     = 10;
rho   = 1000;
drho  = 5 
gp    = g*drho/rho;
D     = 300;
lamb  = sqrt(gp*D)/f0;                     %300km
re    = 6300*10^3;

cnum = input('(1) short wave, (2) stationery (3) long wave ?    ')

switch cnum
  case 1
    nx  = 2*pi*re/10;                      %zoom in
    xgrid = lamb/100;                      %smaller grid to resolve small wavelength
    ndays = 300;
    nt  = ndays*24*3600;                     %10 month
    wavelength = lamb;                     %wavelength < lamb
  case 2
    nx  = 2*pi*re;
    xgrid = lamb/2;
    ndays = 500;
    nt  = ndays*24*3600;                     %15 month
    wavelength = 2*pi*lamb;                %stationary wavelength k = 2*pi/wavelength = 1/lamb 
  case 3
    nx  = 2*pi*re;
    xgrid = lamb;
    ndays = 500;
    nt  = ndays*24*3600;                     %15 month
    wavelength = 3*pi*lamb;                %wavelength > 2*pi*lamb
end
%if nx/500<= lamb
x     = 0:xgrid:nx;
tgrid = nt/500;
t     = 1:tgrid:nt;                        %time step cannot be too big 
beta  = 2*angfreq*cos(lat)/re;
k      = 2*pi/wavelength;
dk     = 10^-1*k;
k1     = k - dk;
k2     = k + dk;
omg1   = -beta * k1 / (k1^2 + lamb^-2);
omg2   = -beta * k2 / (k2^2 + lamb^-2);
phiamp1 = 1;
phiamp2 = 1;

%==============variables================================
for i = 1:length(x)
  for j = 1:length(t)
    phi1 = phiamp1 * cos(k1*x(i) - omg1*t(j));
    phi2 = phiamp2 * cos(k2*x(i) - omg2*t(j));
    phi(j,i) = phi1 + phi2;
  end
end
contourf(0:360*xgrid/(2*pi*re):360*(nx+1)/(2*pi*re),1:ndays*tgrid/nt:ndays+0.5,phi,20,'linestyle','none');
ylabel('ndays')
xlabel('latitude')
