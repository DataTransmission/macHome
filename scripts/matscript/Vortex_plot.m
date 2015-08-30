clear
plev = load('plev.dat');
utc  = load('utc.dat');
vtc  = load('vtc.dat');
ztc  = load('ztc.dat');
ttc  = load('ttc.dat');
vortc = load('vortc.dat');
divtc = load('divtc.dat');
i = 1;
tcxgrid = 5;
tcygrid = 2;
while i<=length(utc)-2*tcygrid
  contourf(vortc(i:i+2*tcygrid,:))
'vortc'
pause
  contourf(divtc(i:i+2*tcygrid,:))
'divtc'
pause
  contourf(ttc(i:i+2*tcygrid,:))
'ttc'
pause
  quiver(utc(i:i+2*tcygrid,:),vtc(i:i+2*tcygrid,:))
'utc'
pause
  contourf(ztc(i:i+2*tcygrid,:))
'ztc'
pause
  i = i+2*tcygrid + 1
end
