*The purpose of analyzing the 2 groups of data is to find 
*whether the coupled model run has different predictability
*than the prescribed model run. whether it were better or worse,
*see if the ensemble spread/variance were different or not. 
reinit
'open obsst_201008_plev.ctl'
'open sys3_201008_plev.ctl'
'set t 2 61'
*'set gxout fwrite'
*'fwrite test.dat'
nr = 1
nlon = 10
nlat = 10
nz = 1
zlev = 3
*================Variance of Ensemble at a point==========
'subplot 2 1 1'
while (np <= zlev)
  while (nr <= 2)
    'set z' nz
    'set lon ' nlon
    'set lat ' nlat
    'define zeave = ave(z.' nr ',e=1,e=41)'
    'd tloop(sqrt(ave(pow(zeave - z.' nr ',2),e=2,e=41)))'
    nr= nr + 1
  endwhile
  nz= nz + 1
endwhile
*================Variance of Ensemble avg'd over globe=====
'subplot 2 2 1'
nr = 1
while (nr <= 2)
  'define zeaave = aave(ave(z.' nr ',e=1,e=41),lon=1,lon=360,lat=-90,lat=90)'
  'define zaave = aave(z.' nr ',lon=1,lon=360,lat=-90,lat=90)'
  'd tloop(sqrt(ave(pow(zeaave - zaave,2),e=2,e=41)))'
  nr= nr + 1
endwhile


