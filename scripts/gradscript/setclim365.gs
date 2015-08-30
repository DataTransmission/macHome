***************************************************************************************
*       $Id: setclim365.gs,v 1.2 2008/07/02 21:36:42 bguan Exp bguan $
*       Copyright (C) 2008 Bin Guan.
*       Distributed under GNU/GPL.
***************************************************************************************
function setclim365(arg)
*
* Set a variable to be climatology.
*
rc=gsfallow('on')

tmpdir='/tmp'
whoamifile='./.whoami.bGASL'
'!whoami>'whoamifile
whoami=sublin(read(whoamifile),2)
rc=close(whoamifile)
'!unlink 'whoamifile
mytmpdir=tmpdir'/bGASL-'whoami
'!mkdir -p 'mytmpdir

input=subwrd(arg,1)
output=subwrd(arg,4)
clim_tims=subwrd(arg,2)
clim_time=subwrd(arg,3)
if(clim_time='')
usage()
return
endif
if(output='')
output=input
endif

qdims()

_ts_old=_ts
_te_old=_te

*
* get number of t-point(s) per year (e.g., 12 for monthly data, 4 for seasonal data...).
*
'set time JAN2000 JAN2001'
'query dims'
line5=sublin(result,5)
tJAN1=subwrd(line5,11);tJAN1=math_int(tJAN1)
tJAN2=subwrd(line5,13);tJAN2=math_int(tJAN2)
NtperYR=tJAN2-tJAN1
say 'SETCLIM365>INFO: 'NtperYR' time grids per calendar year.'

*
* Ensure x-coordinates are integers and there are no redundant grid points.
*
qdims()
_xs_old=_xs
_xe_old=_xe
if(math_int(_xs)!=_xs | math_int(_xe)!=_xe)
xs_new=math_nint(_xs)
xe_new=math_nint(_xe)
'set x 'xs_new' 'xe_new
qdims()
endif
if(_lone-_lons>=360)
rddnt_points=(_lone-_lons-360)/_dlon+1
'set x '_xs' '_xe-rddnt_points
qdims()
endif

*
* Ensure y-coordinates are integers.
*
qdims()
_ys_old=_ys
_ye_old=_ye
if(math_int(_ys)!=_ys | math_int(_ye)!=_ye)
ys_new=math_nint(_ys)
ye_new=math_nint(_ye)
'set y 'ys_new' 'ye_new
qdims()
endif

'set t '_ts_old' '_te_old
'setclimtmp='input

clim_short_ts=_ts_old

'set time 'clim_tims' 'clim_time
qdims()
clim_ts=_ts
clim_te=_te
if(clim_te-clim_ts+1<NtperYR)
say 'SETCLIM365>ERROR: Time period must be >= one year.'
return
endif

offset=math_mod(clim_ts-clim_short_ts,NtperYR)
if(offset<0);offset=offset+NtperYR;endif
cnt_head=NtperYR-offset
cnt_foot=math_mod(clim_te-NtperYR-clim_short_ts,NtperYR)+1
if(cnt_foot<0);cnt_foot=cnt_foot+NtperYR;endif
cnt_body=(clim_te-clim_ts+1)-cnt_head-cnt_foot

'set t 'clim_short_ts

*
* Write .dat file.
*
'set fwrite 'mytmpdir'/setclim.dat~'
'set gxout fwrite'
*
* Head
*
cnt=offset
while(cnt<=cnt_head+offset-1)
zcnt=_zs
while(zcnt<=_ze)
'set z 'zcnt
'display const(setclimtmp(t+'cnt'),'undef()',-u)'
zcnt=zcnt+1
endwhile
cnt=cnt+1
endwhile
*
* Body
*
cnt=1
cnt2=0
while(cnt<=cnt_body)
zcnt=_zs
while(zcnt<=_ze)
'set z 'zcnt
'display const(setclimtmp(t+'cnt2'),'undef()',-u)'
zcnt=zcnt+1
endwhile
cnt=cnt+1
cnt2=cnt2+1
if(cnt2>=NtperYR);cnt2=0;endif
endwhile
*
* Foot
*
cnt=0
while(cnt<=cnt_foot-1)
zcnt=_zs
while(zcnt<=_ze)
'set z 'zcnt
'display const(setclimtmp(t+'cnt'),'undef()',-u)'
zcnt=zcnt+1
endwhile
cnt=cnt+1
endwhile
'disable fwrite'
'set gxout contour'

*
* Write .ctl file.
*
lines=10
line.1='DSET ^setclim.dat~'
line.2='UNDEF 'undef()
line.3='OPTIONS 365_DAY_CALENDAR'
line.4=_xdef
line.5=_ydef
line.6=_zdef
line.7='TDEF 'clim_te-clim_ts+1' LINEAR '_tims' '_dtim
line.8='VARS 1'
line.9=output' '_nz0' 99 'output
line.10='ENDVARS'
cnt=1
while(cnt<=lines)
status=write(mytmpdir'/setclim.ctl~',line.cnt)
cnt=cnt+1
endwhile
status=close(mytmpdir'/setclim.ctl~')

*
* Restore original dimension environment.
*
'set x '_xs_old' '_xe_old
'set y '_ys_old' '_ye_old
'set z '_zs' '_ze
'set time 'clim_tims' 'clim_time
'open 'mytmpdir'/setclim.ctl~'
output'='output'.'file_number()
'close 'file_number()
'set t '_ts_old' '_te_old

return
***************************************************************************************
function file_number()
*
* Get the number of files opened.
*
'q files'
if(result='No Files Open')
return 0
endif

lines=1
while(sublin(result,lines+1)!='')
lines=lines+1
endwhile

return lines/3
***************************************************************************************
function undef()
*
* Get undef value from the default .ctl file.
*
'q ctlinfo'
if(result='No Files Open')
return 'unknown'
endif

lines=1
while(1)
lin=sublin(result,lines)
if(subwrd(lin,1)='undef'|subwrd(lin,1)='UNDEF')
return subwrd(lin,2)
endif
lines=lines+1
endwhile
***************************************************************************************
function usage()
*
* Print usage information.
*
say '  Dependencies: qdims.gsf'
say ''
say '  Copyright (C) 2008 Bin Guan.'
say '  Distributed under GNU/GPL.'
return
