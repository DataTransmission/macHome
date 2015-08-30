***************************************************************************************
*       $Id: fcorr.gs,v 1.6 2010/05/26 22:50:28 bguan Exp bguan $
*       Copyright (C) 2007 Bin Guan.
*       Distributed under GNU/GPL.
***************************************************************************************
function fcorr(arg)
*
* Point-by-point temporal correlations between two fields.
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

input1=subwrd(arg,1)
input2=subwrd(arg,2)
if(input2='')
usage()
return
endif
output=subwrd(arg,3)
if(output='')
output='fcorrout'
endif

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

*
* Calculate correlations and write .dat file.
*
*'set t '_ts' '_te
*'fcorrtmp1='input1
*'fcorrtmp2='input2
'set fwrite 'mytmpdir'/fcorr.dat~'
'set gxout fwrite'
'set t '_ts
zcnt=_zs
while(zcnt<=_ze)
'set z 'zcnt
ycnt=_ys
while(ycnt<=_ye)
'set y 'ycnt
xcnt=_xs
while(xcnt<=_xe)
'set x 'xcnt
*'fcorrtmp3=tcorr(fcorrtmp1,fcorrtmp2,t='_ts',t='_te')'
'fcorrtmp3=tcorr('input1','input2',t='_ts',t='_te')'
'display const(fcorrtmp3,'default_undef()',-u)'
xcnt=xcnt+1
endwhile
ycnt=ycnt+1
endwhile
zcnt=zcnt+1
endwhile
'disable fwrite'
'set gxout contour'

*
* Write .ctl file.
*
lines=9
line.1='DSET ^fcorr.dat~'
line.2='UNDEF 'default_undef()
line.3=_xdef
line.4=_ydef
line.5=_zdef
line.6=_tdef
line.7='VARS 1'
line.8=output' '_nz0' 99 Add description here.'
line.9='ENDVARS'
cnt=1
while(cnt<=lines)
status=write(mytmpdir'/fcorr.ctl~',line.cnt)
cnt=cnt+1
endwhile
status=close(mytmpdir'/fcorr.ctl~',line.cnt)

*
* Open .ctl file and fetch values.
*
'open 'mytmpdir'/fcorr.ctl~'
file_num=file_number()
'set x '_xs_old' '_xe_old
* The above line is needed to ensure that there will not be a gap near the prime meridian in global maps if unintended.
'set y '_ys_old' '_ye_old
'set z '_zs' '_ze
output'='output'.'file_num
'close 'file_num
*'undefine fcorrtmp1'
*'undefine fcorrtmp2'
'undefine fcorrtmp3'

*
* Restore original dimension environment.
*
*'set x '_xs_old' '_xe_old
*'set y '_ys_old' '_ye_old
*'set z '_zs' '_ze
*These were already set above.
'set t '_ts' '_te

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
function default_undef()
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
say '  Point-by-point temporal correlations between two fields.'
say ''
say '  Usage: fcorr <input1> <input2> [<output>]'
say '     <input1>: Input field 1.'
say '     <input2>: Input field 2.'
say '     <output>: Output field. Defaults to "fcorrout".'
say ''
say '  Dependencies: qdims.gsf'
say ''
say '  Copyright (C) 2007 Bin Guan.'
say '  Distributed under GNU/GPL.'
return
