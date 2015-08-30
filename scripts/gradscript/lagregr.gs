***************************************************************************************
*	$Id: lagregr.gs,v 1.11 2008/07/02 01:58:16 bguan Exp bguan $
*	Copyright (C) 2008 Bin Guan.
*	Distributed under GNU/GPL.
***************************************************************************************
function lagregr(arg)
*
* Lag regressions between two variables.
*
rc=gsfallow('on')

switch='regr'

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
output=subwrd(arg,3)
_lags=subwrd(arg,4)
_lage=subwrd(arg,5)
if(input2='')
usage()
return
endif
if(output='')
output='lag'switch'out'
endif
if(_lags='')
_lags=0
endif
if(_lage='')
_lage=0
endif
lags=_lage-_lags+1

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

if(_ts=_te)
say 'LAGCORR>ERROR: Time not varying.'
say 'Use "set t ..." or "set time ..." to set a varying time dimension.'
return
endif

writectl(mytmpdir,lags,output,switch)

'set t '_ts' '_te
'lgcrtmp1='input1
'lgcrtmp2='input2

'set gxout fwrite'
'set fwrite 'mytmpdir'/lag'switch'.dat~'

lag=_lags
while(lag<=_lage)
zcnt=_zs
while(zcnt<=_ze)
'set z 'zcnt
'set t '_ts' '_te
'lgcrtmp2lagged=tloop(lgcrtmp2(t+'lag'))'
'set t 1'
'lgcrtmp3=t'switch'(lgcrtmp1,lgcrtmp2lagged,t='_ts',t='_te')'
'display const(lgcrtmp3,'default_undef()',-u)'
zcnt=zcnt+1
endwhile
lag=lag+1
endwhile

'disable fwrite'
'undefine lgcrtmp1'
'undefine lgcrtmp2'
'undefine lgcrtmp2lagged'
'undefine lgcrtmp3'
'set gxout contour'

'open 'mytmpdir'/lag'switch'.ctl~'
file_num=file_number()
'set x '_xs_old' '_xe_old
* The above line is needed to ensure that there will not be a gap near the prime meridian in global maps if unintended.
'set y '_ys_old' '_ye_old
'set z '_zs' '_ze
'set t '_lags' '_lage
*offset=_ts-_lags
offset=1-_lags
output'=tloop('output'.'file_num'(t+'offset'))'
'close 'file_num

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
function default_tims()
*
* Get the beginning time step of the default file.
*
'set t 1'
'query dims'
lt=sublin(result,5)
tims=subwrd(lt,6)

return tims
***************************************************************************************
function writectl(mytmpdir,lags,output,switch)
*
* Write the .ctl file for the temporary .dat file
*
lines=8
line.1='DSET ^lag'switch'.dat~'
line.2='UNDEF 'default_undef()
line.3=_xdef
line.4=_ydef
line.5=_zdef
*line.6='TDEF 'lags' LINEAR '_tims' '_dtim
line.6='TDEF 'lags' LINEAR 'default_tims()' '_dtim
line.7='VARS 1'
line.8='ENDVARS'
cnt=1
while(cnt<=lines-1)
status=write(mytmpdir'/lag'switch'.ctl~',line.cnt)
cnt=cnt+1
endwhile
cnt=1
while(cnt<=1)
varline=output' '_nz0' 99 Add description here.'
status=write(mytmpdir'/lag'switch'.ctl~',varline)
cnt=cnt+1
endwhile
status=write(mytmpdir'/lag'switch'.ctl~',line.lines)
status=close(mytmpdir'/lag'switch'.ctl~')

return
***************************************************************************************
function usage()
*
* Print usage information.
*
say '  Lag regressions between two variables.'
say ''
say '  Usage: lagregr <input1> <input2> [<output> [<lag_start> [<lag_end>]]]'
say '     <input1>: Independent variable. Can have a vertical dimension. NO horizontal dimensions.'
say '     <input2>: Dependent variable. Can have a vertical dimension and/or up to two horizontal dimensions.'
say '     <output>: Output variable. Defaults to "lagregrout".'
say '               E.g., "set t 0" and then <output>(t-3) will be lag regression when <input2> leads <input1> by 3 time steps.'
say '     <lag_start>: Beginning lag. E.g., <lag1>=-3 for <input2> leading <input1> by 3 time steps. Defaults to 0.'
say '     <lag_end>: Ending lag. E.g., <lag2>=3 for <input2> lagging <input2> by 3 time steps. Defaults to 0.'
say ''
say '  Note: Do NOT use <output>(t=number), instead, "set t 0" and use <output>(t+number) to get the value for Lag(number).'
say '        To use/plot all lags at once, simply "set t <lag1> <lag2>".'
say ''
say '  Dependencies: qdims.gsf'
say ''
say '  Copyright (C) 2008 Bin Guan.'
say '  Distributed under GNU/GPL.'
return
