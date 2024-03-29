***************************************************************************************
*	$Id: ssn2day.gs,v 1.6 2008/07/02 22:20:58 bguan Exp $
*	Copyright (C) 2008 Bin Guan.
*	Distributed under GNU/GPL.
***************************************************************************************
function ssn2day(arg)
*
* Seasonal-to-daily regridding.
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

*
* Parse -v option (variables to be saved).
*
num_var=parseopt(arg,'-','v','var')
if(num_var=0)
usage()
return
endif

*
* Initialize other options.
*
cnt=1
while(cnt<=num_var)
_.name.cnt=_.var.cnt
cnt=cnt+1
endwhile
_.undef.1=default_undef()
_.file.1=''
_.path.1='.'

*
* Parse -n option (name to be used in .ctl).
*
rc=parseopt(arg,'-','n','name')

*
* Parse -u option (undef value to be used in .dat and .ctl).
*
rc=parseopt(arg,'-','u','undef')

*
* Parse -o option (common name for .dat and .ctl pair).
*
rc=parseopt(arg,'-','o','file')

*
* Parse -p option (path to output files).
*
rc=parseopt(arg,'-','p','path')

'set gxout fwrite'
'set fwrite 'mytmpdir'/ssn2dayout.dat~'

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

_tims_old=_tims
_time_old=_time

'set t '_ts' '_te
vcnt=1
while(vcnt<=num_var)
'sn2dytmp'vcnt'='_.var.vcnt
vcnt=vcnt+1
endwhile

'set time JAN'_yrs' OCT'_yre
*Note: Seasonal means were assumed to have been assigned to Jan, Apr, Jul, Oct.
*That is, four seasons were defined as DJF, MAM, JJA, SON, and seasonal means
*were assigned to central months.
qdims()
seasons=(_yre-_yrs+1)*4
days=(_yre-_yrs+1)*365
writectl(mytmpdir'/ssn2dayout.ctl~','^ssn2dayout.dat~',days,num_var,name)

'set time JAN'_yrs

season=1
while(season<=seasons+1)
offset=season-1
flag=math_mod(season,4);
if(season=1); drepeat=59; endif
if(season>1 & season<seasons+1 & flag=1); drepeat=90; endif
if(season>1 & season<seasons+1 & flag=2); drepeat=92; endif
if(season>1 & season<seasons+1 & flag=3); drepeat=92; endif
if(season>1 & season<seasons+1 & flag=0); drepeat=91; endif
if(season=seasons+1); drepeat=31; endif
dcnt=1
while(dcnt<=drepeat)
vcnt=1
while(vcnt<=num_var)
zcnt=_zs
while(zcnt<=_ze)
'set z 'zcnt
'daytmp=sn2dytmp'vcnt'(t+'offset')'
'display const(daytmp,'_.undef.1',-u)'
zcnt=zcnt+1
endwhile
vcnt=vcnt+1
endwhile
dcnt=dcnt+1
endwhile
season=season+1
endwhile

'disable fwrite'
vcnt=1
while(vcnt<=num_var)
'undefine sn2dytmp'vcnt
vcnt=vcnt+1
endwhile
'undefine daytmp'
'set gxout contour'

if(_.file.1!='')
writectl(_.path.1'/'_.file.1'.ctl','^'_.file.1'.dat',days,num_var,name)
'!cp 'mytmpdir'/ssn2dayout.dat~ '_.path.1'/'_.file.1'.dat'
endif

'set x '_xs_old' '_xe_old
'set y '_ys_old' '_ye_old
'set z '_zs' '_ze
'set time '_tims_old' '_time_old

dfile_old=dfile()
'open 'mytmpdir'/ssn2dayout.ctl~'
file_num=file_number()
'set dfile 'file_num
vcnt=1
while(vcnt<=num_var)
_.name.vcnt'='_.name.vcnt'.'file_num
vcnt=vcnt+1
endwhile
'set dfile 'dfile_old

return
***************************************************************************************
function writectl(ctlfile,datfile,nt,nv,var)
*
* Write the .ctl file for the temporary .dat file
*
lines=9
line.1='DSET 'datfile
line.2='UNDEF '_.undef.1
line.3='OPTIONS 365_DAY_CALENDAR'
line.4=_xdef
line.5=_ydef
line.6=_zdef
* Note: 'nt' below is an argument of function 'writectl', not the global variable '_nt'.
line.7='TDEF 'nt' LINEAR '_tims' 1dy'
line.8='VARS 'nv
line.9='ENDVARS'
cnt=1
while(cnt<=lines-1)
status=write(ctlfile,line.cnt)
cnt=cnt+1
endwhile
cnt=1
while(cnt<=nv)
varline=_.var.cnt' '_nz0' 99 '_.var.cnt
status=write(ctlfile,varline)
cnt=cnt+1
endwhile
status=write(ctlfile,line.lines)
status=close(ctlfile)

return
***************************************************************************************
function dfile()
*
* Get the default file number.
*
'q file'

line1=sublin(result,1)
dfile=subwrd(line1,2)

return dfile
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
say '  Seasonal-to-daily regridding.'
say ''
say '  Usage: ssn2day -v <var1> [<var2>...] [-n <name1> [<name2>...]] [-u <undef>] [-o <file>] [-p <path>]'
say '     <var>: seasonal variable.'
say '     <name>: name for daily variable. Same as <var> if unset.'
say '     <undef>: undef value for .dat and .ctl. Defaults to the value found in ctlinfo.'
say '     <file>: common name for output .dat and .ctl pair. No file output if unset.'
say '     <path>: path to output files. Do NOT include trailing "/". Defaults to current path.'
say ''
say '  Note: Input data MUST be 3-month means centered on Jan, Apr, Jul and Oct, and the first (last) time step MUST be Jan (Oct).'
say '        Output data will start from 01Jan and ends at 31Dec (the last Dec always has missing values).'
say ''
say '  Dependencies: parsestr.gsf, parseopt.gsf, qdims.gsf'
say ''
say '  Copyright (C) 2008 Bin Guan.'
say '  Distributed under GNU/GPL.'
return
