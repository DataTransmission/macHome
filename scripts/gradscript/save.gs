***************************************************************************************
*       $Id: save.gs,v 1.30 2008/07/02 21:36:36 bguan Exp $
*       Copyright (C) 2006 Bin Guan.
*       Distributed under GNU/GPL.
***************************************************************************************
function save(arg)
*
* Save variables to .dat and .ctl pair.
*
rc=gsfallow('on')

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
_.file.1='saveout'
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
* Write .dat file.
*
'set fwrite '_.path.1'/'_.file.1'.dat'
'set gxout fwrite'
tcnt=_ts
while(tcnt<=_te)
'set t 'tcnt
vcnt=1
while(vcnt<=num_var)
zcnt=_zs
while(zcnt<=_ze)
'set z 'zcnt
'display const('_.var.vcnt','_.undef.1',-u)'
zcnt=zcnt+1
endwhile
vcnt=vcnt+1
endwhile
say 'T='tcnt': written.'
tcnt=tcnt+1
endwhile
'disable fwrite'
'set gxout contour'

*
* Write .ctl file.
*
lines=8
line.1='DSET ^'_.file.1'.dat'
line.2='UNDEF '_.undef.1
line.3=_xdef
line.4=_ydef
line.5=_zdef
line.6=_tdef
line.7='VARS 'num_var
line.8='ENDVARS'
cnt=1
while(cnt<=lines-1)
status=write(_.path.1'/'_.file.1'.ctl',line.cnt)
cnt=cnt+1
endwhile
cnt=1
while(cnt<=num_var)
varline=_.name.cnt' '_nz0' 99 '_.name.cnt
status=write(_.path.1'/'_.file.1'.ctl',varline)
cnt=cnt+1
endwhile
status=write(_.path.1'/'_.file.1'.ctl',line.lines)
status=close(_.path.1'/'_.file.1'.ctl')

*
* Restore original dimension environment.
*
'set x '_xs_old' '_xe_old
'set y '_ys_old' '_ye_old
'set z '_zs' '_ze
'set t '_ts' '_te

return
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
say '  Save variables to .dat and .ctl pair.'
say ''
say '  USAGE: save -v <var1> [<var2>...] [-n <name1> [<name2>...]] [-u <undef>] [-o <file>] [-p <path>]'
say '     <var>: variable to be saved.'
say '     <name>: name for a variable in saved .ctl file. Original name will be used if unset.'
say '     <undef>: undef value for .dat and .ctl. Defaults to the value found in ctlinfo.'
say '     <file>: common name for .dat and .ctl pair. Defaults to saveout.'
say '     <path>: path to saved files. Do NOT include trailing "/". Defaults to current path.'
say ''
say '  Dependencies: parsestr.gsf, parseopt.gsf, qdims.gsf'
say ''
say '  Copyright (C) 2006 Bin Guan.'
say '  Distributed under GNU/GPL.'
return
