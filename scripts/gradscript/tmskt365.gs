***************************************************************************************
*	$Id: tmskt365.gs,v 1.4 2008/07/03 17:14:50 bguan Exp bguan $
*	Copyright (C) 2004 Bin Guan.
*	Distributed under GNU/GPL.
***************************************************************************************
function tmskt365(arg)
*
* Mask out certain calendar days (or pentads, depending on the time resolution) of a time series.
*
rc=gsfallow('on')

input=subwrd(arg,1)
mon1=subwrd(arg,2)
mon2=subwrd(arg,3)
output=subwrd(arg,4)
if(mon2='')
usage()
return
endif
if(output='')
output=input
endif

if(valnum(mon1)!=1 | mon1<1 | mon1>365)
say 'MONMSKT ERROR: <month1> must be an integer between 1 and 365.'
say ''
return
endif
if(valnum(mon2)!=1 | mon2<1 | mon2>365)
say 'MONMSKT ERROR: <month2> must be an integer between 1 and 365.'
say ''
return
endif

qdims()

'set t '_ts' '_te
'mnmttmp='input

mon1minus1=mon1-1
'set time 01Jan 31Dec'
* In the above, must use calendar time and start from 1st calendar month/season since <start> and <end> are specified relative to 1st calendar month/season.
'mmskII=1'
'mmskII1=const(mmskII(t-'mon1minus1'),0,-u)'
'mmskII2=const(mmskII(t-'mon2'),0,-u)'
'mmskII=mmskII1-mmskII2'
if(mon1<=mon2)
'mmskII=maskout(mmskII,mmskII-0.5)'
else
'mmskII=maskout(mmskII+1,mmskII+0.5)'
endif

'setclim365 mmskII '_tims' '_time
'set t '_ts' '_te
if(output='display'|output='DISPLAY')
'display mnmttmp*mmskII'
else
output'=mnmttmp*mmskII'
endif
'undefine mnmttmp'
'undefine mmskII'
'undefine mmskII1'
'undefine mmskII2'

return
***************************************************************************************
function usage()
*
* Print usage information.
*
say '  Mask out certain calendar days (or pentads, depending on the time resolution) of a time series.'
say ''
say '  Usage 1 (for daily time series): tmskt <input> <start> <end> [<output>]'
say '     <input>: input daily time series.'
say '     <start>: the first Julian day not to mask out, e.g., 32 for Feb 1.'
say '     <end>: the last Julian day not to mask out, e.g., 59 for Feb 28.'
say '     <output>: output daily time series. Defaults to <input>.'
say '               If <output> = DISPLAY, only a graph will be displayed.'
say ''
say '  Example: tmskt sst 32 59'
say '           This would set all days but Feb 1-28 to missing values for the variable sst.'
say ''
say '  Usage 2 (for pentad time series): tmskt <input> <start> <end> [<output>]'
say '     <input>: input pentad time series.'
say '     <start>: the first Julian pentad not to mask out.'
say '     <end>: the last Julian pentad not to mask out.'
say '     <output>: output pentad time series. Defaults to <input>.'
say '               If <output> = DISPLAY, only a graph will be displayed.'
say ''
say '  Dependencies: setclim365.gs, qdims.gsf'
say ''
say '  Copyright (C) 2008 Bin Guan.'
say '  Distributed under GNU/GPL.'
return
