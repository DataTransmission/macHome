***************************************************************************************
*       $Id: plot.gs,v 1.6 2010/03/28 05:38:51 bguan Exp bguan $
*       Copyright (C) 2008 Bin Guan.
*       Distributed under GNU/GPL.
***************************************************************************************
function plot(arg)
*
* Draw a line graph.
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
* Parse -v option (variable).
*
num_var=parseopt(arg,'-','v','variable')

if(num_var=0)
usage()
return
endif

*
* Parse -r option (range).
*
range_rc=parseopt(arg,'-','r','range')

*
* Initialize other options.
*
cnt=1
while(cnt<=num_var)
_.mark.cnt=cnt+1
_.style.cnt=-1
_.color.cnt=cnt
_.thick.cnt=-1
_.text.cnt='Variable 'cnt
cnt=cnt+1
endwhile

*
* Parse -m option (mark).
*
rc=parseopt(arg,'-','m','mark')

*
* Parse -s option (style).
*
rc=parseopt(arg,'-','s','style')

*
* Parse -c option (color).
*
rc=parseopt(arg,'-','c','color')

*
* Parse -k option (thick).
*
rc=parseopt(arg,'-','k','thick')

*
* Parse -t option (text).
*
rc=parseopt(arg,'-','t','text')

*
* get lower/upper boundaries
*
if(range_rc<=1)
qdims()
cnt=1
while(cnt<=num_var)
say 'PLOT>INFO: getting range...'
'varmintmp'cnt'=min(min(min(min('_.variable.cnt',t='_ts',t='_te'),x='_xs',x='_xe'),y='_ys',y='_ye'),z='_zs',z='_ze')'
'varmaxtmp'cnt'=max(max(max(max('_.variable.cnt',t='_ts',t='_te'),x='_xs',x='_xe'),y='_ys',y='_ye'),z='_zs',z='_ze')'
* In the above two lines time dimension is done first since that's much quicker.
'query defval varmintmp'cnt' 1 1'
varmin.cnt=subwrd(result,3)
'query defval varmaxtmp'cnt' 1 1'
varmax.cnt=subwrd(result,3)
cnt=cnt+1
endwhile
lowerbnd=varmin.1
upperbnd=varmax.1
cnt=2
while(cnt<=num_var)
if(lowerbnd>varmin.cnt)
lowerbnd=varmin.cnt
endif
if(upperbnd<varmax.cnt)
upperbnd=varmax.cnt
endif
cnt=cnt+1
endwhile
else
lowerbnd=_.range.1
upperbnd=_.range.2
endif

'set vrange 'lowerbnd' 'upperbnd
say 'PLOT>INFO: range = 'lowerbnd' to 'upperbnd'.'

cnt=1
while(cnt<=num_var)
if(_.mark.cnt!=-1);'set cmark '_.mark.cnt;endif
if(_.style.cnt!=-1);'set cstyle '_.style.cnt;endif
if(_.color.cnt!=-1);'set ccolor '_.color.cnt;endif
if(_.thick.cnt!=-1);'set cthick '_.thick.cnt;endif
'display '_.variable.cnt
cnt=cnt+1
endwhile

cnt=1
while(cnt<=num_var)
line=_.mark.cnt' '_.style.cnt' '_.color.cnt' '_.thick.cnt' '_.text.cnt
rc=write(mytmpdir'/plot.txt~',line)
cnt=cnt+1
endwhile
rc=close(mytmpdir'/plot.txt~')

return
***************************************************************************************
function usage()
*
* Print usage information.
*
say '  Draw a line graph.'
say ''
say '  Usage: plot -v <var1> [<var2>...] [-r <range_from> <range_to>] [-m <mark1> [<mark2>...]] [-s <style1> [<style2>...]] [-c <color1> [<color2>...]] [-k <thick1> [<thick2>...]] [-t <text1> [<text2>...]]'
say '     <var>: variable to be plotted.'
say '     <range_from>, <range_to>: set the axis limit. Defaults to the minimum and maximum values.'
say '     <mark>: Defaults to "2 3 4...", i.e., open circle, closed circle, open square, closed square, and so on.'
say '     <style>: Defaults to current GrADS setting.'
say '     <color>: Defaults to "1 2 3...", i.e., foreground color, red, green, dark blue, and so on.'
say '     <thick>: Defaults to current GrADS setting.'
say '     <text>: Text to be shown in the legend (use "legend.gs"). Text beginning with a minus sign or containing spaces must be double quoted.'
say ''
say '  Dependencies: parsestr.gsf, parseopt.gsf, qdims.gsf'
say ''
say '  See also: legend.gs'
say ''
say '  Copyright (C) 2008 Bin Guan.'
say '  Distributed under GNU/GPL.'
return
