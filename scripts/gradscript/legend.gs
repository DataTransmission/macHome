***************************************************************************************
*       $Id: legend.gs,v 1.12 2010/03/28 05:07:34 bguan Exp $
*       Copyright (C) 2009 Bin Guan.
*       Distributed under GNU/GPL.
***************************************************************************************
function legend(arg)
*
* Display the legend for a line graph (the latter must have been produced using plot.gs).
*
rc=gsfallow('on')

if(arg='-h')
usage()
return
endif

tmpdir='/tmp'
whoamifile='./.whoami.bGASL'
'!whoami>'whoamifile
whoami=sublin(read(whoamifile),2)
rc=close(whoamifile)
'!unlink 'whoamifile
mytmpdir=tmpdir'/bGASL-'whoami
'!mkdir -p 'mytmpdir

*
* Read legend data.
*
flag=1
cnt=0
while(flag)
result=read(mytmpdir'/plot.txt~')
status=sublin(result,1)
if(status!=0)
flag=0
else
cnt=cnt+1
line.cnt=sublin(result,2)
endif
endwhile
num_var=cnt

*
* Initialize other options.
*
x=0
y=0
_.orientation.1='v'
_.position.1='tl'
_.xoffset.1=0
_.yoffset.1=0
_.scalefactor.1=1

*
* Parse -orient option (orientation).
*
rc=parseopt(arg,'-','orient','orientation')

*
* Parse -xo option (x offset).
*
rc=parseopt(arg,'-','xo','xoffset')

*
* Parse -yo option (y offset).
*
rc=parseopt(arg,'-','yo','yoffset')

*
* Parse -scale option (scaling factor).
*
rc=parseopt(arg,'-','scale','scalefactor')

*
* Get plot area
*
'query gxinfo'
line3=sublin(result,3)
line4=sublin(result,4)
x1=subwrd(line3,4)
x2=subwrd(line3,6)
y1=subwrd(line4,4)
y2=subwrd(line4,6)

*
* Define sizes and spacing.
*
mark_size=0.11
line_length=0.55
small_spacing=0.11
'query pp2xy 0 0'
tmpxa=subwrd(result,3)
'query pp2xy 1 1'
tmpxb=subwrd(result,3)
rvratio=tmpxb-tmpxa
mark_size=mark_size*rvratio*_.scalefactor.1
line_length=line_length*math_sqrt(_.scalefactor.1)
small_spacing=small_spacing*rvratio*_.scalefactor.1

*
* Draw legend.
*
if(_.position.1='tl')
xx=x1+(0.5*mark_size+small_spacing)
yy=y2-(0.5*mark_size+small_spacing)
endif
xx=xx+_.xoffset.1
yy=yy+_.yoffset.1
x=xx
y=yy
cnt=1
while(cnt<=num_var)
mark=subwrd(line.cnt,1)
style=subwrd(line.cnt,2)
color=subwrd(line.cnt,3)
thick=subwrd(line.cnt,4)
wrdcnt=5
text=subwrd(line.cnt,wrdcnt)
while(subwrd(line.cnt,wrdcnt)!='')
wrdcnt=wrdcnt+1
text=text' 'subwrd(line.cnt,wrdcnt)
endwhile
'query string 'text
string_width=subwrd(result,4)
if(_.orientation.1='h'&x+line_length+small_spacing+string_width>x2)
x=xx
y=y-(mark_size+small_spacing)
endif
'set line 'color' 'style' 'thick
'draw mark 'mark' 'x' 'y' 'mark_size
if(style=0)
'draw mark 'mark' 'x+0.5*line_length' 'y' 'mark_size
endif
* Note: above: if no line then draw a third mark in the middle
if(style!=0)
if(mark=2|mark=4|mark=10|mark=11)
'draw line 'x+0.5*mark_size' 'y' 'x+line_length-0.5*mark_size' 'y
endif
if(mark=7)
'draw line 'x+0.34*mark_size' 'y' 'x+line_length-0.34*mark_size' 'y
endif
if(mark=8)
'draw line 'x+0.27*mark_size' 'y' 'x+line_length-0.27*mark_size' 'y
endif
if(mark!=2&mark!=4&mark!=10&mark!=11&mark!=7&mark!=8)
'draw line 'x' 'y' 'x+line_length' 'y
endif
endif
'draw mark 'mark' 'x+line_length' 'y' 'mark_size
'set string 1 l'
'draw string 'x+line_length+small_spacing' 'y' 'text
if(_.orientation.1='v')
y=y-(mark_size+small_spacing)
endif
if(_.orientation.1='h')
x=x+(line_length+small_spacing+string_width+2*small_spacing)
endif
cnt=cnt+1
endwhile
'set string 1 bl'

return
***************************************************************************************
function usage()
*
* Print usage information.
*
say '  Display the legend for a line graph (the latter must have been produced using plot.gs).'
say ''
say '  Usage: legend [-orient v|h] [-xo <xoffset>] [-yo <yoffset>] [-scale <scalingfactor>]'
say '     -orient h: use for horizontal oriention. Default is vertical.'
say '     <xoffset>: Horizontal offset to default position (i.e., top-left). Defaults to 0.'
say '     <yoffset>: Vertical offset to default position (i.e., top-left). Defaults to 0.'
say '     <scalingfactor>: scaling factor for mark size. Will NOT affect string size (Use "set strsiz" for that). Defaults to 1.'
say ''
say '  "legend -h" prints this usage information.'
say ''
say '  Dependencies: parsestr.gsf, parseopt.gsf'
say ''
say '  See also: plot.gs'
say ''
say '  Copyright (C) 2009 Bin Guan.'
say '  Distributed under GNU/GPL.'
return
