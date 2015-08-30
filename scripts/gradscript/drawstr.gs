***************************************************************************************
*       $Id: drawstr.gs,v 1.31 2009/02/23 17:45:04 bguan Exp bguan $
*       Copyright (C) 2009 Bin Guan.
*       Distributed under GNU/GPL.
***************************************************************************************
function drawstr(arg)
*
* Draw string in designated position.
*
rc=gsfallow('on')

*
* Parse -T option (main title).
*
num_TXT=parseopt(arg,'-','T','TEXT')

*
* Parse -t option (text).
*
num_txt=parseopt(arg,'-','t','text')
if(num_TXT+num_txt=0 | num_TXT*num_txt>0)
usage()
return
endif

*
* Initialize other options.
*
if(num_TXT>0)
cnt=1
while(cnt<=num_TXT)
_.size.cnt=0.18
_.thickness.cnt=5
_.xoffset.cnt=0
_.yoffset.cnt=0
cnt=cnt+1
endwhile
endif
if(num_txt>0)
x=0
y=0
cnt=1
while(cnt<=num_txt)
_.position.cnt=cnt
_.thickness.cnt=5
_.xoffset.cnt=0
_.yoffset.cnt=0
cnt=cnt+1
endwhile
endif

*
* Parse -p option (position).
*
ps=12
p.1='tl'
p.2='tc'
p.3='tr'
p.4='bl'
p.5='br'
p.6='bc'
p.7='b25'
p.8='b75'
p.9='l'
p.10='r'
p.11='tl2'
p.12='tr2'
rc=parseopt(arg,'-','p','position')
cnt=1
while(cnt<=num_txt)
if(valnum(_.position.cnt)=0)
p_cnt=1
flag=0
while(p_cnt<=ps)
flag=flag | (_.position.cnt=p.p_cnt)
p_cnt=p_cnt+1
endwhile
  if(!flag)
  say 'DRAWSTR>ERROR: invalid <position>.'
  say ''
  return
  endif
endif
if(valnum(_.position.cnt)=2)
say 'DRAWSTR>ERROR: <position> must be an integer.'
say ''
return
endif
cnt=cnt+1
endwhile

*
* Parse -s option (size).
*
sizerc=parseopt(arg,'-','s','size')

*
* Parse -k option (thickness).
*
rc=parseopt(arg,'-','k','thickness')

*
* Parse -xo option (x offset).
*
rc=parseopt(arg,'-','xo','xoffset')

*
* Parse -yo option (y offset).
*
rc=parseopt(arg,'-','yo','yoffset')

*
* Draw main title.
*
if(num_TXT>0)
'set vpage off'
'set parea off'
cnt=1
while(cnt<=num_TXT)
while(_.TEXT.cnt='')
cnt=cnt+1
endwhile
if(cnt>num_TXT)
return
endif
'query defval vpagexa'_.position.cnt' 1 1'
vpagexa=subwrd(result,3)
'query defval vpagexb'_.position.cnt' 1 1'
vpagexb=subwrd(result,3)
'query defval vpageya'_.position.cnt' 1 1'
vpageya=subwrd(result,3)
TXT_x=vpagexa+(vpagexb-vpagexa)/2
TXT_y=vpageya-0.05
* Note: vpage top is 0.22 higher than parea top given default y padding; therefore, default spacing between title and parea is 0.22-0.05=0.17.
'set string 1 bc '_.thickness.cnt
'set strsiz '_.size.cnt
'draw string 'TXT_x+_.xoffset.cnt' 'TXT_y+_.yoffset.cnt' '_.TEXT.cnt
cnt=cnt+1
endwhile
return
endif

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
x25=x1+(x2-x1)/4
x50=x1+(x2-x1)/2
x75=x1+(x2-x1)/4*3
y50=y1+(y2-y1)/2

*
* Define spacing.
*
small_spacing=0.05
big_spacing=0.5
'query pp2xy 0 0'
tmpxa=subwrd(result,3)
'query pp2xy 1 1'
tmpxb=subwrd(result,3)
rvratio=tmpxb-tmpxa
small_spacing=small_spacing*rvratio
big_spacing=big_spacing*rvratio

*
* Draw string.
*
cnt=1
while(cnt<=num_txt)
while(_.text.cnt='')
cnt=cnt+1
endwhile
if(cnt>num_txt)
return
endif
if(_.position.cnt=1 | _.position.cnt=p.1)
x=x1
y=y2+small_spacing
'set string 1 bl '_.thickness.cnt' 0'
endif
if(_.position.cnt=2 | _.position.cnt=p.2)
x=x50
y=y2+small_spacing
'set string 1 bc '_.thickness.cnt' 0'
endif
if(_.position.cnt=3 | _.position.cnt=p.3)
x=x2
y=y2+small_spacing
'set string 1 br '_.thickness.cnt' 0'
endif
if(_.position.cnt=4 | _.position.cnt=p.4)
x=x1
y=y1+small_spacing
'set string 1 bl '_.thickness.cnt' 0'
endif
if(_.position.cnt=5 | _.position.cnt=p.5)
x=x2
y=y1+small_spacing
'set string 1 br '_.thickness.cnt' 0'
endif
if(_.position.cnt=6 | _.position.cnt=p.6)
x=x50
y=y1-big_spacing
'set string 1 tc '_.thickness.cnt' 0'
endif
if(_.position.cnt=7 | _.position.cnt=p.7)
x=x25
y=y1-1.7*big_spacing
'set string 1 tc '_.thickness.cnt' 0'
endif
if(_.position.cnt=8 | _.position.cnt=p.8)
x=x75
y=y1-1.7*big_spacing
'set string 1 tc '_.thickness.cnt' 0'
endif
if(_.position.cnt=9 | _.position.cnt=p.9)
x=x1-2*big_spacing
y=y50
'set string 1 c '_.thickness.cnt' 90'
endif
if(_.position.cnt=10 | _.position.cnt=p.10)
x=x2+big_spacing
y=y50
'set string 1 c '_.thickness.cnt' 270'
endif
if(_.position.cnt=11 | _.position.cnt=p.11)
x=x1
y=y2-small_spacing
'set string 1 tl '_.thickness.cnt' 0'
endif
if(_.position.cnt=12 | _.position.cnt=p.12)
x=x2
y=y2-small_spacing
'set string 1 tr '_.thickness.cnt' 0'
endif
x=x+_.xoffset.cnt
y=y+_.yoffset.cnt
if(sizerc=num_txt)
* Need the above if to prevent possible errors since _.size.cnt has no initial value defined in this script.
'set strsiz '_.size.cnt
endif
'draw string 'x' 'y' '_.text.cnt
cnt=cnt+1
endwhile
'set string 1 bl 5 0'
return

return
***************************************************************************************
function usage()
*
* Print usage information.
*
say '  Annotate a plot.'
say ''
say '  Usage 1: drawstr -t <text1> [<text2>...] [-p <position1> [<position2>...]] [-s <size1> [<size2>...]] [-k <thickness1> [<thickness2>...]] [-xo <xoffset1> [<xoffset2>...]] [-yo <yoffset1> [<yoffset2>...]]'
say '  Usage 2: drawstr -T <TEXT1> [<TEXT2>...] [-p <position1> [<position2>...]] [-s <size1> [<size2>...]] [-k <thickness1> [<thickness2>...]] [-xo <xoffset1> [<xoffset2>...]] [-yo <yoffset1> [<yoffset2>...]]'
say '     <text>|<TEXT>: Panel labels | main titles. Text beginning with a minus sign or containing spaces must be double quoted.'
say '     <position>: Position of <text>|<TEXT>. For <text>, refer to schematic below. For <TEXT>, use the panel index to specify the location. Defaults to "1 2 3...".'
say '     <size>: size of <text>|<TEXT>. Defaults to 0.18 for <TEXT>.'
say '     <thickness>: Thickness of <text>|<TEXT>.'
say '     <xoffset>: Horizontal offset to default position. Defaults to 0.'
say '     <yoffset>: Vertical offset to default position. Defaults to 0.'
say ''
say '                <TEXT>'
say '   1               2               3'
say '   ---------------------------------'
say '   |11                           12|'
say '   |                               |'
say '   |                               |'
say '  9|             plot              |10'
say '   |                               |'
say '   |                               |'
say '   |4                             5|'
say '   ---------------------------------'
say '                   6                '
say '           7               8        '
say ''
say '  Note: The "-T" and "-t" options cannot be used together.'
say ''
say '  Dependencies: parsestr.gsf, parseopt.gsf'
say ''
say '  Copyright (C) 2009 Bin Guan.'
say '  Distributed under GNU/GPL.'
return
