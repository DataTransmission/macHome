***************************************************************************************
*       $Id: drawmark.gs,v 1.12 2010/03/28 05:39:02 bguan Exp $
*       Copyright (C) 2009 Bin Guan.
*       Distributed under GNU/GPL.
***************************************************************************************
function drawmark(arg)
*
* Draw marks at grid points with mark size (i.e., area) proportional to the magnitude of data.
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

var=subwrd(arg,1)
mark=subwrd(arg,2)
color=subwrd(arg,3)
size=subwrd(arg,4)
mag=subwrd(arg,5)
*text=subwrd(arg,6)
text=parsestr(arg,6)
if(size='')
usage()
return
endif
if(mag='')
mag=0
endif
if(text='')
text='Variable'
endif

qdims()

'drmrktmp='var

'query undef'
undef=subwrd(result,7)

*
* Run display to get scaling environment for gr2xy.
*
'set cmax -1e30'
*'display drmrktmp'
'display lon'
'query gxinfo'
line3=sublin(result,3)
line4=sublin(result,4)
line5=sublin(result,5)
xa=subwrd(line3,4)
ya=subwrd(line4,4)
xaxis=subwrd(line5,3)
* Note: if xyrev off then xaxis='Lon'; if xyrev on then xaxis='Lat'; to be used later.

*
* Draw mark.
*
xcnt=math_int(_xs)
while(xcnt<=_xe)
ycnt=math_int(_ys)
while(ycnt<=_ye)
'set x 'xcnt
'set y 'ycnt
'display drmrktmp'
value=subwrd(result,4)
if(value!=undef)
if(mag!=0); sizescaled=size*math_sqrt(math_abs(value)/mag); else; sizescaled=size; endif
if(xaxis='Lon')
'query gr2xy 'xcnt' 'ycnt
else
'query gr2xy 'ycnt' 'xcnt
endif
x=subwrd(result,3)
y=subwrd(result,6)
'set line 'color
'draw mark 'mark' 'x' 'y' 'sizescaled
endif
ycnt=ycnt+1
endwhile
xcnt=xcnt+1
endwhile

'set x '_xs' '_xe
'set y '_ys' '_ye

*
* Define spacing.
*
small_spacing=0.05
'query pp2xy 0 0'
tmpxa=subwrd(result,3)
'query pp2xy 1 1'
tmpxb=subwrd(result,3)
rvratio=tmpxb-tmpxa
small_spacing=small_spacing*rvratio

*
* Draw legend.
*
if(mag!=0)
mag_bigger=mag*2
mag_smaller=mag/2
size_bigger=size*math_sqrt(mag_bigger/mag)
size_smaller=size*math_sqrt(mag_smaller/mag)
'set line 'color
* Bigger mark
'draw mark 'mark' 'xa+small_spacing+size_bigger/2' 'ya+small_spacing+size_bigger/2' 'size_bigger
'set string 1 l'
'draw string 'xa+2*small_spacing+size_bigger' 'ya+small_spacing+size_bigger/2' 'mag_bigger
'query string 'mag_bigger
strwid_bigger=subwrd(result,4)
* Middle mark
'draw mark 'mark' 'xa+3*small_spacing+size_bigger+strwid_bigger+size/2' 'ya+small_spacing+size_bigger/2' 'size
'set string 1 l'
'draw string 'xa+4*small_spacing+size_bigger+strwid_bigger+size' 'ya+small_spacing+size_bigger/2' 'mag
'query string 'mag
strwid_middle=subwrd(result,4)
* Smaller mark
'draw mark 'mark' 'xa+5*small_spacing+size_bigger+strwid_bigger+size+strwid_middle+size_smaller/2' 'ya+small_spacing+size_bigger/2' 'size_smaller
'set string 1 l'
'draw string 'xa+6*small_spacing+size_bigger+strwid_bigger+size+strwid_middle+size_smaller' 'ya+small_spacing+size_bigger/2' 'mag_smaller
'query string 'mag_smaller
strwid_smaller=subwrd(result,4)
'set line 1'
'draw rec 'xa' 'ya' 'xa+7*small_spacing+size_bigger+strwid_bigger+size+strwid_middle+size_smaller+strwid_smaller' 'ya+2*small_spacing+size_bigger
endif

*
* Save mark information for use by legend.gs (in case several variables are plotted on same figure)
*
num_var=1
cnt=1
while(cnt<=num_var)
line=mark' 'color' 'size' 'text
rc=write(mytmpdir'/drawmark.txt~',line)
cnt=cnt+1
endwhile
rc=close(mytmpdir'/drawmark.txt~')

return
***************************************************************************************
function usage()
*
* Print usage information.
*
say '  Draw marks at grid points with mark size (i.e., area) proportional to the magnitude of data.'
say ''
say '  Usage: drawmark <var> <mark> <color> <size> [<magnitude> [<text>]]'
say '     <var>: variable name.'
say '     <mark>: mark type.'
say '     <color>: mark color.'
say '     <size>: mark size.'
say '     <magnitude>: magnitude of <var> corresponding to <size>. If =0 then all marks will have the same size of <size>. Defaults to 0.'
say '     <text>: Text to be shown in the legend (see "legend2.gs"). Text beginning with a minus sign or containing spaces must be double quoted.'
say ''
say '  Note 1: Must have varying X and Y dimensios, and fixed Z/T/E dimensions.'
say '  Note 2: <var> must be on a grid compatible with the default file. If not, use "set dfile" to change the default file.'
say ''
say '  Dependencies: qdims.gsf'
say ''
say '  See also: legend2.gs'
say ''
say '  Copyright (C) 2009 Bin Guan.'
say '  Distributed under GNU/GPL.'
return
