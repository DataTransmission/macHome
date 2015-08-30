***************************************************************************************
*       $Id: subplot.gs,v 1.26 2010/03/28 05:31:11 bguan Exp bguan $
*       Copyright (C) 2009 Bin Guan.
*       Distributed under GNU/GPL.
***************************************************************************************
function subplot(arg)
*
* Prepare GrADS for a multi-panel plot.
*
rc=gsfallow('on')

*
* Read input/initialize.
*
ntot=subwrd(arg,1)
idx=subwrd(arg,2)
if(idx='')
usage()
return
endif
wrd3=subwrd(arg,3)
if(valnum(wrd3))
*evaluate wrd3
*0 is not a #
*1 is an integer
*2 is not an integer
ncol=wrd3
else
ncol=2
endif
nrow=ntot/ncol
if(nrow!=math_int(nrow))
nrow=math_int(nrow)+1
endif
*_.islandscape.1=''
_.istall.1=0
_.istight.1=0
_.isxtight.1=0
_.isytight.1=0
_.scalefactor.1=1
_.xyratio.1=0
_.morespacex.1=0
_.morespacey.1=0
_.morepadx.1=0
_.morepady.1=0
_.sameh2left.1=0
_.samey2left.1=0
*rc=parseopt(arg,'-','landscape','islandscape')
rc=parseopt(arg,'-','tall','istall')
rc=parseopt(arg,'-','tight','istight')
rc=parseopt(arg,'-','xtight','isxtight')
rc=parseopt(arg,'-','ytight','isytight')
rc=parseopt(arg,'-','scale','scalefactor')
rc=parseopt(arg,'-','xy','xyratio')
rc=parseopt(arg,'-','xs','morespacex')
rc=parseopt(arg,'-','ys','morespacey')
rc=parseopt(arg,'-','xp','morepadx')
rc=parseopt(arg,'-','yp','morepady')
rc=parseopt(arg,'-','sameh2left','sameh2left')
rc=parseopt(arg,'-','samey2left','samey2left')

*
* Get aspect ratio.
*
qdims()
'set z 1'
'set t 1'
* If z/t dimension was varying, then variables defined would only be available to the specific z/t segment.
* Setting z/t to non-varying will make variables defined available to all z's and t's. This is not an issue for x/y dimension.
xy_ratio=5/3
'query gxinfo'
line6=sublin(result,6)
mproj=subwrd(line6,3)
if(mproj=2 & _lone!=_lons & _late!=_lats)
multi_factor=(5/3)/(360/180)
lonlat_ratio=(_lone-_lons)/(_late-_lats)
xy_ratio=multi_factor*lonlat_ratio
endif
if(mproj!=2&mproj!=1&_lone!=_lons&_late!=_lats)
if(idx=1)
'display lon'
'clear'
endif
'query gxinfo'
line3=sublin(result,3)
line4=sublin(result,4)
tmpa=subwrd(line3,4)
tmpb=subwrd(line3,6)
tmpc=subwrd(line4,4)
tmpd=subwrd(line4,6)
xy_ratio=(tmpb-tmpa)/(tmpd-tmpc)
endif
if(_.xyratio.1)
xy_ratio=_.xyratio.1
endif

*
* Set up margins/spacing/padding.
*
'set vpage off'
'set parea off'
'query gxinfo'
line2=sublin(result,2)
realpagewid=subwrd(line2,4)
realpagehgt=subwrd(line2,6)
*if(_.islandscape.1=0&realpagewid=11)
*realpagewid=(8.5/11)*realpagehgt
*endif
*if(_.islandscape.1=1&realpagewid=8.5)
*realpagehgt=realpagewid/(11/8.5)
*endif
scaledpagewid=realpagewid*_.scalefactor.1
scaledpagehgt=realpagehgt*_.scalefactor.1
marginleft=0.25
marginright=0.25
margintop=0.5
marginbottom=0.5
spacex0=0
spacey0=0
padx0=0.33
pady0=0.22
spacex=spacex0+_.morespacex.1
spacey=spacey0+_.morespacey.1
padx=padx0+_.morepadx.1
pady=pady0+_.morepady.1
if(_.istight.1)
spacex=-2*padx
spacey=-2*pady
endif
if(_.isxtight.1)
spacex=-2*padx
endif
if(_.isytight.1)
spacey=-2*pady
endif

*
* Get virtual page width and height.
*
if(_.istall.1=0)
tmpwid=(scaledpagewid-marginleft-marginright-2*padx*ncol-spacex*(ncol-1))/ncol
'define vpagewid='tmpwid
'define vpagehgt=vpagewid/'xy_ratio
else
tmphgt=(scaledpagehgt-margintop-marginbottom-2*pady*nrow-spacey*(nrow-1))/nrow
'define vpagehgt='tmphgt
'define vpagewid=vpagehgt*'xy_ratio
endif
'define vpagewidpadded=vpagewid+2*'padx
'define vpagehgtpadded=vpagehgt+2*'pady
if(_.sameh2left.1)
idx_left=idx-nrow
'query defval vpageya'idx_left' 1 1'
vpage_ya_left=subwrd(result,3)
'query defval vpageyb'idx_left' 1 1'
vpage_yb_left=subwrd(result,3)
'define vpagehgtpadded='vpage_ya_left'-'vpage_yb_left
'define vpagehgt=vpagehgtpadded-2*'pady
'define vpagewid=vpagehgt*'xy_ratio
'define vpagewidpadded=vpagewid+2*'padx
endif

*
* Get virtual page boundaries.
*
row_coordinate=idx-math_int((idx-1)/nrow)*nrow
col_coordinate=math_int((idx-1)/nrow)+1
if(idx=1)
'define vpagexa'idx'=0+'marginleft
'define vpageya'idx'='realpagehgt'-'margintop
endif
if(idx>1&idx<=nrow)
idx_above=idx-1
'define vpagexa'idx'=vpagexa'idx_above
'define vpageya'idx'=vpageyb'idx_above'-'spacey
endif
if(idx>nrow&row_coordinate=1)
idx_left=idx-nrow
'define vpagexa'idx'=vpagexb'idx_left'+'spacex
'define vpageya'idx'=vpageya'idx_left
endif
if(idx>nrow&row_coordinate!=1&_.samey2left.1=0)
idx_above=idx-1
idx_left=idx-nrow
'define vpagexa'idx'=vpagexb'idx_left'+'spacex
'define vpageya'idx'=vpageyb'idx_above'-'spacey
endif
if(idx>nrow&row_coordinate!=1&_.samey2left.1=1)
idx_left=idx-nrow
'define vpagexa'idx'=vpagexb'idx_left'+'spacex
'define vpageya'idx'=vpageya'idx_left
endif
'define vpagexb'idx'=vpagexa'idx'+vpagewidpadded'
'define vpageyb'idx'=vpageya'idx'-vpagehgtpadded'

*
* Set virtual page boundaries.
*
'query defval vpagexa'idx' 1 1'
vpage_xa=subwrd(result,3)
'query defval vpagexb'idx' 1 1'
vpage_xb=subwrd(result,3)
'query defval vpageya'idx' 1 1'
vpage_ya=subwrd(result,3)
'query defval vpageyb'idx' 1 1'
vpage_yb=subwrd(result,3)
'set vpage 'vpage_xa' 'vpage_xb' 'vpage_yb' 'vpage_ya

*
* Set plotting area.
*
'query gxinfo'
line2=sublin(result,2)
psudopagewid=subwrd(line2,4)
psudopagehgt=subwrd(line2,6)
'query defval vpagewidpadded 1 1'
vpagewidpadded=subwrd(result,3)
'query defval vpagehgtpadded 1 1'
vpagehgtpadded=subwrd(result,3)
'query defval vpagewid 1 1'
vpagewid=subwrd(result,3)
'query defval vpagehgt 1 1'
vpagehgt=subwrd(result,3)
if(psudopagewid=realpagewid)
rvratio=realpagewid/vpagewidpadded
else
rvratio=realpagehgt/vpagehgtpadded
endif
parea_wid=vpagewid*rvratio
parea_hgt=vpagehgt*rvratio
parea_xa=(psudopagewid-parea_wid)/2
parea_xb=psudopagewid-(psudopagewid-parea_wid)/2
parea_ya=psudopagehgt-(psudopagehgt-parea_hgt)/2
parea_yb=(psudopagehgt-parea_hgt)/2
if(parea_xa<0);parea_xa=0;endif
if(parea_xb>=psudopagewid);parea_xb=math_format('%8.6f',psudopagewid-1e-6);endif
if(parea_ya>=psudopagehgt);parea_ya=math_format('%8.6f',psudopagehgt-1e-6);endif
if(parea_yb<0);parea_yb=0;endif
'set parea 'parea_xa' 'parea_xb' 'parea_yb' 'parea_ya

*
* Set label sizes (optional).
*
'set clopts -1 -1 '0.07*rvratio
'set xlopts 1 4 '0.09*rvratio
'set ylopts 1 4 '0.09*rvratio
'set strsiz '0.11*rvratio
say 'SUBPLOT>INFO: set clopts -1 -1 '0.07*rvratio
say 'SUBPLOT>INFO: set xlopts 1 4 '0.09*rvratio
say 'SUBPLOT>INFO: set ylopts 1 4 '0.09*rvratio
say 'SUBPLOT>INFO: set strsiz '0.11*rvratio

*
* Set xlab, ylab (optional).
*
if(_.istight.1=1 | _.isxtight.1=1)
if(col_coordinate=1)
'set ylab on'
else
'set ylab off'
endif
endif
if(_.istight.1=1 | _.isytight.1=1)
if(row_coordinate=nrow)
'set xlab on'
else
'set xlab off'
endif
endif

*
* Reset vertical/time dimension.
*
'set z '_zs' '_ze
'set t '_ts' '_te

return
***************************************************************************************
function usage()
*
* Print usage information.
*
say '  Prepare GrADS for a multi-panel plot.'
say ''
say '  USAGE: subplot <ntot> <idx> [<ncol>] [-tall 0|1] [-tight 0|1] [-xtight 0|1] [-ytight 0|1] [-scale <scalingfactor>] [-xy <xyratio>] [-xs <xspacing>] [-ys <yspacing>] [-xp <xpadding>] [-yp <ypadding>] [-sameh2left 0|1] [-samey2left 0|1]'
say '     <ntot>: total number of panels to be plotted; no preset limit. Do NOT have to be # of rows times # of columns; will be rounded up to that value.'
say '     <idx>: index of the panel to be plotted, numbered column-wise, from top to bottom, then left to right.'
say '            In each row/column, panels with smaller <idx> MUST be plotted earlier; no constraints otherwise.'
say '     <ncol>: number of columns; no preset limit. Defaults to 2 (even if <ntot>=1).'
say '     -tall 1: use if the plot is so tall that cannot otherwise be fitted into one page.'
say '     -tight 1: use if no spaces are wanted between panels.'
say '     -xtight 1: use if no horizontal spaces are wanted between panels.'
say '     -ytight 1: use if no vertical spaces are wanted between panels.'
say '     <scalingfactor>: scaling factor. Defaults to 1.'
say '     <xyratio>: aspect ratio of the plotting area. Defaults to 1.667. An optimal value will be calculated for latlon map projections.'
say '     <xspacing>: horizontal spacing (in addition to initial value 0). Defaults to 0.'
say '     <yspacing>: vertical spacing (in addition to initial value 0). Defaults to 0.'
say '     <xpadding>: horizontal padding (in addition to initial value 0.33). Defaults to 0.'
say '     <ypadding>: vertical padding (in addition to initial value 0.22). Defaults to 0.'
say '     -sameh2left 1: to set panel height the same to the immediate left panel.'
say '     -samey2left 1: to align panel top to the immediate left panel.'
say ''
say '  Note: Spacing refers to blank space between virtual pages; can be any value.'
say '        Padding refers to space between virtual page boundaries and plotting area; cannot be negative values.'
say ''
say '  EXAMPLE (2 rows by 2 columns):'
say '     set lon 120 300'
say '     set lat -20 60'
say '     subplot 4 1'
say '     display sst1'
say '     subplot 4 2'
say '     display sst2'
say '     subplot 4 3'
say '     display sst3'
say '     subplot 4 4'
say '     display sst4'
say ''
say '  Dependencies: parsestr.gsf, parseopt.gsf, qdims.gsf'
say ''
say '  Copyright (C) 2009 Bin Guan.'
say '  Distributed under GNU/GPL.'
return
