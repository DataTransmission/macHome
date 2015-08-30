***************************************************************************************
*	$Id: ppp.gs,v 1.26 2009/11/17 23:41:13 bguan Exp bguan $
*	Copyright (C) 2006 Bin Guan.
*	Distributed under GNU/GPL.
***************************************************************************************
function ppp(arg)
*
* Produce image output in EPS and other formats.
*
outfile=subwrd(arg,1)
fmt1=subwrd(arg,2)
fmt2=subwrd(arg,3)
fmt3=subwrd(arg,4)
fmt4=subwrd(arg,5)
fmt5=subwrd(arg,6)
if(outfile='')
usage()
return
endif
if(fmt1='')
fmt1='eps'
endif

tmpdir='/tmp'
whoamifile='./.whoami.bGASL'
'!whoami>'whoamifile
whoami=sublin(read(whoamifile),2)
rc=close(whoamifile)
'!unlink 'whoamifile
mytmpdir=tmpdir'/bGASL-'whoami
'!mkdir -p 'mytmpdir

'enable print 'mytmpdir'/pppout.gmf~'
'print'
'disable print'

*
* Produce EPS file.
*
if(fmt1='eps' | fmt2='eps' | fmt3='eps' | fmt4='eps' | fmt5='eps')
'!gxeps -c -i 'mytmpdir'/pppout.gmf~ -o 'mytmpdir'/pppout.eps~'
'!gs -dBATCH -dNOPAUSE -q -sDEVICE=bbox 'mytmpdir'/pppout.eps~ 2>'mytmpdir'/pppout.epsbb~'
epsbb=sublin(read(mytmpdir'/pppout.epsbb~'),2)
rc=close(mytmpdir'/pppout.epsbb~')
'!gxeps -c -s -i 'mytmpdir'/pppout.gmf~ -o 'outfile'.eps'
'!sed -i -e s/^%%BoundingBox:.\*\$/"'epsbb'"/ 'outfile'.eps'
say 'PPP>INFO: 'outfile'.eps generated.'
endif

*
* Produce EPS clip.
*
if(fmt1='eps.clip' | fmt2='eps.clip' | fmt3='eps.clip' | fmt4='eps.clip' | fmt5='eps.clip')
'!gxeps -c -i 'mytmpdir'/pppout.gmf~ -o 'outfile'.eps.clip'
'!gs -dBATCH -dNOPAUSE -q -sDEVICE=bbox 'outfile'.eps.clip 2>'mytmpdir'/pppout.epsbb~'
epsbb=sublin(read(mytmpdir'/pppout.epsbb~'),2)
rc=close(mytmpdir'/pppout.epsbb~')
'!sed -i -e s/^%%BoundingBox:.\*\$/"'epsbb'"/ 'outfile'.eps.clip'
say 'PPP>INFO: 'outfile'.eps.clip generated.'
endif

*
* Produce PDF (clip).
*
if(fmt1='pdf' | fmt2='pdf' | fmt3='pdf' | fmt4='pdf' | fmt5='pdf')
'!gxeps -c -i 'mytmpdir'/pppout.gmf~ -o 'mytmpdir'/pppout.eps.clip~'
'!gs -dBATCH -dNOPAUSE -q -sDEVICE=bbox 'mytmpdir'/pppout.eps.clip~ 2>'mytmpdir'/pppout.epsbb~'
epsbb=sublin(read(mytmpdir'/pppout.epsbb~'),2)
rc=close(mytmpdir'/pppout.epsbb~')
'!sed -i -e s/^%%BoundingBox:.\*\$/"'epsbb'"/ 'mytmpdir'/pppout.eps.clip~'
'!gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dEPSCrop -sOutputFile='outfile'.pdf 'mytmpdir'/pppout.eps.clip~'
say 'PPP>INFO: 'outfile'.pdf generated.'
endif

*
* Produce PNG (clip).
*
if(fmt1='png' | fmt2='png' | fmt3='png' | fmt4='png' | fmt5='png')
'!gxeps -c -i 'mytmpdir'/pppout.gmf~ -o 'mytmpdir'/pppout.eps.clip~'
'!gs -dBATCH -dNOPAUSE -q -sDEVICE=bbox 'mytmpdir'/pppout.eps.clip~ 2>'mytmpdir'/pppout.epsbb~'
epsbb=sublin(read(mytmpdir'/pppout.epsbb~'),2)
rc=close(mytmpdir'/pppout.epsbb~')
'!sed -i -e s/^%%BoundingBox:.\*\$/"'epsbb'"/ 'mytmpdir'/pppout.eps.clip~'
'!eps2png -png256 -scale 1.2 -output 'outfile'.png 'mytmpdir'/pppout.eps.clip~'
say 'PPP>INFO: 'outfile'.png generated.'
endif
***************************************************************************************
function usage()
*
* Print usage information.
*
say '  Produce image output in EPS and other formats.'
say ''
say '  Usage: ppp <outfile> [<format1>] [<format2>]...'
say '    <outfile>: full path of output file. Do NOT include the suffix (e.g., use mypath/myfile, instead of mypath/myfile.eps).'
say '    <format>: currently eps, eps.clip, pdf and png are supported; defaults to eps. The only difference between eps and eps.clip is that the former has a timestamp in the header (for printing purpose).'
say ''
say '  Note: png format requires eps2png. See http://search.cpan.org/~jv/eps2png/ for more information.'
say ''
say '  Copyright (C) 2006 Bin Guan.'
say '  Distributed under GNU/GPL.'
