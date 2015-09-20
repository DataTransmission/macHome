%---------------------------------------------------------------------------
%Please  be  aware that the use of this LaTeX class file is governed by the
%following conditions:
%
% based on the original LaTeX ARTICLE DOCUMENT STYLE
% Copyright (C) 1988, 1989 by Leslie Lamport
%
% Copyright (c) 2007 John Wiley & Sons, Ltd, The Atrium, Southern Gate, Chichester,
% West Sussex, PO19 8SQ UK.   All rights reserved.
%
% Rules of Use
%
% You are NOT ALLOWED to change this file.
%
%This class file is made available for use by authors who wish to prepare an
%article for publication in the
%Quarterly Journal of the Royal Meteorological Society
%published by John Wiley & Sons Ltd.  The user may not exploit any part of
%the class file commercially.
%
%This class file is provided on an `as is'  basis, without warranties of any
%kind, either expressed or implied, including but not limited to warranties of
%title,   or  implied  warranties  of  merchantablility  or  fitness  for  a
%particular purpose.  There will be no duty on the author[s] of the software
%or  John Wiley & Sons Ltd to correct any errors or defects in the software.
%Any  statutory  rights you may have remain unaffected by your acceptance of
%these rules of use.
%---------------------------------------------------------------------------
%
% Created by Alistair Smith, Sunrise Setting Ltd, 30 July 2010
%
% 14/10/10 Small change to work with BibTeX
% 18/5/11 Minor style changes.
% 4/10/13 Font sizes and margin changes, plus other minire style changes. 
%
% qjrms4.cls --- For Q. J. R. Meteorol. Soc

\def\update{2013/10/14 v1.1}

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{qjrms4}[\update\ Q. J. R. Meteorol. Soc]

\newif\if@timesfont
\DeclareOption{times}{%
\@timesfonttrue}

\newif\if@doublespace
\DeclareOption{doublespace}{%
\@doublespacetrue}

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions
\LoadClass{article}
\if@timesfont
\RequirePackage{times}
\fi
\if@doublespace
\RequirePackage[doublespacing]{setspace}
\fi

\RequirePackage{graphicx}
\RequirePackage{pifont,latexsym,ifthen,theorem,rotating,calc,textcase,booktabs,color}
\RequirePackage{amsfonts,amssymb,amsbsy,amsmath}
\RequirePackage[errorshow]{tracefnt}

%\usepackage{natbib}
%\bibpunct{(}{)}{;}{}{}{,}

\@twosidetrue
\flushbottom
\frenchspacing

%papersize={width,height}
\usepackage[text={179mm,268mm},%
papersize={210mm,297mm},%
columnsep=4mm,%
headsep=29pt,%
centering]{geometry}

\headheight 19pt
\parindent 1em
\footskip 18pt

\hyphenpenalty=1000
\pretolerance=8000
\tolerance=9500
\hbadness=8000
\vbadness=9000
\displaywidowpenalty=0
\clubpenalty=10000
\widowpenalty=10000
\lefthyphenmin=3%
\righthyphenmin=3%
\brokenpenalty=10000%
\hbadness=10000
\vbadness=10000

\thinmuskip = 3mu
\medmuskip = 4mu
\thickmuskip = 5mu

\setcounter{topnumber}{5}
\def\topfraction{1}
\setcounter{bottomnumber}{5}
\def\bottomfraction{.8}
\setcounter{totalnumber}{10}
\def\textfraction{0}
\renewcommand{\floatpagefraction}{0.95}
\renewcommand{\dblfloatpagefraction}{0.95}
\renewcommand{\dbltopfraction}{1}

\renewcommand\normalsize{%
\fontsize{9.5}{11.5pt plus 0.10pt minus 0.10pt}\selectfont%
   %\@setfontsize\normalsize\@xipt{12pt plus 0.10pt minus 0.10pt}%
   \abovedisplayskip 11\p@ \@plus3\p@ \@minus6\p@
   \abovedisplayshortskip \z@ \@plus3\p@
   \belowdisplayshortskip 6.5\p@ \@plus3.5\p@ \@minus3\p@
   \belowdisplayskip \abovedisplayskip
   \let\@listi\@listI}
\normalsize
\renewcommand\small{%
\fontsize{8.5}{9.5pt plus 0.10pt minus 0.10pt}\selectfont%
   %\@setfontsize\small\@xpt{11pt plus 0.10pt minus 0.10pt}
   \abovedisplayskip 10\p@ \@plus2\p@ \@minus5\p@
   \abovedisplayshortskip \z@ \@plus3\p@
   \belowdisplayshortskip 6\p@ \@plus3\p@ \@minus3\p@
   \def\@listi{\leftmargin\leftmargini
               \topsep 6\p@ \@plus2\p@ \@minus2\p@
               \parsep 3\p@ \@plus2\p@ \@minus\p@
               \itemsep \parsep}%
   \belowdisplayskip \abovedisplayskip
}
\renewcommand\footnotesize{%
\fontsize{7.5}{8.5}\selectfont%
   %\@setfontsize\footnotesize\@viiipt{9}%
   \abovedisplayskip 8\p@ \@plus2\p@ \@minus4\p@
   \abovedisplayshortskip \z@ \@plus\p@
   \belowdisplayshortskip 4\p@ \@plus2\p@ \@minus2\p@
   \def\@listi{\leftmargin\leftmargini
               \topsep 4\p@ \@plus2\p@ \@minus2\p@
               \parsep 2\p@ \@plus\p@ \@minus\p@
               \itemsep \parsep}%
   \belowdisplayskip \abovedisplayskip
}
\renewcommand{\scriptsize}{\fontsize{7.0}{8.5pt}\selectfont}
\renewcommand\tiny{\@setfontsize\tiny\@vipt\@viipt}
%\renewcommand\large{\@setfontsize\large\@xiipt{14}}
\renewcommand\Large{\@setfontsize\Large\@xivpt{18}}
\renewcommand\LARGE{\@setfontsize\LARGE\@xviipt{22}}
\renewcommand\huge{\@setfontsize\huge\@xxpt{25}}
\renewcommand\Huge{\@setfontsize\Huge\@xxvpt{30}}

\renewcommand{\large}{\fontsize{12.5}{15pt}\selectfont}
\newcommand{\refsize}{\fontsize{8.0}{9.5pt}\selectfont}
\newcommand{\rhsize}{\fontsize{11.0}{13pt}\selectfont}
\newcommand{\titlesize}{\fontsize{17.8}{19.5pt}\selectfont}
\newcommand{\abssize}{\fontsize{11.0}{12.5pt}\selectfont}

\def\journalname{QUARTERLY JOURNAL OF THE ROYAL METEOROLOGICAL SOCIETY}
\def\journalabb{Q.~J.~R. Meteorol. Soc.}

\def\address#1{\gdef\@address{#1}}

\newbox\absbox
\def\abstract{\lrbox\absbox\minipage{143mm}%
  \abssize\bfseries%
  }
\def\endabstract{\endminipage\endlrbox}

\def\keywords#1{%
  \gdef\@keywords{%
    \noindent\emph{Key Words:}\hspace*{1.0em}\ignorespaces#1}}
\let\@keywords\@empty

\def\received#1{%
\gdef\@received{\noindent Received #1}}
\newcommand\revised[2][Revised]{%
\g@addto@macro\@revised{{; #1 #2}}}
\let\@revised\@empty
\newcommand\accepted[2][Accepted]{%
\g@addto@macro\@accepted{{; #1 #2}}}
\let\@accepted\@empty

\setlength{\skip\footins}{17pt}
%\renewcommand{\footnoterule}{\rule{0pc}{0pt}}
\def\footnoterule{\kern-3\p@
  \hrule \@width 60pt \kern 2.6\p@}

\renewcommand{\thefootnote}{\fnsymbol{footnote}}
\long\def\@makefntext#1{\parindent 1em%
\noindent{$\m@th^{\@thefnmark}$}#1}

\def\corraddr#1{%
  \gdef\@corraddr{%
   {$^{*}$Correspondence to: #1}}}
\let\@corraddr\@empty
\def\corrauth{$^{*}$}

\def\affilnum#1{${}^{\text{\normalfont#1}}$}
\def\affil#1{${}^{\text{\normalfont#1}}$}
\def\comma{$^{\text{,}}$}

\renewcommand\maketitle{\par
  \begingroup
     \if@twocolumn
      \ifnum \col@number=\@ne
        \@maketitle
      \else
        \twocolumn[\@maketitle]%
       \fi
    \else
      \newpage
      \global\@topnum\z@   % Prevents figures from going at top of page.
      \@maketitle
    \fi
    \thispagestyle{title}\label{FirstPage}%\@corraddr
  \endgroup
  %\setcounter{footnote}{0}%
  \global\let\address\relax
  \global\let\thanks\relax
  \global\let\maketitle\relax
  \global\let\@maketitle\relax
  \global\let\@thanks\@empty
  \global\let\@author\@empty
  \global\let\@date\@empty
  \global\let\@title\@empty
  \global\let\@address\@empty
  \global\let\corraddr\relax
  \global\let\title\relax
  \global\let\author\relax
  \global\let\date\relax
  \global\let\and\relax
}
\def\@maketitle{%
  \null
  \if@doublespace
\vskip 0pt
\else  
  \vskip 4pt %
\fi
\mbox{}\hfill\framebox[43mm]{\rule{0pt}{15mm}}\hfill\mbox{}
\rule{\textwidth}{1.0pt}
  \begin{center}%
  \vskip 2pt
    {\titlesize\bfseries \@title \par}%
    \vskip 1.75em%
    {\large
      \lineskip .5em%
      \@author
      \par}%
      \vskip 2pt
 {\normalsize
      \lineskip .5em%
      %
        {\emph\@address}
      \par}
      \vskip 6pt
      {\@corraddr}
      % 
    \vskip 6pt
 \rule{\textwidth}{1.0pt}    
    \vskip 17pt%
    {\usebox\absbox\par}%
    \begin{minipage}{143mm}
   {\lineskip .5em%
   \vskip 14pt
      {\small\@keywords}
      \par}
  \vskip 7pt
  {\lineskip .5em%
     {\small\emph{\@received\@revised\@accepted}}
      \par}
  \vskip 7pt
  \end{minipage}
      \end{center}\par
  \vskip 1.5em}
 
\def\startpage{\pageref{FirstPage}}
\def\endpage{\pageref{LastPage}}
\def\volumeyear{0000}
\def\volumenumber{00}
\def\DOI{qj.000}

\gdef\copyrightline{\copyright\ \volumeyear\ Royal Meteorological Society}
\def\runningheads#1#2{\markboth{{#1}}{{#2}}}

\definecolor{QJGray}{gray}{0.937}

\def\ps@title{\setlength{\fboxsep}{0pt}%
\def\@oddhead{\parbox[t]{\textwidth}{%
\noindent\colorbox{QJGray}{\parbox{\textwidth}{\mbox{}\\[3.5mm]%
\noindent\small\bfseries Quarterly Journal of the Royal Meteorological Society\hfill%
\emph{\journalabb} \textbf{\volumenumber}:\ \startpage--\endpage\ (\volumeyear)\\[1.4mm]}}}}%
\let\@evenhead\@oddhead
\def\@oddfoot{\parbox[t]{\textwidth}{%
\noindent\parbox{\columnwidth}{\footnotesize\mbox{}\hfill\copyrightline\hfill\mbox{}}\hfill\mbox{}%\\
\noindent\parbox{\columnwidth}{\footnotesize\mbox{}\hfill\textit{Prepared using \textsf{qjrms4.cls} [Version: \update]}\hfill\mbox{}}\hfill\mbox{}}%
}}

\def\ps@qjrms{\setlength{\fboxsep}{0pt}%
\def\@evenhead{\parbox[t]{\textwidth}{%
\noindent\colorbox{QJGray}{\parbox{\textwidth}{\mbox{}\\[2.5mm]%
\noindent{\rhsize\bfseries\thepage\hfill\leftmark\hfill\phantom{\thepage}}\\[0.5mm]}}}}
\def\@oddhead{\parbox[t]{\textwidth}{%
\noindent\colorbox{QJGray}{\parbox{\textwidth}{\mbox{}\\[2.5mm]%
\noindent{\rhsize\bfseries\phantom{\thepage}\hfill\rightmark\hfill\thepage}\\[0.5mm]}}}}
\def\@evenfoot{\parbox[t]{\textwidth}{%
\noindent\parbox{\columnwidth}{\footnotesize\mbox{}\hfill\copyrightline\hfill\mbox{}}\hfill%
\noindent\parbox{\columnwidth}{\footnotesize\mbox{}\hfill\emph{Prepared using \textsf{qjrms4.cls}}\hfill\mbox{}}}}%
\def\@oddfoot{\@evenfoot}}

\let\@hangfrom\noindent %Remove indent on section turnovers
\renewcommand{\@seccntformat}[1]{{\csname the#1\endcsname.}\hspace{1.0em}}

\renewcommand\section{\@startsection {section}{1}{\z@}%
                                   {-12pt \@plus -6pt \@minus -4pt}%
                                   {10pt \@plus2pt \@minus 2pt}%
                                   {\normalfont\normalsize\raggedright\bfseries}}
\renewcommand\subsection{\@startsection{subsection}{2}{\z@}%
                                     {-12pt \@plus -6pt \@minus -4pt}%
                                     {10pt \@plus2pt \@minus 2pt}%
                                     {\normalfont\normalsize\raggedright\itshape}}
\renewcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}%
                                     {-12pt \@plus -6pt \@minus -4pt}%
                                     {10pt \@plus2pt \@minus 2pt}%
                                     {\normalfont\normalsize\raggedright\itshape}}

\def\enumerate{\ifnum \@enumdepth >3 \@toodeep\else
      \advance\@enumdepth \@ne
      \edef\@enumctr{enum\romannumeral\the\@enumdepth}\list
      {\csname label\@enumctr\endcsname}{\usecounter
        {\@enumctr}\itemsep 0pt\parsep 0pt
         \def\makelabel##1{\hss\llap{##1}}}\fi}

\let\endenumerate =\endlist

\def\itemize{\ifnum \@itemdepth >3 \@toodeep\else \advance\@itemdepth \@ne
\edef\@itemitem{labelitem\romannumeral\the\@itemdepth}%
\list{\csname\@itemitem\endcsname}{\itemsep 0pt\parsep 0pt
\def\makelabel##1{\hss\llap{##1}}}\fi}

\let\enditemize =\endlist

%\renewcommand{\thetable}{\Roman{table}}

\usepackage{caption}
\DeclareCaptionFont{refsize}{\refsize}
\DeclareCaptionLabelSeparator{jwperiod}{.\hspace*{0.75ex}}
\captionsetup[figure]{font=footnotesize,labelfont=bf,labelsep=jwperiod,justification=justified,singlelinecheck=true}
\captionsetup[table]{position=top,font=refsize,labelfont=rm,labelsep=jwperiod,justification=justified,singlelinecheck=true}

\newcommand\ack{\section*{Acknowledgement}}
\newcommand\acks{\section*{Acknowledgements}}

\renewenvironment{thebibliography}{%
\section*{\refname}
\list{}{\topsep=0\p@\parsep=0\p@
\partopsep=0\p@\itemsep=1\p@\labelsep=0\p@\itemindent=-8\p@
\labelwidth=0\p@\leftmargin=8\p@
}\refsize\rm
\def\newblock{\ }
\sloppy\clubpenalty4000\widowpenalty4000
\sfcode`\.=1000\relax}{\endlist}

\usepackage{natbib}
\bibpunct{(}{)}{;}{}{}{,}
\setlength{\bibsep}{1pt}
\renewcommand{\bibfont}{\refsize}

\hyphenation{never-the-less}
\hyphenation{atmos-phere atmos-pheric}
\hyphenation{fila-ment}
\hyphenation{kin-etic}
\hyphenation{meso-scale}
\hyphenation{micro-wave}
\hyphenation{par-ameter parame-trize}
\hyphenation{plan-etary}
\hyphenation{spheri-cal}
\hyphenation{reas-on reas-ons}
\hyphenation{for-mula for-mulae for-mulas}
\hyphenation{meteorol}
\hyphenation{thermo-dynamics thermo-dynamic}
\hyphenation{follows}
\hyphenation{radi-ation}
\hyphenation{per-iph-ery}
\hyphenation{meteor-ology meteoro-logi-cal meteor-ol-ogist}
\hyphenation{other-wise}
\hyphenation{down-slope}
\hyphenation{wave-guide wave-guides}
\hyphenation{extra-tropical}

\AtEndDocument{%
\label{LastPage}}

\newcommand\re{\mbox{\rm e}}
\newcommand\rd{\mbox{\rm d}}
\newcommand\en{\phantom{0}}
\newcommand\etal{\emph{et~al.}}

\pagestyle{qjrms}
\twocolumn
\normalsize
\sloppy
\received{\dots}

