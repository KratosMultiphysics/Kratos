\NeedsTeXFormat{LaTeX2e}[1995/12/01]
\ProvidesClass{deltares_manual}
              [2014/04/01 v1.0 Standard LaTeX document class]
\DeclareOption{onecolumn}{\OptionNotUsed}
\makeatletter

\newif\if@ddraft% deltares draft
\@ddraftfalse
\DeclareOption{ddraft}{%
    \@ddrafttrue%
}
\newif\if@dutch
\@dutchfalse
\DeclareOption{dutch}{%
    \@dutchtrue%
}
\newif\if@spanish
\@spanishfalse
\DeclareOption{spanish}{%
    \@spanishtrue%
}
\newif\if@signature
\@signaturefalse
\DeclareOption{signature}{%
    \@signaturetrue%
}
\newif\if@druler% deltares ruler
\@drulerfalse
\DeclareOption{druler}{%
    \@drulertrue%
}
\newif\if@dbiblatex% deltares biblatex
\@dbiblatexfalse
\DeclareOption{biblatex}{%
    \@dbiblatextrue%
}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{report}}
\ProcessOptions\relax
%\LoadClass[11pt,twoside,openright,a4paper,titlepage]{report}
\LoadClass[11pt,twoside,a4paper,titlepage]{report}
%------------------------------------------------------------------------------
%
% Change to Helvetica font
%
%\usepackage{hvmaths} % commercial
%
\usepackage[scaled=.92]{helvet}
\def\rmdefault{phv}
\DeclareMathSizes{\@xipt}{\@xiipt}{8}{6}
\renewcommand{\ttdefault}{pcr} % courier lettertype, needed for bold face in verbatim environment
%
%\usepackage{frutiger} % font van geo-systems manuals
%------------------------------------------------------------------------------
% Declaring the letters as OML/phv/m/it does not change the font.
% However, upon changing it to OT1/phv/m/it characters such as
% <, > and greek characters are lost.
%
%\DeclareSymbolFont{letters}     {OML}{phv}{m}{it}
%\DeclareSymbolFont{letters}     {OT1}{phv}{m}{it}
%
%\DeclareSymbolFont{operators}   {OT1}{phv}{m}{n}
%\DeclareSymbolFont{symbols}     {OMS}{phv}{m}{n}
%\DeclareSymbolFont{largesymbols}{OMX}{phv}{m}{n}
%\DeclareSymbolFont{bold}        {OT1}{phv}{bx}{n}
%\DeclareSymbolFont{italic}      {OT1}{phv}{m}{it}
%\renewcommand{\familydefault}{phv}

%------------------------------------------------------------------------------
\usepackage{sectsty}  % to set section fonts
%\allsectionsfont{\fontfamily{cmss}\selectfont}
%\allsectionsfont{\fontfamily{pag}\selectfont} % AvantGarde
%\allsectionsfont{\fontfamily{pbk}\selectfont} % bookman
%\allsectionsfont{\fontfamily{pch}\selectfont} % Charter
%\allsectionsfont{\fontfamily{pcr}\selectfont} % Courier
%\allsectionsfont{\fontfamily{pnc}\selectfont} % New Century Schoolbook
%\allsectionsfont{\fontfamily{ppl}\selectfont} % palatino
%\allsectionsfont{\fontfamily{ptm}\selectfont} % times
%\allsectionsfont{\fontfamily{put}\selectfont} % utopia
%\allsectionsfont{\fontfamily{phv}\selectfont} % helvetica

\allsectionsfont{\fontfamily{phv}\selectfont} % helvetica

\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage[pdftex]{graphicx}
\usepackage{epstopdf}
\usepackage{color}
\usepackage{frcursive} % needed for 'write' letter (French Cursive), use in processes library
\usepackage{float} % needed to force the location of figures \begin{figure}[H]
\usepackage{textcomp} % needed for the promille sign: \textperthousand
\usepackage[final]{pdfpages} % needed for cover pages

\if@dutch
    \usepackage[dutch]{babel}
\else
    \if@spanish
        \usepackage[spanish]{babel}
    \else
        \usepackage[english]{babel}
    \fi
\fi

\usepackage{svn-multi}
\usepackage{tabularx}
\usepackage{longtable}
\usepackage{verbatim}
\usepackage{fancyvrb} % needed for framed verbatim environment
\usepackage{calc} % adding, substracting, multiplying of lengths
\usepackage{wrapfig}
\usepackage{caption}
\usepackage{subcaption} % two or more figures in a figure-environment
\usepackage{url}
\usepackage{moreverb}
\usepackage[smallscripts]{moresize}
\usepackage{rotating} % needed for rotating text in tables
\usepackage{nameref}
\usepackage{import} % needed for nested \input{...} -> \import{...}{...}
\usepackage{marginnote} % needed to place the word TODO in the margin, also when it is used in tables

\usepackage{amstext} % for \text command in math environment
%
% to align equations left add [fleqn] as option to amsmath
%
\usepackage[fleqn]{amsmath}
\usepackage{amssymb}
\usepackage{paralist} % inline enumeration

\usepackage[detect-all=true, product-units=power, detect-family=false, math-micro=\mu, math-rm=\selectfont, text-rm=\sffamily]{siunitx}
\AtBeginDocument{\sisetup{math-rm=\selectfont, text-rm=\sffamily}}

%\usepackage[noabbrev]{cleveref}

\usepackage{xspace} % needed for spaces behind a macro, but no space before punctuation
\usepackage{colortbl} % colour the cells and rows
\usepackage{lastpage}
\usepackage{etoolbox}% http://ctan.org/pkg/etoolbox: needed for conditional list of figures and list of tables
\usepackage{nomencl}
\usepackage[intoc]{nomentbl}
\usepackage{bookmark}
\usepackage{currfile} % needed to detect the directory of the current tex-file 
\hypersetup{destlabel, pagebackref, pdfborder={0 0 0}, bookmarksnumbered=true, colorlinks=true, linkcolor=dbluelink, citecolor=dbluelink, urlcolor=dbluelink}
%------------------------------------------------------------------------------

\hoffset -1in \voffset -1in
\setlength{\baselineskip}{4.5mm}
\linespread{0.941}

\topmargin 2\baselineskip
\headheight 2\baselineskip
\headsep 2\baselineskip
% until text = \topmargin + \headheight + \headsep == 6\baselineskip
\setlength{\textheight}{297mm-6\baselineskip-6\baselineskip}
\footskip 3\baselineskip
%foothoogte 2\baselineskip
% 1.5cm remaining

\oddsidemargin 35mm
\evensidemargin 30mm
\setlength{\textwidth}{210mm-65mm} % 145mm

\setlength{\marginparwidth}{18mm}
\setlength{\marginparsep}{4.5mm}
\let\oldmarginpar\marginpar
\renewcommand\marginpar[1]{\-\oldmarginpar[\raggedleft\footnotesize #1]%
{\raggedright\footnotesize #1}}

\widowpenalty=3000
\clubpenalty=3000

%------------------------------------------------------------------------------
%
% front- and mainmatter commands as known from book.cls, now extended with
% pagestyle, such that page numbers are: i,ii, and 1 of 10, 2 of 10, etc.
% respectively. front- and mainmatter commands automatically set in
% \reporttitle
\newif\if@mainmatter \@mainmattertrue
\newcommand\frontmatter{%
    \cleardoublepage
    \@mainmatterfalse
    \pagenumbering{roman}%
    \pagestyle{frontmatterplain}%
    \if@druler
        \AddToShipoutPicture{\BackgroundRuler}
    \fi
}
\newcommand\mainmatter{%
    \cleardoublepage
    \@mainmattertrue
    \pagenumbering{arabic}
    \pagestyle{headings}
    \if@druler
        \AddToShipoutPicture{\BackgroundRuler}
    \fi
}
\if@druler
    \usepackage{deltares_ruler}
\fi

%------------------------------------------------------------------------------
%
% copy of tweaklist.sty, to adjust itemize, enumerate and description
%
\def\enumhook{}
\def\enumhooki{}
\def\enumhookii{}
\def\enumhookiii{}
\def\enumhookiv{}
\def\itemhook{}
\def\itemhooki{}
\def\itemhookii{}
\def\itemhookiii{}
\def\itemhookiv{}
\def\descripthook{}

\def\enumerate{%
  \ifnum \@enumdepth >\thr@@\@toodeep\else
    \advance\@enumdepth\@ne
    \edef\@enumctr{enum\romannumeral\the\@enumdepth}%
      \expandafter
      \list
        \csname label\@enumctr\endcsname
        {\usecounter\@enumctr\def\makelabel##1{\hss\llap{##1}}%
          \enumhook \csname enumhook\romannumeral\the\@enumdepth\endcsname}%
  \fi}
\def\itemize{%
  \ifnum \@itemdepth >\thr@@\@toodeep\else
    \advance\@itemdepth\@ne
    \edef\@itemitem{labelitem\romannumeral\the\@itemdepth}%
    \expandafter
    \list
      \csname\@itemitem\endcsname
      {\def\makelabel##1{\hss\llap{##1}}%
        \itemhook \csname itemhook\romannumeral\the\@itemdepth\endcsname}%
  \fi}
\renewenvironment{description}
                 {\list{}{\labelwidth\z@ \itemindent-\leftmargin
                          \let\makelabel\descriptionlabel\descripthook}}
                 {\endlist}

\renewcommand{\itemhook}{%
  \setlength{\topsep}{\baselineskip/2}%
  \setlength{\itemsep}{0pt}%
  \setlength{\labelwidth}{1em}%
  \setlength{\labelsep}{2.25mm}%
  \setlength{\leftmargin}{1em+2.25mm}%
  \setlength{\parskip}{0pt}%
}%
\renewcommand{\enumhooki}{%
  \setlength{\topsep}{\baselineskip/2}%
  \setlength{\itemsep}{0pt}%
  \setlength{\labelwidth}{1em}%
  \setlength{\labelsep}{2.25mm}%
  \setlength{\leftmargin}{1em+2.25mm}%
  \setlength{\parskip}{0pt}%
}
\renewcommand{\enumhookii}{%
  \setlength{\topsep}{\baselineskip/2}%
  \setlength{\itemsep}{0pt}%
  \setlength{\labelwidth}{2em}%
  \setlength{\labelsep}{2.25mm}%
  \setlength{\leftmargin}{2em+2.25mm}%
  \setlength{\parskip}{0pt}%
}
\renewcommand{\enumhookiii}{%
  \setlength{\topsep}{\baselineskip/2}%
  \setlength{\itemsep}{0pt}%
  \setlength{\labelwidth}{3em}%
  \setlength{\labelsep}{2.25mm}%
  \setlength{\leftmargin}{3em+2.25mm}%
  \setlength{\parskip}{0pt}%
}
\renewcommand{\enumhookiv}{%
  \setlength{\topsep}{\baselineskip/2}%
  \setlength{\itemsep}{0pt}%
  \setlength{\labelwidth}{4em}%
  \setlength{\labelsep}{2.25mm}%
  \setlength{\leftmargin}{4em+2.25mm}%
  \setlength{\parskip}{0pt}%
}
\renewcommand{\descripthook}{%
  \setlength{\topsep}{\baselineskip/2}%
  \setlength{\itemsep}{0pt}%
  \setlength{\labelwidth}{4.5mm}%
  \setlength{\labelsep}{4.5mm}%
  \setlength{\leftmargin}{9.0mm}%
  \setlength{\parskip}{0pt}%
}
%------------------------------------------------------------------------------
%
% newcommand for tableofcontents, listoffigures, listotables; taken from report.sty and macro MakeUppercase removed
% 4 nov 2008 Jan Mooiman
%
\setcounter{tocdepth}{3} % i.e. chapter, section, subsection and subsubsection

\renewcommand*\l@section{\@dottedtocline{1}{1.5em}{2.7em}}
\renewcommand*\l@subsection{\@dottedtocline{2}{4.2em}{3.6em}}
\renewcommand*\l@subsubsection{\@dottedtocline{3}{7.8em}{4.5em}}
\renewcommand*\l@paragraph{\@dottedtocline{4}{12.3em}{5.4em}}
\renewcommand*\l@subparagraph{\@dottedtocline{5}{17.7em}{6.3em}}
\renewcommand*\l@figure{\l@section}
\renewcommand*\l@table{\l@figure}

\AtEndEnvironment{figure}{\gdef\there@is@a@figure{}\label{fig:was:used:in:doc}}%
\AtEndDocument{\ifdefined\there@is@a@figure\label{fig:was:used:in:doc}\fi}%
\newcommand{\conditionalLoF}{\@ifundefined{r@fig:was:used:in:doc}{}{\listofFigures}}%

\AtBeginEnvironment{longtable}{\gdef\there@is@a@table{}\label{table:was:used:in:doc}}%
\AtEndEnvironment{table}{\gdef\there@is@a@table{}\label{table:was:used:in:doc}}%
\AtEndDocument{\ifdefined\there@is@a@table\label{table:was:used:in:doc}\fi}%
\newcommand{\conditionalLoT}{\@ifundefined{r@table:was:used:in:doc}{}{\listofTables}}%

\renewcommand\tableofcontents{%
    \if@twocolumn
      \@restonecoltrue\onecolumn
    \else
      \@restonecolfalse
    \fi
    \begingroup
        \chapter*{\contentsname
            \@mkboth{%
               \contentsname}{\contentsname}}%
               \parskip\z@%
        \@starttoc{toc}%
        \if@restonecol\twocolumn\fi
    \endgroup
    }
\newcommand\listofFigures{%
    \if@twocolumn
      \@restonecoltrue\onecolumn
    \else
      \@restonecolfalse
    \fi
    \begingroup
        \nonumchapter{\listfigurename}%
          \@mkboth{\listfigurename}%
                  {\listfigurename}%
          \parskip\z@%
        \@starttoc{lof}%
        \if@restonecol\twocolumn\fi
    \endgroup
    }
\newcommand\listofTables{%
    \if@twocolumn
      \@restonecoltrue\onecolumn
    \else
      \@restonecolfalse
    \fi
    \begingroup
        \nonumchapter{\listtablename}%
          \@mkboth{%
              \listtablename}%
             {\listtablename}%
          \parskip\z@%
        \@starttoc{lot}%
        \if@restonecol\twocolumn\fi
    \endgroup
    }
\newcommand\listofPhotos{%
    \if@twocolumn
      \@restonecoltrue\onecolumn
    \else
      \@restonecolfalse
    \fi
    \begingroup
        \nonumchapter{\listtablename}%
          \@mkboth{%
              \listtablename}%
             {\listtablename}%
          \parskip\z@%
        \@starttoc{lot}%
        \if@restonecol\twocolumn\fi
    \endgroup
    }
%------------------------------------------------------------------------------
\if@ddraft
  \def\draftname{}%
  \AtBeginDocument{%
  \RequirePackage{draftwatermark}%
  \SetWatermarkAngle{54.7}%
  \SetWatermarkScale{5}%
  \SetWatermarkFontSize{1cm}%
  \SetWatermarkLightness{0.92}%
  \SetWatermarkText{DRAFT}%

    \renewcommand\sc@watermark{%
        \setlength{\@tempdimb}{.63\paperwidth}%
        \setlength{\@tempdimc}{-.58\paperheight}%
        \put(\strip@pt\@tempdimb,\strip@pt\@tempdimc){%
        \makebox(0,0){\rotatebox{\sc@wm@angle}{%
        \scalebox{\sc@wm@scale}{%
        \expandafter\expandafter\expandafter\textcolor
        \expandafter\sc@wm@colormodel\expandafter{\sc@wm@colorspecs}{%
        \fontsize{\sc@wm@fontsize}{10.0\sc@wm@fontsize}
        \selectfont
        \sc@wm@text}}}}}}
  }
\else
  \def\draftname{}
\fi
\def\ds@draft{\overfullrule 5pt \@drafttrue}

\if@dutch
    \def\@Version{Versie: \@version.\svnrev}
\else
    \def\@Version{Version: \@version.\svnrev}
\fi
\def\svn@revinit{00}
\let\svnrev\svn@revinit\let\@svn@rev\svn@revinit

\if@dutch
    \newcommand{\monthname}{%
        \ifcase\month\or januari \or februari%
        \or maart\or april \or mei%
        \or juni\or juli \or augustus%
        \or september \or oktober \or november%
        \or december\fi}%
\else
    \newcommand{\monthname}{%
        \ifcase\month\or January\or February%
        \or March\or April\or May%
        \or June\or July\or August%
        \or September\or October\or November%
        \or December\fi}%
\fi

\raggedbottom

\def\title#1{\gdef\@title{#1}}
\def\subtitle#1{\gdef\@subtitle{#1}}
\def\manualtype#1{\gdef\@manualtype{#1}}
\def\distribution#1{\gdef\@distribution{#1}}
\def\version#1{\gdef\@version{#1}}
\def\author#1{\gdef\@author{#1}}

\def\@title{use \texttt{$\backslash$title\{\ldots\}}}
\def\@subtitle{use \texttt{$\backslash$subtitle\{\ldots\}}}
\def\@manualtype{use \texttt{$\backslash$manualtype\{\ldots\}}}
\def\@version{use \texttt{$\backslash$version\{\ldots\}}}
\def\@distribution{}
\def\@author{}

\newcommand{\thechap}{$\;$\leftmark{}$\;$}
\renewcommand\today{%
    \number\day\space\monthname\space\number\year%
}
% New column type, to assure parskip in tabular environment (see macro \summary)
\newcolumntype{P}[1]{>{%
    \parskip=1\baselineskip%
    \advance\parskip by 0pt plus 2pt%
    \setlength{\parfillskip}{30pt plus 1fil}}p{#1}%
}
%------------------------------------------------------------------------------

%\def\tiny{\@setsize\tiny{7pt}\vipt\@vipt}
%\def\scriptsize{\@setsize\scriptsize{9.5pt}\viiipt\@viiipt}
%\def\large{\@setsize\large{14pt}\xiipt\@xiipt}
%\def\Large{\@setsize\Large{18pt}\xivpt\@xivpt}
%\def\LARGE{\@setsize\LARGE{22pt}\xviipt\@xviipt}
%\def\huge{\@setsize\huge{25pt}\xxpt\@xxpt}
%\def\Huge{\@setsize\Huge{30pt}\xxvpt\@xxvpt}
%\def\HUGE{\@setsize\Huge{40pt}\xxvpt\@xxvpt}

%\def\tiny{\@setsize\tiny{7.33pt}\vipt\@vipt}
%\def\scriptsize{\@setsize\scriptsize{8.5pt}\viiipt\@viiipt}
%\def\footnotesize{\@setsize\footnotesize{9.25pt}\viiipt\@viiipt}
%\def\small{\@setsize\small{10.00pt}\viiipt\@viiipt}
%%\def\normalsize{\@setsize\normalsize{10.95pt}\viiipt\@viiipt}
%\def\large{\@setsize\large{11.75pt}\xiipt\@xiipt}
%\def\Large{\@setsize\Large{14.10pt}\xivpt\@xivpt}
%\def\LARGE{\@setsize\LARGE{15.85}\xviipt\@xviipt}
%\def\huge{\@setsize\huge{19.02}\xxpt\@xxpt}
%\def\Huge{\@setsize\Huge{22.82}\xxvpt\@xxvpt}
%\def\HUGE{\@setsize\Huge{30pt}\xxvpt\@xxvpt}

\def\title#1{\gdef\@title{\selectfont #1}}
\def\subtitle#1{\gdef\@subtitle{\selectfont #1}}
\def\author#1{\gdef\@author{\selectfont #1}}
%\def\reporttype#1{\gdef\@reporttype{\selectfont #1}}
%\def\reportdate#1{\gdef\@reportdate{\selectfont #1}}

\long\def\approvali#1{\gdef\@approvali{#1}}
\long\def\approvalii#1{\gdef\@approvalii{#1}}
\long\def\approvaliii#1{\gdef\@approvaliii{#1}}
\long\def\approvaliv#1{\gdef\@approvaliv{#1}}
\long\def\authori#1{\gdef\@authori{#1}}
\long\def\authorii#1{\gdef\@authorii{#1}}
\long\def\authoriii#1{\gdef\@authoriii{#1}}
\long\def\authoriv#1{\gdef\@authoriv{#1}}
\long\def\classification#1{\gdef\@classification{#1}}
\long\def\client#1{\gdef\@client{#1}}
\long\def\datei#1{\gdef\@datei{#1}}
\long\def\dateii#1{\gdef\@dateii{#1}}
\long\def\dateiii#1{\gdef\@dateiii{#1}}
\long\def\dateiv#1{\gdef\@dateiv{#1}}
\long\def\keywords#1{\gdef\@keywords{#1}}
\long\def\projectnumber#1{\gdef\@projectnumber{#1}}
\long\def\reference#1{\gdef\@reference{#1}}
\long\def\references#1{\gdef\@references{#1}}
\long\def\revieweri#1{\gdef\@revieweri{#1}}
\long\def\reviewerii#1{\gdef\@reviewerii{#1}}
\long\def\revieweriii#1{\gdef\@revieweriii{#1}}
\long\def\revieweriv#1{\gdef\@revieweriv{#1}}
\long\def\status#1{\gdef\@status{#1}}
\long\def\summary#1{\gdef\@summary{#1}}
\long\def\versioni#1{\gdef\@versioni{#1}}
\long\def\versionii#1{\gdef\@versionii{#1}}
\long\def\versioniii#1{\gdef\@versioniii{#1}}
\long\def\versioniv#1{\gdef\@versioniv{#1}}

\long\def\@approvali{\tiny use \texttt{$\backslash$approvali\{\ldots\}}}
\long\def\@approvalii{}
\long\def\@approvaliii{}
\long\def\@approvaliv{}
\long\def\@authori{\tiny use \texttt{$\backslash$authori\{\ldots\}}}
\long\def\@authorii{}
\long\def\@authoriii{}
\long\def\@authoriv{}
\long\def\@classification{use \texttt{$\backslash$classification\{\ldots\}}}
\long\def\@client{use \texttt{$\backslash$client\{\ldots\}}}
\long\def\@datei{\tiny use \texttt{$\backslash$datei\{\ldots\}}}
\long\def\@dateii{}
\long\def\@dateiii{}
\long\def\@dateiv{}
\long\def\@keywords{use \texttt{$\backslash$keywords\{\ldots\}}}
\long\def\@projectnumber{use \texttt{$\backslash$projectnumber\{\ldots\}}}
\long\def\@reference{use \texttt{$\backslash$reference\{\ldots\}}}
\long\def\@references{use \texttt{$\backslash$references\{\ldots\}}}
\long\def\@revieweri{\tiny use \texttt{$\backslash$revieweri\{\ldots\}}}
\long\def\@reviewerii{}
\long\def\@revieweriii{}
\long\def\@revieweriv{}
\long\def\@status{use \texttt{$\backslash$status\{\ldots\}}}
\long\def\@summary{use \texttt{$\backslash$summary\{\ldots\}}}
\long\def\@versioni{\tiny use \texttt{$\backslash$versioni\{\ldots\}}}
\long\def\@versionii{}
\long\def\@versioniii{}
\long\def\@versioniv{}
%
%------------------------------------------------------------------------------
%
%  TITLE PAGE
%
%------------------------------------------------------------------------------
%
\newcommand{\helv}{%
\fontfamily{phv}\selectfont}

\RequirePackage{transparent}
\definecolor{dbluelink}{RGB}{000,000,100} % the colour of hyper links in the pdf-document (daaark blue)
\definecolor{dblue}{RGB}{000,139,191} % the colour of deltares-logo.pdf
\definecolor{dname}{RGB}{150,183,197} % the colour of deltares-logo.pdf
\definecolor{dmetallic}{RGB}{128,169,189}
\definecolor{dsand}{RGB}{160,158,141}
\definecolor{dgreen}{RGB}{180,211,0}
\definecolor{dorange}{RGB}{255,114,0}
\definecolor{dcyclaam}{RGB}{235,0,115}
\definecolor{dcyan}{RGB}{0,196,214}
\definecolor{dgrey1}{RGB}{150,150,150}
\definecolor{dgrey2}{RGB}{180,180,180}
\definecolor{dgrey3}{RGB}{210,210,210}
\definecolor{dgrey4}{RGB}{232,232,232}


\def\maketitle{
    \newpage
    \thispagestyle{empty}
    \setcounter{page}{1}

    \hspace*{80pt}{\LARGE\bf \draftname}

    \vskip 200pt plus \z@ minus \z@
    \hangindent 80pt \hspace*{80pt}{\helv \HUGE\bf \helv \@title \par}
    \vfil\vfil
    \hangindent 80pt \hspace*{80pt}{\baselineskip 14pt \Large\bf \helv \@subtitle \par}
    \vfil
    \hangindent 80pt \hspace*{80pt}{\Large \bf \helv \@manualtype \par}
    \vfil\vfil\vfil
    \hangindent 80pt \hspace*{80pt}{\Large \bf \helv \@distribution \par}
    \vfil\vfil\vfil
    \hangindent 80pt \hspace*{80pt}{\helv \@author \par}
    \vfil
\if@dutch
    \hangindent 80pt \hspace*{80pt}{\helv Versie: \@version\newline
                                          SVN Revisie: \svnrev \par}
\else
    \if@spanish
        \hangindent 80pt \hspace*{80pt}{\helv Versi\'on: \@version\newline
                                              SVN Revisi\'on: \svnrev \par}
    \else
        \hangindent 80pt \hspace*{80pt}{\helv Version: \@version\newline
                                              SVN Revision: \svnrev \par}
    \fi
\fi
    \hangindent 80pt \hspace*{80pt}{\helv \@date  \par}
    \vfil

    \clearpage
    \let\maketitle\relax
}
%
%------------------------------------------------------------------------------
% Signature GB
%------------------------------------------------------------------------------
%
\newcommand{\signaturegb}{
    {\footnotesize\textbf{Title}}  \newline
    \@title

    \begin{tabular}{@{}p{39.5mm}p{39.5mm}p{41.0mm}p{\textwidth-120.0mm-42pt}}
        \footnotesize{\textbf{Client}} & \footnotesize{\textbf{Project}} &
        \footnotesize{\textbf{Reference}} & \footnotesize{\textbf{Pages}} \\
        \@client     & \@projectnumber     & \@reference     & \pageref*{LastPage}
    \end{tabular}

    {\footnotesize\textbf{Classification}} \newline
    \@classification

    {\footnotesize\textbf{Keywords}} \newline
    \@keywords

    {\footnotesize\textbf{Summary}} \newline
    \@summary

    {\footnotesize\textbf{References}} \newline
    \@references
    
    \vspace{3\baselineskip}
        {{\footnotesize
          %\begin{tabular}{p{18mm}@{}p{20mm}@{}p{25mm}@{}p{12mm}@{}p{25mm}@{}p{12mm}@{}p{25mm}@{}p{12mm}} \hline
           \begin{tabular}{p{15mm-2pt}@{}p{20mm-2pt}@{}p{25mm-2pt}@{}p{12mm}@{}p{25mm-2pt}@{}p{12mm}@{}p{25mm-2pt}@{}p{\textwidth-133mm-1pt}} \hline
            \rowcolor{dgrey4}\rule{0in}{2ex}\!\!
             \textbf{Version} & \textbf{Date} & \textbf{Author} & \textbf{Initials} & \textbf{Review} & \textbf{Initials} & \textbf{Approval} & \textbf{Initials}  \\ \hline
                \@versioni & \@datei & \@authori & & \@revieweri && \@approvali &   \\ \hline
                \@versionii & \@dateii & \@authorii & & \@reviewerii && \@approvalii &   \\ \hline
                \@versioniii & \@dateiii & \@authoriii & & \@revieweriii && \@approvaliii &   \\ \hline
                \@versioniv  & \@dateiv  & \@authoriv  & & \@revieweriv  && \@approvaliv  &   \\ \hline
          \end{tabular}
        }}

    \vfill
    {\footnotesize\textbf{Status}}  \newline
    \@status
}
%
%------------------------------------------------------------------------------
% Signature NL
%------------------------------------------------------------------------------
%
\newcommand{\signaturenl}{
    {\footnotesize\textbf{Titel}} \newline
    \@title

    \begin{tabular}{@{}p{39.5mm}p{39.5mm}p{41.0mm}p{\textwidth-120.0mm-42pt}}
        \footnotesize{\textbf{Opdrachtgever}} & \footnotesize{\textbf{Project}} &
        \footnotesize{\textbf{Kenmerk}} & \footnotesize{\textbf{Pagina's}} \\
        \@client     & \@projectnumber     & \@reference     & \pageref*{LastPage}
    \end{tabular}

    {\footnotesize\textbf{Classificatie}} \newline
    \@classification

    {\footnotesize\textbf{Trefwoorden}} \newline
    \@keywords

    {\footnotesize\textbf{Samenvatting}} \newline
    \@summary

    {\footnotesize\textbf{Referenties}}  \newline
    \@references

    \vspace{3\baselineskip}
    {{\footnotesize
      \begin{tabular}{p{14mm-2pt}@{}p{20mm-2pt}@{}p{25mm-2pt}@{}p{12mm}@{}p{25mm-2pt}@{}p{12mm}@{}p{25mm-2pt}@{}p{\textwidth-133mm-1pt}} \hline
      \rowcolor{dgrey4} \rule{0in}{2ex}\!\!
       \textbf{Versie} & \textbf{Datum} & \textbf{Auteur} & \textbf{Paraaf} & \textbf{Review} & \textbf{Paraaf} & \textbf{Goedkeuring} & \textbf{Paraaf}  \\ \hline
          \@versioni   & \@datei   & \@authori   & & \@revieweri   && \@approvali   &   \\ \hline
          \@versionii  & \@dateii  & \@authorii  & & \@reviewerii  && \@approvalii  &   \\ \hline
          \@versioniii & \@dateiii & \@authoriii & & \@revieweriii && \@approvaliii &   \\ \hline
          \@versioniv  & \@dateiv  & \@authoriv  & & \@revieweriv  && \@approvaliv  &   \\ \hline
    \end{tabular}
    }}

    \vfill
    {\footnotesize\textbf{Status}}  \newline
    \@status
}
%
%------------------------------------------------------------------------------
% Copyright GB, NL, ESP
%------------------------------------------------------------------------------
%
\newcommand{\copyrightGB}{
    Copyright {\copyright} \number\year~\Deltares\\
    All rights reserved. No part of this document may be reproduced in any form by print,
    photo print, photo copy, microfilm or any other means, without written permission from the
    publisher: Deltares.
}
\newcommand{\copyrightNL}{
    Copyright {\copyright} \number\year~\Deltares\\
    Alle rechten voorbehouden. Niets uit deze uitgave mag worden verveelvoudigd in enige vorm door middel van druk, fotokopie, 
    microfilm of op welke andere wijze dan ook, zonder voorafgaande schriftelijke toestemming van de uitgever: Deltares.
}
\newcommand{\copyrightESP}{
    Copyright {\copyright} \number\year~\Deltares\\
    Todos los derechos reservados. Queda prohibida la reproducci\'on total o parcial,
    ya sea de forma impresa, reprografiada, microfilmada o de cualquier otra forma,
    salvo autorizaci\'on previa expresa por escrito del editor: Deltares.
}
%
%------------------------------------------------------------------------------
%
%  STATE PAGE
%
%------------------------------------------------------------------------------
%
\long\def\LastPage{%
    \cleardoublepage%
    \pagestyle{empty}%
    \newpage%
    \phantom{m}%
    \newpage%
    \phantom{m}%
    %\AddToShipoutPicture*{\BottomLastPage}%
}

\if@dutch
    \newcommand{\contactsalesandsupport}{
        \begin{tabular}[t]{@{}p{0.50\textwidth}p{0.50\textwidth}}
        \begin{tabbing}
        \textbf{Verkoop:} \\
        telefoon: \= +31\,88\,335\,81\,88 \\
        fax:       \> +31\,88\,335\,81\,11 \\
        e-mail:    \> software@deltares.nl \\
        www:       \> https://www.deltares.nl/software 
        \end{tabbing}
        %
        &
        %
        \begin{tabbing}
        \textbf{Ondersteuning:} \\
        telefoon: \= +31\,88\,335\,81\,00 \\
        fax:       \> +31\,88\,335\,81\,11 \\
        e-mail:    \> software.support@deltares.nl \\
        www:       \> https://www.deltares.nl/software
        \end{tabbing}
        \end{tabular}
    }
\else
    \newcommand{\contactsalesandsupport}{
        \begin{tabular}[t]{@{}p{0.50\textwidth}p{0.50\textwidth}}
        \begin{tabbing}
        \textbf{For sales contact:} \\
        telephone: \= +31\,88\,335\,81\,88 \\
        fax:       \> +31\,88\,335\,81\,11 \\
        e-mail:    \> software@deltares.nl\\
        www:       \> https://www.deltares.nl/software
        \end{tabbing}
        %
        &
        %
        \begin{tabbing}
        \textbf{For support contact:} \\
        telephone: \= +31\,88\,335\,81\,00 \\
        fax:       \> +31\,88\,335\,81\,11 \\
        e-mail:    \> software.support@deltares.nl\\
        www:       \> https://www.deltares.nl/software
        \end{tabbing}
        \end{tabular}
    }
\fi
\newcommand{\manualtitle}{%
    \frontmatter
    \pagenumbering{roman}%
    \maketitle%
    \pagestyle{empty}%
    \textbf{\@title, \@manualtype}

    \vskip 1cm
    \vfill\vfill\vfill
    \if@dutch
        \begin{tabular}{@{}p{0.50\textwidth}p{0.50\textwidth}}
            \begin{tabbing}
                \textbf{Gepubliceerd en gedrukt door:} \\
                Deltares\\
                Boussinesqweg 1\\
                2629 HV Delft\\
                Postbus 177\\
                2600 MH Delft\\
                Nederland
            \end{tabbing}
            &
            \begin{tabbing}
            \\
            telefoon:  \=  +31\,88\,335\,82\,73\\
            fax:       \>  +31\,88\,335\,85\,82\\
            e-mail:    \>  info@deltares.nl\\
            www:       \>  https://www.deltares.nl
            \end{tabbing}
        \end{tabular}
    \else
        \if@spanish
            \begin{tabular}{@{}p{0.50\textwidth}p{0.50\textwidth}}
                \begin{tabbing}
                    \textbf{Publicado e impreso por:} \\
                    Deltares\\
                    Boussinesqweg 1\\
                    2629 HV Delft\\
                    buz\'on 177\\
                    2600 MH Delft\\
                    Pa\'ises Bajos
                \end{tabbing}
                &
                \begin{tabbing}
                    \\
                    tel\'efono: \=  +31\,88\,335\,82\,73\\
                    fax:       \>  +31\,88\,335\,85\,82\\
                    e-mail:    \>  info@deltares.nl\\
                    www:       \>  https://www.deltares.nl
                \end{tabbing}
            \end{tabular}
        \else
            \begin{tabular}{@{}p{0.50\textwidth}p{0.50\textwidth}}
                \begin{tabbing}
                    \textbf{Published and printed by:} \\
                    Deltares\\
                    Boussinesqweg 1\\
                    2629 HV Delft\\
                    P.O.\ 177\\
                    2600 MH Delft\\
                    The Netherlands
                \end{tabbing}
                &
                \begin{tabbing}
                    \\
                    telephone: \=  +31\,88\,335\,82\,73\\
                    fax:       \>  +31\,88\,335\,85\,82\\
                    e-mail:    \>  info@deltares.nl\\
                    www:       \>  https://www.deltares.nl
                \end{tabbing}
            \end{tabular}
        \fi
    \fi

    \contactsalesandsupport
    %\vskip 10mm
    %The screen pictures shown in this manual were produced during the development stages and may differ from the actual product

    \if@dutch
        \copyrightNL
    \else
        \if@spanish
            \copyrightESP
        \else
            \copyrightGB
        \fi
    \fi
%
%    Print the signatures
%
    \cleardoublepage%
    \if@signature
        \if@dutch
            \signaturenl
       \else
            \signaturegb
        \fi% end if@dutch
    \else
    \fi% end if@signature

    \pagestyle{frontmatterplain}%
    \tableofcontents%
    \conditionalLoF% Conditionally insert List of Figures
    \conditionalLoT% Conditionally insert List of Tables
    \printnomenclature% Conditionally insert List of Symbols, package: nomentbl
    \mainmatter

    % At this point document title is known for sure.
    % Now measure the width of title in page footers (later used in fancy page styles).
    % When larger than max width, reduce it to max width.
    % Result: parbox in footer has exactly right with, such that alignment and
    % possible line break of title in footer is fully automatic.
    \setlength{\maxfoottitlewidth}{125mm}
    \settowidth{\foottitlewidth}{\headerfont \@title}
    \ifthenelse{\lengthtest{\foottitlewidth > \maxfoottitlewidth}}{%
        \setlength{\foottitlewidth}{\maxfoottitlewidth}}

    %\global\let\manualtitle\relax%
}%
%------------------------------------------------------------------------------
%
% Change header and footer
%
\usepackage{fancyhdr}
\newcommand{\headerfont}{\selectfont}

\pagestyle{fancy}
\renewcommand{\chaptermark}[1]{\markboth{#1}{}}

% parbox width for title in footer: lengths are set in \manualtitle.
\newlength{\maxfoottitlewidth}
\newlength{\foottitlewidth}
% Used for title page
%\fancypagestyle{empty}{%
    %\fancyhead{}%
    %\fancyfoot{}%
    %\fancyfoot[CO,CE]{\headerfont pagestyle: empty}
%}

% Used for toc, lof, lot, preface (frontmatter)
\fancypagestyle{frontmatterplain}{%
    \fancyhead{}%
    \fancyhead[LO]{\headerfont }
    \fancyhead[CO]{\headerfont }
    \fancyhead[RO]{\headerfont \thechap} % chaptername
    \fancyhead[RE]{\headerfont }
    \fancyhead[CE]{\headerfont }
    \fancyhead[LE]{\headerfont \@title, \@manualtype}
    \fancyfoot{}
    \fancyfoot[LO,RE]{\headerfont \Deltares}
%    \fancyfoot[CO,CE]{\headerfont pagestyle: frontmatterplain}
    \fancyfoot[RO,LE]{\headerfont \thepage}
}

% Used for other pages
\fancypagestyle{headings}{%
    \fancyhead{}
    \fancyhead[LO]{\headerfont }
    \fancyhead[CO]{\headerfont }
    \fancyhead[RO]{\headerfont \thechap} % chaptername
    \fancyhead[RE]{\headerfont }
    \fancyhead[CE]{\headerfont }
    \fancyhead[LE]{\headerfont \@title, \@manualtype}
    \fancyfoot{}
    \fancyfoot[LO,RE]{\headerfont \Deltares}
%    \fancyfoot[CO,CE]{\headerfont pagestyle: headings}
    \if@dutch
        \fancyfoot[RO,LE]{\headerfont \thepage\ van \pageref*{LastPage}}%
    \else
        \fancyfoot[RO,LE]{\headerfont \thepage\ of \pageref*{LastPage}}%
    \fi
}

% Used for first page of chapter
\fancypagestyle{plainChapter}{%
    \renewcommand{\headrulewidth}{0pt}
    \fancyhead{}%
    \fancyfoot{}%
    \fancyfoot[LO,RE]{\headerfont \Deltares}
%    \fancyfoot[CO,CE]{\headerfont \thepage}%pagestyle: plainChapter}
    \if@dutch
        \fancyfoot[RO,LE]{\headerfont \thepage\ van \pageref*{LastPage}}%
    \else
        \fancyfoot[RO,LE]{\headerfont \thepage\ of \pageref*{LastPage}}%
    \fi
}

% Used for first page of part
\fancypagestyle{plainPart}{%
    \AddToShipoutPicture*{\BackgroundPicPart}
    \renewcommand{\headrulewidth}{0pt}
    \fancyhead{}%
    \fancyfoot{}%
%    \fancyfoot[LO,RE]{\headerfont \Deltares}
%    \fancyfoot[CO,CE]{\headerfont pagestyle: plainPart}
%    \fancyfoot[RO,LE]{\headerfont \thepage}%
}

\addtocounter{secnumdepth}{1}
\setlength{\parskip}{\baselineskip}
\setlength{\parindent}{0pt}
\setlength{\parsep}{1pt}
\setlength{\partopsep}{1pt}
%\setlength{\rightskip}{0pt plus 2cm} % Niet rechts uitlijnen tekst
%\setlength{\@rightskip}{0pt plus 2cm} % Niet rechts uitlijnen lists

%------------------------------------------------------------------------------
%
% To change the caption of tables and figures.
% Format is defautl or `hang' wich aligns caption to the right of the figure statement
\captionsetup{format=hang,margin=0.9cm,font=small,labelfont={sl,bf},textfont=sl}
%
%------------------------------------------------------------------------------
%
% Define a more advanced remark environment that writes Remark if there is only
% one item and that writes Remarks if there are multiple items.
%
\usepackage{ifthen}

%
\newcounter{n_sp_blocks}
\newcounter{n_sp_items}
\newenvironment{singleplural}[2]%
  {\addtocounter{n_sp_blocks}{1}%
   \ifthenelse{\equal{a\ref{sp_block\arabic{n_sp_blocks}}}{a1}}%
     {\textbf{#1:}}%
     {\textbf{#2:}}%
   \begin{list}{\labelitemi}{%
        \setlength{\parskip}{0pt}%
        \setlength{\parsep}{0pt}%
        \setlength{\topsep}{0pt}%
        \setlength{\itemsep}{0pt}%
        \usecounter{n_sp_items}%
   }}%
  {\label{sp_block\arabic{n_sp_blocks}}\end{list}}
%
\if@dutch
    \newenvironment{Remark}%
        {\marginnote{{\includegraphics[height=1.5\baselineskip]{pictures/margin_remark.pdf}}}[1.50ex]%
         \begin{singleplural}{Opmerking}{Opmerkingen}}
        {\end{singleplural}}
\else
    \newenvironment{Remark}%
        {\marginnote{{\includegraphics[height=1.5\baselineskip]{pictures/margin_remark.pdf}}}[1.50ex]%
         \begin{singleplural}{Remark}{Remarks}}
        {\end{singleplural}}
\fi
%
\if@dutch
    \newenvironment{Restriction}%
        {\marginnote{{\includegraphics[height=1.5\baselineskip]{pictures/margin_restriction.pdf}}}[1.50ex]%
         \begin{singleplural}{Restrictie}{Restricties}}
        {\end{singleplural}}
\else
    \newenvironment{Restriction}%
        {\marginnote{{\includegraphics[height=1.5\baselineskip]{pictures/margin_restriction.pdf}}}[1.50ex]%
         \begin{singleplural}{Restriction}{Restrictions}}
        {\end{singleplural}}
\fi
%
\if@dutch
    \newenvironment{Warning}%
        {\marginnote{{\includegraphics[height=1.5\baselineskip]{pictures/margin_warning.pdf}}}[1.25ex]%
         \begin{singleplural}{Waarschuwing}{Waarschuwingen}}
        {\end{singleplural}}
\else
    \newenvironment{Warning}%
        {\marginnote{{\includegraphics[height=1.5\baselineskip]{pictures/margin_warning.pdf}}}[1.25ex]%
         \begin{singleplural}{Warning}{Warnings}}
        {\end{singleplural}}
\fi
%
\newcommand{\Note}%
{%
    \textbf{Note:\ }
    \marginnote{ {\includegraphics[height=1.5\baselineskip]{pictures/margin_note.pdf}} }[-0.75ex]%
}
%
\newcommand{\Tip}%
{%
    \textbf{Tip:\ }%
    \marginnote[ {\includegraphics[height=1.5\baselineskip]{pictures/margin_tip_left.pdf}} ]%
               { {\includegraphics[height=1.5\baselineskip]{pictures/margin_tip_right.pdf}} }[-0.75ex]%
}
% todo macro: usage of package marginnote
\renewcommand*{\raggedleftmarginnote}{\raggedleft}
\renewcommand*{\raggedrightmarginnote}{\raggedright}
\newcommand{\todo}[1]{
    {\color[rgb]{0,0.2,0.5}\textbf{TODO(??):} }%
    \protect\marginnote{\textbf{TODO}}%
    {\color[rgb]{0,0.2,0.5}\slshape #1}%
}
%
\RecustomVerbatimEnvironment
  {Verbatim}{Verbatim}
  {fontsize=\footnotesize}

%
% Nomenclature table
\if@dutch
    \def\nomname{Lijst van symbolen}
    \def\nomAname{Latijnse letters}%
    \def\nomGname{Griekse letters}%
    \def\nomXname{Superscripts}%
    \def\nomZname{Subscripts}%
\else
    \def\nomname{List of Symbols}
    \def\nomAname{Latin Letters}%
    \def\nomGname{Greek Letters}%
    \def\nomXname{Superscripts}%
    \def\nomZname{Subscripts}%
\fi
%
%------------------------------------------------------------------------------
% Bring items closer together in list environments
%
\let\orig@Itemize =\itemize
\let\orig@Enumerate =\enumerate
\let\orig@Description =\description
% Zero the vertical spacing parameters
\def\Nospacing{\itemsep=0pt\topsep=0pt\partopsep=0pt\parskip=0pt\parsep=0pt}
\renewenvironment{itemize}{\orig@Itemize\Nospacing}{\endlist}
\renewenvironment{enumerate}{\orig@Enumerate\Nospacing}{\endlist}
\renewenvironment{description}{\orig@Description\Nospacing}{\endlist}

\renewcommand\labelitemi{$\boldsymbol\diamond$}
\renewcommand\labelitemii{\rotatebox{45}{$\boldsymbol\diamond$}} %
\renewcommand\labelitemiii{$\boldsymbol\circ$}
\renewcommand\labelitemiv{$\boldsymbol\triangleright$}

\renewcommand{\labelenumi}{\arabic{enumi}}
\renewcommand{\labelenumii}{\arabic{enumi}.\arabic{enumii}}
\renewcommand{\labelenumiii}{\arabic{enumi}.\arabic{enumii}.\arabic{enumiii}}
\renewcommand{\labelenumiv}{\arabic{enumi}.\arabic{enumii}.\arabic{enumiii}.\arabic{enumiv}}

%
%------------------------------------------------------------------------------

\newenvironment{guilist}%
  {\begin{list}{}%
    {%
      \renewcommand{\makelabel}[1]%
        {\raisebox{0pt}[1ex][0pt]{\makebox[\labelwidth][l]%
            {\parbox[t]{\labelwidth}{\raggedright\hspace{0pt}##1}}}}%
        \setlength{\parskip}{\baselineskip/2}%
        \setlength{\labelwidth}{3.7cm}%
        \setlength{\leftmargin}{\labelwidth+\labelsep}%
        \setlength{\topsep}{0pt}%
        \setlength{\topsep}{0pt}%
        \setlength{\parsep}{0pt}%
        \setlength{\itemsep}{0pt}%
    }%
  }%
  {\end{list}}

\if@dutch
    \newenvironment{domain} % Added more spacing around the text in the header, 19/11/2008 RMC Spiteri
    {%
        \small\begin{tabular}{|p{4.4cm}|p{2.2cm}|p{2.2cm}|p{2.2cm}|p{\textwidth-63pt-11.0cm}|}
        \multicolumn{5}{@{}l}{\STRUT \textbf{Domein:}} \\ [1ex] \hline  \STRUT
        \textbf{Parameter} & \textbf{Onder limiet} & \textbf{Boven limiet} & \textbf{Default} & \textbf{Eenheid} \\[1ex] \hline
    }
    {\end{tabular}}
\else
    \newenvironment{domain} % Added more spacing around the text in the header, 19/11/2008 RMC Spiteri
    {%
        \small\begin{tabular}{|p{4.4cm}|p{2.2cm}|p{2.2cm}|p{2.2cm}|p{\textwidth-63pt-11.0cm}|}
        \multicolumn{5}{@{}l}{\STRUT \textbf{Domain:}} \\ [1ex] \hline  \STRUT
        \textbf{Parameter} & \textbf{Lower limit} & \textbf{Upper limit} & \textbf{Default} & \textbf{Unit} \\[1ex] \hline
    }
    {\end{tabular}}
\fi

\newenvironment{action}%
  {\begin{list}{}{%
    \renewcommand{\makelabel}[1]{\includegraphics[height=0.5\baselineskip]{pictures/action_arrow.pdf}\hfill\varlistcolon}%
    \setlength{\labelwidth}{1em+1pt}% Changed to 3.7cm, 19/11/2008 RMC Spiteri
    \setlength{\leftmargin}{\labelwidth+\labelsep}%
    \setlength{\parskip}{\baselineskip/2}%
    \setlength{\topsep}{0pt}%
    \setlength{\parsep}{0pt}%
    \setlength{\itemsep}{0pt}}}
  {\end{list}}
%
%                             definition of varlist environment for listing
%                               variables and their meaning after a formula
%
\newcommand{\varlistcolon}{} % originally :
\newenvironment{varlist}%
  {\begin{list}{}{%
    \renewcommand{\makelabel}[1]{\quad{##1}\hfill\varlistcolon}%
    \setlength{\labelwidth}{3.7cm}% Changed to 3.7cm, 19/11/2008 RMC Spiteri
    \setlength{\leftmargin}{\labelwidth+\labelsep}%
    \setlength{\parskip}{\baselineskip/2}%
    \setlength{\topsep}{0pt}%
    \setlength{\parsep}{0pt}%
    \setlength{\itemsep}{0pt}}}
  {\end{list}}

%------------------------------------------------------------------------------

\newenvironment{symbollist}%
  {\begin{list}{}{%
    \renewcommand{\makelabel}[1]{\quad{##1}\hfill\varlistcolon}%
    \setlength{\labelwidth}{2.0cm}
    \setlength{\leftmargin}{\labelwidth+\labelsep}%
    \setlength{\parskip}{\baselineskip/2}%
    \setlength{\topsep}{0pt}%
    \setlength{\parsep}{0pt}%
    \setlength{\itemsep}{0pt}}}
  {\end{list}}

%------------------------------------------------------------------------------

\newenvironment{options}%
  {\underline{Options:}
  \begin{list}{}%
    {%
      \renewcommand{\makelabel}[1]%
        {\raisebox{0pt}[1ex][0pt]{\makebox[\labelwidth][l]%
            {\parbox[t]{\labelwidth}{\raggedright\hspace{0pt}##1}}}}%
        \setlength{\labelwidth}{3.7cm}%
        \setlength{\leftmargin}{\labelwidth+\labelsep}%
        \setlength{\parskip}{\baselineskip/2}%
        \setlength{\topsep}{0pt}%
        \setlength{\parsep}{0pt}%
        \setlength{\itemsep}{0pt}
    }%
  }%
  {\end{list}}

%------------------------------------------------------------------------------

\def\nonumchapter#1{%
    \chapter*{#1}%
    \addcontentsline{toc}{chapter}{#1}%
    \markboth{#1}{}}

%------------------------------------------------------------------------------
\newcommand{\clearemptydoublepage}{\newpage{\pagestyle{empty}\cleardoublepage}}
\newcommand{\STRUT}{\rule{0in}{3ex}}


\newcommand{\WL}{{WL\,$|$\,Delft Hydraulics}\xspace}
\newcommand{\Deltares}{Deltares\xspace}
\newcommand{\degr}{\ensuremath{^\circ}}
\newcommand{\tenlog}{ \ensuremath{\ {}^{10}\!\log} }
\newcommand{\twolog}{ \ensuremath{\ {}^{2}\!\log} }
\newcommand{\xml}{xml\xspace}

\newcommand{\window}[1]{\textbf{#1}}
\newcommand{\menu}[1]{\emph{#1}}
\newcommand{\button}[1]{\emph{#1}}
\newcommand{\file}[1]{$<${#1}$>$} 
\newcommand{\dir}[1]{$<${#1}$>$}  
\newcommand{\ginput}[1]{``#1''}
\newcommand{\command}[1]{{\normalfont\ttfamily #1}}
\newcommand{\block}[1]{\textbf{\rmfamily #1}}
\newcommand{\keyw}[1]{{\normalfont\ttfamily #1}}
\newcommand{\unitbrackets}[1]{[#1]}
\newcommand{\ext}[1]{$<$#1$>$}
\newcommand{\key}[1]{\texttt{#1}}
\newcommand{\function}[1]{\emph{#1}}

\DeclareRobustCommand\dotp{\mathpalette\dotp@{.5}} % bullet for the dot product
\newcommand\dotp@[2]{\mathbin{\vcenter{\hbox{\scalebox{#2}{$\m@th#1\bullet$}}}}}
\newcommand{\ssquare}{\stackrel{\scriptscriptstyle\square}{{}}} % small square for itemize list

\newcommand{\option}[1]{\vspace{0.5cm}\underline{#1}}  % Emphasize an option, mainly used in WAQ TRM

%% Referencing
\let\origtheequation\theequation
\def\tagform@#1{\maketag@@@{\ignorespaces#1\unskip\@@italiccorr}}
\renewcommand{\theequation}{(\origtheequation)}
\newcommand*\Autoref[1]{% First letter is a capital, also from macro \Autoref
    \begingroup
        \if@dutch
            \def\appendixautorefname{Appendix}%
            \def\chapterautorefname{Hoofdstuk}%
            \def\equationautorefname{Vergelijking}%
            \def\figureautorefname{Figuur}%
            \def\footnoteautorefname{Voetnoot}%
            \def\itemautorefname{Item}%
            \def\pageautorefname{Bladzijde}%
            \def\paragraphautorefname{Paragraaf}%
            \def\partautorefname{Deel}%
            \def\sectionautorefname{Paragraaf}%
            \def\subparagraphautorefname{Paragraaf}%
            \def\subsectionautorefname{Paragraaf}%
            \def\subsubsectionautorefname{Paragraaf}%
            \def\tableautorefname{Tabel}%
            \def\theoremautorefname{Theorema}%
        \else
            \def\appendixautorefname{Appendix}%
            \def\chapterautorefname{Chapter}%
            \def\equationautorefname{Equation}%
            \def\figureautorefname{Figure}%
            \def\footnoteautorefname{Footnote}%
            \def\itemautorefname{Item}%
            \def\pageautorefname{Page}%
            \def\paragraphautorefname{Section}%
            \def\partautorefname{Part}%
            \def\sectionautorefname{Section}%
            \def\subparagraphautorefname{Section}%
            \def\subsectionautorefname{Section}%
            \def\subsubsectionautorefname{Section}%
            \def\tableautorefname{Table}%
            \def\theoremautorefname{Theorem}%
        \fi
        % ...
        \autoref{#1}%
    \endgroup%
}

\if@dutch % First letter is a lower case, also from macro \autoref
    \def\appendixautorefname{appendix}%
    \def\chapterautorefname{hoofdstuk}%
    \def\equationautorefname{vergelijking}%
    \def\figureautorefname{figuur}%
    \def\footnoteautorefname{voetnoot}%
    \def\itemautorefname{item}%
    \def\pageautorefname{bladzijde}%
    \def\paragraphautorefname{paragraaf}%
    \def\partautorefname{deel}%
    \def\sectionautorefname{paragraaf}%
    \def\subparagraphautorefname{paragraaf}%
    \def\subsectionautorefname{paragraaf}%
    \def\subsubsectionautorefname{paragraaf}%
    \def\tableautorefname{tabel}%
    \def\theoremautorefname{theorema}%
\else
    \addto\extrasenglish{
        \def\appendixautorefname{Appendix}%
        \def\chapterautorefname{chapter}%
        \def\equationautorefname{Equation}%
        \def\figureautorefname{Figure}%
        \def\footnoteautorefname{footnote}%
        \def\itemautorefname{item}%
        \def\pageautorefname{page}%
        \def\paragraphautorefname{section}%
        \def\partautorefname{Part}%
        \def\sectionautorefname{section}%
        \def\subparagraphautorefname{section}%
        \def\subsectionautorefname{section}%
        \def\subsubsectionautorefname{section}%
        \def\tableautorefname{Table}%
        \def\theoremautorefname{Theorem}%
    }
\fi

\newcommand{\Cref}[1]{\autoref{#1}}
\newcommand{\Sref}[1]{\autoref{#1}}
\newcommand{\Fref}[1]{\autoref{#1}}
\if@dutch
    \newcommand{\Frefs}[2]{Figuren~\ref{#1} en \ref{#2}}
    \newcommand{\Freft}[2]{Figuren~\ref{#1} tot \ref{#2}}
\else
    \newcommand{\Frefs}[2]{Figures~\ref{#1} and \ref{#2}}
    \newcommand{\Freft}[2]{Figures~\ref{#1} to \ref{#2}}
\fi
\newcommand{\Aref}[1]{\autoref{#1}}
\newcommand{\Apref}[1]{\autoref{#1}}
\newcommand{\Tref}[1]{\autoref{#1}}
\newcommand{\Eref}[1]{\autoref{#1}}


\hyphenation{bathy-metry fa-vour-able}

\def\execute{%
    \begingroup%
        \catcode`\%=12%
        \catcode`\\=12%
        \executeaux}%
        \def\executeaux#1{\immediate\write18{#1}%
    \endgroup%
}
%------------------------------------------------------------------------------
%
% PART
%
% Copied from report.cls, but pagestyle now depends on front/mainmatter:
\newcommand\BackgroundPicPart{}
\usepackage{eso-pic}
\renewcommand\chapter{\if@openright\cleardoublepage\else\clearpage\fi
                    \if@mainmatter\thispagestyle{plainChapter}\else\thispagestyle{frontmatterplain}\fi%
                    \global\@topnum\z@
                    \@afterindentfalse
                    \secdef\@chapter\@schapter
                    }
\def\@part[#1]#2{%
    \ifnum \c@secnumdepth >-2\relax
      \refstepcounter{part}%
      \addcontentsline{toc}{part}{\thepart\hspace{1em}#1}%
    \else
      \addcontentsline{toc}{part}{#1}%
    \fi
    \markboth{}{}%
    {\centering
     \interlinepenalty \@M
     \normalfont
     \thispagestyle{plainPart}
     \ifnum \c@secnumdepth >-2\relax
       \leavevmode\HUGE\bfseries\raggedright #2
       \par
       \vskip 20\p@
     \fi
     \par}%
     \renewcommand{\partname}[1]{#1}
    \@endpart}
%------------------------------------------------------------------------------
%
% CHAPTER redefinition Hoofdstuk nummmer en titel
%
\newcommand\BackgroundPicChapter{}
\def\@makechapterhead#1{%
    \AddToShipoutPicture*{\BackgroundPicChapter}  %
    \vspace{0\baselineskip}%
    {\parindent \z@%
        \leavevmode\Large\bfseries\raggedright%
          \llap{\thechapter\hspace{4.5mm}}#1%
        \par\nobreak
        \vskip 0\baselineskip%
    }
}
\def\@makeschapterhead#1{%
    \AddToShipoutPicture*{\BackgroundPicChapter}  %
    \vspace{0\baselineskip}%
    {\parindent \z@
        \leavevmode\Large\bfseries\raggedright%
        \llap{}#1
        \par\nobreak
        \vskip 0\baselineskip%
    }
}
%------------------------------------------------------------------------------
%
% SECTION redefinition
%
\def\section{\@ifstar\unnumberedsection\numberedsection}
\def\numberedsection{\@ifnextchar[%]
  \numberedsectionwithtwoarguments\numberedsectionwithoneargument}
\def\unnumberedsection{\@ifnextchar[%]
  \unnumberedsectionwithtwoarguments\unnumberedsectionwithoneargument}
\def\numberedsectionwithoneargument#1{\numberedsectionwithtwoarguments[#1]{#1}}
\def\unnumberedsectionwithoneargument#1{\unnumberedsectionwithtwoarguments[#1]{#1}}

\def\numberedsectionwithtwoarguments[#1]#2{%
  \def\@currentlabelname{#1}
  \ifhmode\par\fi
  \removelastskip
  \vskip -0.5\baselineskip\nopagebreak
  \refstepcounter{section}%
  \noindent
  \begingroup
  % Taken from package needspaces
  \setlength{\dimen@}{3\baselineskip}%
  \vskip\z@\@plus\dimen@
  \penalty -100\vskip\z@\@plus -\dimen@
  \vskip\dimen@
  \penalty 9999%
  \vskip -\dimen@
  \vskip\z@skip % hide the previous \vskip from \addvspace
  % end needspaces
  \leavevmode\normalsize\bfseries\raggedright
  \llap{\thesection\hspace{4.5mm}}%
  #2
  \par
  \endgroup
  \vskip -0.5\parskip\nobreak
  \addcontentsline{toc}{section}{%
    \protect\numberline{\thesection}%
    #1}%
  }
\def\unnumberedsectionwithtwoarguments[#1]#2{%
  \def\@currentlabelname{#1}
  \ifhmode\par\fi
  \removelastskip
  \vskip -0.5\baselineskip\nopagebreak
%  \refstepcounter{section}%
  \noindent
  \begingroup
  % Taken from package needspaces
  \setlength{\dimen@}{3\baselineskip}%
  \vskip\z@\@plus\dimen@
  \penalty -100\vskip\z@\@plus -\dimen@
  \vskip\dimen@
  \penalty 9999%
  \vskip -\dimen@
  \vskip\z@skip % hide the previous \vskip from \addvspace
  % end needspaces
  \leavevmode\normalsize\bfseries\raggedright
%  \llap{\thesection\hspace{4.5mm}}%
  #2
  \par
  \endgroup
  \vskip -0.5\parskip\nobreak
%  \addcontentsline{toc}{section}{%
%     \protect\numberline{\thesection}%
%     #1}%
  }
%------------------------------------------------------------------------------
%
% SUBSECTION redefinition
%
\def\subsection{\@ifstar\unnumberedsubsection\numberedsubsection}
\def\numberedsubsection{\@ifnextchar[%]
  \numberedsubsectionwithtwoarguments\numberedsubsectionwithoneargument}
\def\unnumberedsubsection{\@ifnextchar[%]
  \unnumberedsubsectionwithtwoarguments\unnumberedsubsectionwithoneargument}
\def\numberedsubsectionwithoneargument#1{\numberedsubsectionwithtwoarguments[#1]{#1}}
\def\unnumberedsubsectionwithoneargument#1{\unnumberedsubsectionwithtwoarguments[#1]{#1}}

\def\numberedsubsectionwithtwoarguments[#1]#2{%
  \def\@currentlabelname{#1}
  \ifhmode\par\fi
  \removelastskip
  \vskip -0.5\baselineskip\nopagebreak
  \refstepcounter{subsection}%
  \noindent
  \begingroup
  % Taken from package needspaces
  \setlength{\dimen@}{3\baselineskip}%
  \vskip\z@\@plus\dimen@
  \penalty -100\vskip\z@\@plus -\dimen@
  \vskip\dimen@
  \penalty 9999%
  \vskip -\dimen@
  \vskip\z@skip % hide the previous \vskip from \addvspace
  % end needspaces
  \leavevmode\normalsize\bfseries\raggedright
  \llap{\thesubsection\hspace{4.5mm}}%
  #2
\par
  \endgroup
  \vskip -0.5\parskip\nobreak
  \addcontentsline{toc}{subsection}{%
    \protect\numberline{\thesubsection}%
    #1}%
  }
\def\unnumberedsubsectionwithtwoarguments[#1]#2{%
  \def\@currentlabelname{#1}
  \ifhmode\par\fi
  \removelastskip
  \vskip -0.5\baselineskip\nopagebreak
%  \refstepcounter{subsection}%
  \noindent
  \begingroup
  % Taken from package needspaces
  \setlength{\dimen@}{3\baselineskip}%
  \vskip\z@\@plus\dimen@
  \penalty -100\vskip\z@\@plus -\dimen@
  \vskip\dimen@
  \penalty 9999%
  \vskip -\dimen@
  \vskip\z@skip % hide the previous \vskip from \addvspace
  % end needspaces
  \leavevmode\normalsize\bfseries\raggedright
%  \llap{\thesubsection\hspace{4.5mm}}%
  #2
  \par
  \endgroup
  \vskip -0.5\parskip\nobreak
%    \addcontentsline{toc}{subsection}{%
%    \protect\numberline{\thesubsection}%
%    #1}%
  }
%------------------------------------------------------------------------------
%
% SUBSUBSECTION redefinition
%
\def\subsubsection{\@ifstar\unnumberedsubsubsection\numberedsubsubsection}
\def\numberedsubsubsection{\@ifnextchar[%]
  \numberedsubsubsectionwithtwoarguments\numberedsubsubsectionwithoneargument}
\def\unnumberedsubsubsection{\@ifnextchar[%]
  \unnumberedsubsubsectionwithtwoarguments\unnumberedsubsubsectionwithoneargument}
\def\numberedsubsubsectionwithoneargument#1{\numberedsubsubsectionwithtwoarguments[#1]{#1}}
\def\unnumberedsubsubsectionwithoneargument#1{\unnumberedsubsubsectionwithtwoarguments[#1]{#1}}

\def\numberedsubsubsectionwithtwoarguments[#1]#2{%
  \def\@currentlabelname{#1}
  \ifhmode\par\fi
  \removelastskip
  \vskip -0.5\baselineskip\nopagebreak
  \refstepcounter{subsubsection}%
  \noindent
  \begingroup
  % Taken from package needspaces
  \setlength{\dimen@}{3\baselineskip}%
  \vskip\z@\@plus\dimen@
  \penalty -100\vskip\z@\@plus -\dimen@
  \vskip\dimen@
  \penalty 9999%
  \vskip -\dimen@
  \vskip\z@skip % hide the previous \vskip from \addvspace
  % end needspaces
  \leavevmode\normalsize\bfseries\raggedright
  \llap{\thesubsubsection\hspace{4.5mm}}%
  #2
  \par
  \endgroup
  \vskip -0.5\parskip\nobreak
  \addcontentsline{toc}{subsubsection}{%
    \protect\numberline{\thesubsubsection}%
    #1}%
  }
\def\unnumberedsubsubsectionwithtwoarguments[#1]#2{%
  \def\@currentlabelname{#1}
  \ifhmode\par\fi
  \removelastskip
  \vskip -0.5\baselineskip\nopagebreak
%  \refstepcounter{subsubsection}%
  \noindent
  \begingroup
  % Taken from package needspaces
  \setlength{\dimen@}{3\baselineskip}%
  \vskip\z@\@plus\dimen@
  \penalty -100\vskip\z@\@plus -\dimen@
  \vskip\dimen@
  \penalty 9999%
  \vskip -\dimen@
  \vskip\z@skip % hide the previous \vskip from \addvspace
  % end needspaces
  \leavevmode\normalsize\bfseries\raggedright
%  \llap{\thesubsubsection\hspace{4.5mm}}%
  #2
  \par
  \endgroup
  \vskip -0.5\parskip\nobreak
%    \addcontentsline{toc}{subsubsection}{%
%    \protect\numberline{\thesubsubsection}%
%    #1}%
  }
%------------------------------------------------------------------------------
%
% PARAGRAPH redefinition
%
\def\paragraph{\@ifstar\unnumberedparagraph\numberedparagraph}
\def\numberedparagraph{\@ifnextchar[%]
  \numberedparagraphwithtwoarguments\numberedparagraphwithoneargument}
\def\unnumberedparagraph{\@ifnextchar[%]
  \unnumberedparagraphwithtwoarguments\unnumberedparagraphwithoneargument}
\def\numberedparagraphwithoneargument#1{\numberedparagraphwithtwoarguments[#1]{#1}}
\def\unnumberedparagraphwithoneargument#1{\unnumberedparagraphwithtwoarguments[#1]{#1}}

\def\numberedparagraphwithtwoarguments[#1]#2{%
  \def\@currentlabelname{#1}
  \ifhmode\par\fi
  \removelastskip
  \vskip -0.5\baselineskip\nopagebreak
  \refstepcounter{paragraph}%
  \noindent
  \begingroup
  % Taken from package needspaces
  \setlength{\dimen@}{3\baselineskip}%
  \vskip\z@\@plus\dimen@
  \penalty -100\vskip\z@\@plus -\dimen@
  \vskip\dimen@
  \penalty 9999%
  \vskip -\dimen@
  \vskip\z@skip % hide the previous \vskip from \addvspace
  % end needspaces
  \leavevmode\normalsize\bfseries\raggedright
  \llap{\theparagraph\hspace{4.5mm}}%
  #2
  \par
  \endgroup
  \vskip -0.5\parskip\nobreak
  \addcontentsline{toc}{paragraph}{%
    \protect\numberline{\theparagraph}%
    #1}%
  }
\def\unnumberedparagraphwithtwoarguments[#1]#2{%
  \def\@currentlabelname{#1}
  \ifhmode\par\fi
  \removelastskip
  \vskip -0.5\baselineskip\nopagebreak
%  \refstepcounter{paragraph}%
  \noindent
  \begingroup
  % Taken from package needspaces
  \setlength{\dimen@}{3\baselineskip}%
  \vskip\z@\@plus\dimen@
  \penalty -100\vskip\z@\@plus -\dimen@
  \vskip\dimen@
  \penalty 9999%
  \vskip -\dimen@
  \vskip\z@skip % hide the previous \vskip from \addvspace
  % end needspaces
  \leavevmode\normalsize\bfseries\raggedright
%  \llap{\theparagraph\hspace{4.5mm}}%
  #2
  \par
  \endgroup
  \vskip -0.5\parskip\nobreak
%    \addcontentsline{toc}{paragraph}{%
%    \protect\numberline{\theparagraph}%
%    #1}%
  }
%------------------------------------------------------------------------------
%
% SUBPARAGRAPH redefinition
%
\def\subparagraph{\@ifstar\unnumberedsubparagraph\numberedsubparagraph}
\def\numberedsubparagraph{\@ifnextchar[%]
  \numberedsubparagraphwithtwoarguments\numberedsubparagraphwithoneargument}
\def\unnumberedsubparagraph{\@ifnextchar[%]
  \unnumberedsubparagraphwithtwoarguments\unnumberedsubparagraphwithoneargument}
\def\numberedsubparagraphwithoneargument#1{\numberedsubparagraphwithtwoarguments[#1]{#1}}
\def\unnumberedsubparagraphwithoneargument#1{\unnumberedsubparagraphwithtwoarguments[#1]{#1}}

\def\numberedsubparagraphwithtwoarguments[#1]#2{%
  \def\@currentlabelname{#1}
  \ifhmode\par\fi
  \removelastskip
  \vskip -0.5\baselineskip\nopagebreak
  \refstepcounter{subparagraph}%
  \noindent
  \begingroup
  % Taken from package needspaces
  \setlength{\dimen@}{3\baselineskip}%
  \vskip\z@\@plus\dimen@
  \penalty -100\vskip\z@\@plus -\dimen@
  \vskip\dimen@
  \penalty 9999%
  \vskip -\dimen@
  \vskip\z@skip % hide the previous \vskip from \addvspace
  % end needspaces
  \leavevmode\normalsize\bfseries\itshape\raggedright
  \llap{\thesubparagraph\hspace{4.5mm}}%
  #2
  \par
  \endgroup
  \vskip -0.5\parskip\nobreak
  \addcontentsline{toc}{subparagraph}{%
    \protect\numberline{\thesubparagraph}%
    #1}%
  }
\def\unnumberedsubparagraphwithtwoarguments[#1]#2{%
  \def\@currentlabelname{#1}
  \ifhmode\par\fi
  \removelastskip
  \vskip -0.5\baselineskip\nopagebreak
%  \refstepcounter{subparagraph}%
  \noindent
  \begingroup
  % Taken from package needspaces
  \setlength{\dimen@}{3\baselineskip}%
  \vskip\z@\@plus\dimen@
  \penalty -100\vskip\z@\@plus -\dimen@
  \vskip\dimen@
  \penalty 9999%
  \vskip -\dimen@
  \vskip\z@skip % hide the previous \vskip from \addvspace
  % end needspaces
  \leavevmode\normalsize\bfseries\itshape\raggedright
%  \llap{\thesubparagraph\hspace{4.5mm}}%
  #2
  \par
  \endgroup
  \vskip -0.5\parskip\nobreak
%    \addcontentsline{toc}{subparagraph}{%
%    \protect\numberline{\thesubparagraph}%
%    #1}%
  }

%------------------------------------------------------------------------------
%
% INDEX
%
\usepackage{makeidx}
\usepackage[totoc=true]{idxlayout}
%------------------------------------------------------------------------------
%
%\def\promille{\leavevmode\kern-.05em
%\raise.5ex\hbox{\the\scriptfont0 0}\kern-.15em/\kern-.15em
%\lower.25ex\hbox{\the\scriptfont0 00}}

\newcommand{\promille}{%
  \relax\ifmmode\promillezeichen
        \else\leavevmode\(\mathsurround=0pt\promillezeichen\)\fi}
\newcommand{\promillezeichen}{%
  \kern-.05em%
  \raise.5ex\hbox{\the\scriptfont0 0}%
  \kern-.15em/\kern-.15em%
  \lower.25ex\hbox{\the\scriptfont0 00}}

%% Longtable improvement
% \nobreakhline is a \hline in a longtable that disallows a page break after it.
% From: http://tex.stackexchange.com/questions/6350/how-to-disable-pagebreak-on-hline-in-longtable
\def\nobreakhline{%
  \noalign{\ifnum0=`}\fi
    \penalty\@M
    \futurelet\@let@token\LT@@nobreakhline%
}
\def\LT@@nobreakhline{%
  \ifx\@let@token\hline
    \global\let\@gtempa\@gobble
    \gdef\LT@sep{\penalty\@M\vskip\doublerulesep}% <-- change here
  \else
    \global\let\@gtempa\@empty
    \gdef\LT@sep{\penalty\@M\vskip-\arrayrulewidth}% <-- change here
  \fi
  \ifnum0=`{\fi}%
  \multispan\LT@cols
     \unskip\leaders\hrule\@height\arrayrulewidth\hfill\cr
  \noalign{\LT@sep}%
  \multispan\LT@cols
     \unskip\leaders\hrule\@height\arrayrulewidth\hfill\cr
  \noalign{\penalty\@M}%
  \@gtempa
}

%------------------------------------------------------------------------------
% footnote references
%
\newcommand{\footnoteref}[1]{%
% requires a label to be present in the footnote to be referenced!
\ltx@ifpackageloaded{hyperref}{% hyperref package loaded
  \ifHy@hyperfootnotes% option hyperfootnotes=true
    \hbox{\hyperref[#1]{%
            % once upon a time this was not "#1" but "\footnote@reference",
            % but somewhere something changed at some time
            \@textsuperscript {\normalfont \ref*{#1}}}}%
  \else% option hyperfootnotes=false
    \hbox{\@textsuperscript {\normalfont \ref*{#1}}}%
  \fi%
}{% hyperref package not loaded
    \hbox{\@textsuperscript {\normalfont \ref{#1}}}%
 }%
}
%------------------------------------------------------------------------------

\newcommand{\diff}[3][1]{\ensuremath{\frac{d \ifthenelse{\equal{#1}{1}}{}{^#1} {#2}}{d {#3} \ifthenelse{\equal{#1}{1}}{}{^#1}}}}
\newcommand{\pdiff}[3][1]{\ensuremath{\frac{\partial \ifthenelse{\equal{#1}{1}}{}{^#1} {#2}}{\partial {#3} \ifthenelse{\equal{#1}{1}}{}{^#1}}}}

\newcommand{\perscomm}[1]{\citetext{#1, pers.\ comm.}}
\newcommand{\deltarestitle}{\manualtitle}%

\if@dbiblatex
    \usepackage[backend=bibtex, style=deltares_authoryear, maxcitenames=2, maxbibnames=10]{biblatex}
    \newcommand{\printrefsegment}{\printbibliography[segment=\therefsegment, heading=subbibheading]}
    \newcommand{\printallbibliography}{\printbibliography[heading=bibintoc]}
    \if@dutch
        \DefineBibliographyStrings{dutch}{%
            bibliography = {Referenties}
        }
        \defbibheading{subbibheading}{
            \subparagraph*{Referenties} % for \Autoref{refsegment:\therefsection\therefsegment}}
        }
        \DeclareFieldFormat{urldate}{%
            (Benaderd %
            \thefield{urlday} \addspace%
            \mkbibmonth{\thefield{urlmonth}}\addspace%
            \thefield{urlyear}\isdot)}
    \else %  i.e. English
        \DefineBibliographyStrings{english}{%
            bibliography = {References}
        }
        \defbibheading{subbibheading}{
            \subparagraph*{References} % for \Autoref{refsegment:\therefsection\therefsegment}}
        }
        \DeclareFieldFormat{urldate}{%
            (Accessed %
            \thefield{urlday} \addspace%
            \mkbibmonth{\thefield{urlmonth}}\addspace%
            \thefield{urlyear}\isdot)}
    \fi
    \setlength\bibitemsep{0.5\baselineskip}
\else
    \usepackage{natbib}
    \bibpunct[, ]{(}{)}{;}{a}{,}{,}
    \renewcommand\bibsection{}
    \if@dutch
        \bibliographystyle{deltares_chicago_like_dutch}
    \else
        \if@spanish
            \bibliographystyle{deltares_chicago_like_spanish}
        \else
            \bibliographystyle{deltares_chicago_like}
        \fi
    \fi
    %\usepackage{chapterbib}
\fi

%
% Nomenclature; run makeindex with appropriate style file
%
\IfFileExists{./\jobname.nlo}{
    \if@dutch
        \execute{makeindex -q \jobname.nlo -s deltares_nomentbl_nl.ist -o \jobname.nls}
    \else
        \execute{makeindex -q \jobname.nlo -s deltares_nomentbl.ist -o \jobname.nls}
    \fi
}{}

\def\@@nomenclature[#1]#2#3#4{\endgroup\@esphack}
\def\@@@nomenclature[#1]#2#3#4{%
\def\@tempa{#2}\def\@tempb{#3}%
\protected@write\@nomenclaturefile{}%
{\string\nomenclatureentry{#1\nom@verb\@tempa @{\nom@verb\@tempa}&%
\begingroup\nom@verb\@tempb\endgroup &\begingroup#4\endgroup&%
\begingroup\protect\nompageref}{0}% nompageref==0, so all entries of the page number are the same, needed to skip double entries
}%
\endgroup
\@esphack}


\makenomenclature

\makeatother

