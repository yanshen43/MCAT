###############################################################################
# Makefile for a generic LaTeX document with figures and a bibliography.
#   The figures are in a directory Figures, which must have a Makefile.
#   Web of Science bibliography files must be in a directory WoS and
#   should be named *_ISIWoS.txt.
#   The PROJECT variable below can be changed to be more descriptive
#   of the LaTeX document.
#
# Copyright, 2010	Lenwood S. Heath
#
###############################################################################

LATEX_BATCHMODE  = latex
LATEX_SCROLLMODE = latex '\scrollmode \input '
BIBTEX = bibtex
DVIPS  = dvips -t letter -P pdf -GO
PS2PDF = ps2pdf

PROJECT = Paper1_V1
LATEX = $(PROJECT).tex
AUX = $(PROJECT).aux
BIB = $(PROJECT).bib
BBL = $(PROJECT).bbl
DVI = $(PROJECT).dvi
TOC = $(PROJECT).toc
PS = $(PROJECT).ps
PDF = $(PROJECT).pdf

PROTECT = keep_capital.txt
WOS = WoS/*_ISIWoS.txt

all: $(PDF)

$(PDF): $(PS)
	$(PS2PDF) $(PS)

$(PS): $(DVI)
	$(DVIPS) $(DVI)

$(BBL): figures $(BIB) $(AUX)
	$(BIBTEX) $(PROJECT)
	$(LATEX_BATCHMODE) $(PROJECT)

$(AUX): figures $(LATEX)
	$(LATEX_BATCHMODE) $(PROJECT)

$(DVI): figures $(LATEX) $(BBL)
	cd Figures;make
	$(LATEX_BATCHMODE) $(PROJECT)
	$(LATEX_SCROLLMODE) $(PROJECT)

figures:
	cd ./Figures;make

realclean: clean
	rm -f $(PS) $(PDF) $(TOC)
	cd ./Figures;make realclean
	rm -f *.pstex*

clean:
	rm -f *.log
	rm -f *.aux *.blg *.bbl *.dvi
	rm -f *.pdf
	rm -f *.tmp *.bak
	rm -f *.eps
