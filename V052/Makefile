# This work is licensed under the Creative Commons Attribution-NonCommercial
# 3.0 Unported License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/3.0/.

.PHONY: all clean

all : ersatzschaltung.eps protokoll.pdf

clean:
	rm -f protokoll.aux protokoll.log protokoll.pdf	\
	protokoll.out protokoll.run.xml protokoll.toc	\
	protokoll.bcf ersatzschaltung.eps

%.pdf : %.tex
	pdflatex $<

%.eps : %.svgz
	inkscape --export-eps=$@ $<
