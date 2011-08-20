#
# Makefile per convertire automaticamente il file LaTeX nei formati 
# usati per la pubblicazione.
#
# $Id: Makefile,v 1.1.1.1 2002/01/10 10:01:51 ri Exp $
#

NAME = fdgc
TEX = index.tex intro.tex altgc.tex feynman.tex \
	hermitian-model.tex appoperads.tex
STY = miscmath.sty rg.sty paralist.sty prettyref.sty
BIB = index.bbl
FMT = 
DVI = index.dvi

clean:
	rm -f *~ *.aux *.dvi *.toc *.lof *.lot *.lol *.lox *.log *.out \
		*.bbl *.blg *.glo *.gls

tar: index.bbl
	tar -cf $(NAME).tar $(TEX) $(STY) $(BIB) $(FMT) $(DVI)

tar.gz:
	gzip $(NAME).tar

%.aux: %.tex
	latex $<

%.bbl:
	bibtex $<

%.gls: %.glo
	makeindex $< -s nomencl.ist -o $<

%.bbl: %.aux
	bibtex $<

%.dvi: %.aux
	latex $< || latex $<

%.ps: %.dvi alleps
	dvips -o $@ $<

%.eps: %.xwd
	convert $< $@
	-cvs add $@

%.png: %.xwd
	convert $< $@
	-cvs add $@

alleps:
	for x in `find -name '*.xwd'`; do\
		make -C `dirname $$x` `basename $$x .xwd`.eps;\
	done

epshere:
	make `find . -name '*.xwd' -print`

allpng:
	for x in `find -name '*.xwd'`; do\
		make -C `dirname $$x` `basename $$x .xwd`.png;\
	done

pnghere:
	make `find . -name '*.xwd' -print`

new:
	test "$(dir)"
	mkdir $(dir)
	cp Makefile $(dir)
	cp body.tex $(dir).tex
	cvs new $(dir)
	cvs new $(dir).tex

newdir:
	test "$(dir)"
	mkdir $(dir)
	mkdir $(dir)/auto
	touch $(dir)/auto/body.el
	cp Makefile $(dir)
	sed -e 's|index|../index|' <body.tex >$(dir)/body.tex
	cvs new $(dir)
	cvs new $(dir)/auto
	cvs new $(dir)/auto/body.el
	cvs new $(dir)/body.tex

propagate:
	test "$(file)"
	for x in `find -type d \
		-not -name . -and -not -name CVS -and -not -name auto`; do\
		cp -iv $(file) $$x;\
	done

html:
	latex2html -latex -tex_defs -white -address "Riccardo Murri <r.murri@oltrelinux.com>" -split 0 -mkdir -dir HTML -no_navigation -html_version 2.0 index.tex

%.txt: HTML/%.html
	lynx -dump HTML/index.html > HTML/index.raw
	sed -e 's/LATEX/LaTeX/;s/TEX/TeX/' < HTML/index.raw > HTML/index.txt

images:
	latex2html -latex -tex_defs -white -address "Riccardo Murri <r.murri@oltrelinux.com>" -split 0 -mkdir -dir HTML -no_navigation -html_version 2.0 -images_only index.tex


