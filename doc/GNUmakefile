target = UsingSOMAR

src_files = $(shell find -maxdepth 16 -name \*.tex -type f | tr '\n' ' ')
bib_files = $(shell find -maxdepth 16 -name \*.bib -type f | tr '\n' ' ')
misc_files = ./somar.sty

tex_compiler = pdflatex -synctex=1 -interaction=nonstopmode
bibtex_compiler = bibtex

.phony: clean deps

all: $(target).pdf

$(target).pdf:$(src_files) $(bib_files) $(misc_files)
	$(tex_compiler) $(target).tex
	$(bibtex_compiler) $(target).aux
	$(tex_compiler) $(target).tex
	$(tex_compiler) $(target).tex

clean:
	rm -f *.aux  *.dvi  *.log  *.out  *.pdf  *.synctex.gz  *.toc *.bbl *.blg

deps:
	@echo $(src_files) $(bib_files) $(misc_files)
