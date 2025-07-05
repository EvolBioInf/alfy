date = $(shell git log | grep -m 1 Date | sed -r 's/Date: +[A-Z][a-z]+ ([A-Z][a-z]+) ([0-9]+) [^ ]+ ([0-9]+) .+/\2_\1_\3/')
version = $(shell git describe)

all: alfyDoc.pdf

alfyDoc.pdf: alfyDoc.tex header.tex homol.tex
	echo $(date) | tr '_' ' ' > date.txt
	echo $(version) | tr '-' ' ' | awk '{printf "%s", $$1; if ($$2) printf "-%s", $$2; printf "\n"}' > version.txt
	latex alfyDoc
	bibtex alfyDoc
	latex alfyDoc
	latex alfyDoc
	dvips -o -q alfyDoc
	ps2pdf -dALLOWPSTRANSPARENCY alfyDoc.ps
publish: alfyDoc.pdf
	if [ -d ~/WinHome/ownCloud ]; then \
		cp alfyDoc.pdf ~/WinHome/ownCloud/docs/; \
	fi
clean:
	rm -f alfyDoc.aux alfyDoc.bbl alfyDoc.bib alfyDoc.blg alfyDoc.d alfyDoc.dvi \
	alfyDoc.log alfyDoc.pdf alfyDoc.ps
