date = $(shell git log | grep -m 1 Date | sed -r 's/Date: +[A-Z][a-z]+ ([A-Z][a-z]+) ([0-9]+) [^ ]+ ([0-9]+) .+/\2_\1_\3/')
version = $(shell git describe)

src1 = header.tex ../intro/pic.tex
src2 = intro.tex alfy.tex sir.tex util.tex

all: alfyDoc.pdf

%.tex:	../*/%.org alfyDoc.tex
	bash ../scripts/org2nw $< | awk -f ../scripts/preWeave.awk | noweave -n -x | sed 's/_test/\\_test/g' > $@

alfyDoc.pdf: alfyDoc.tex $(src1) $(src2)
	echo $(date) | tr '_' ' ' > date.txt
	echo $(version) | tr '-' ' ' | awk '{printf "%s", $$1; if ($$2) printf "-%s", $$2; printf "\n"}' > version.txt
	latex alfyDoc
	bibtex alfyDoc
	latex alfyDoc
	latex alfyDoc
	dvips -o -q alfyDoc
	ps2pdf -dALLOWPSTRANSPARENCY alfyDoc.ps
publish: alfyDoc.pdf
	if [ -d ~/win/ownCloud ]; then \
		cp alfyDoc.pdf ~/WinHome/ownCloud/docs/alfyGoDoc.pdf; \
	fi
clean:
	rm -f alfyDoc.aux alfyDoc.bbl alfyDoc.bib alfyDoc.blg alfyDoc.d alfyDoc.dvi \
	alfyDoc.log alfyDoc.pdf alfyDoc.ps $(src2)
