date = $(shell git log | grep -m 1 Date | sed -r 's/Date: +[A-Z][a-z]+ ([A-Z][a-z]+) ([0-9]+) [^ ]+ ([0-9]+) .+/\2_\1_\3/')
version = $(shell git describe)
ownCloud = ~/WinHome/ownCloud

all: alfyDoc.pdf

alfyDoc.pdf: alfy2bed.tex alfyDoc.tex cross.tex header.tex homol.tex ms2nn.tex neiSim.tex senSpec.tex tut.tex util.tex ../ms2nn/ms.tex ../ms2nn/tree.tex ../ms2nn/t1.tex ../ms2nn/t2.tex ../ms2nn/t3.tex
	echo $(date) | tr '_' ' ' > date.txt
	echo $(version) | tr '-' ' ' | awk '{printf "%s", $$1; if ($$2) printf "-%s", $$2; printf "\n"}' > version.txt
	latex alfyDoc
	bibtex alfyDoc
	latex alfyDoc
	latex alfyDoc
	dvipdf -dALLOWPSTRANSPARENCY alfyDoc
%.tex:	../*/%.org
	bash ../scripts/org2nw $< | awk -f ../scripts/preWeave.awk | noweave -n -x | sed 's/_test/\\_test/g' | sed 's/,\*/*/' > $@
publish: $(ownCloud)/docs/alfyDoc.pdf
$(ownCloud)/docs/alfyDoc.pdf: alfyDoc.pdf
	if [ -d $(ownCloud) ]; then \
		cp alfyDoc.pdf $(ownCloud)/docs/; \
	fi
clean:
	rm -f alfyDoc.aux alfyDoc.bbl alfyDoc.bib alfyDoc.blg alfyDoc.d alfyDoc.dvi \
	alfyDoc.log alfyDoc.pdf ms2nn.tex tut.tex util.tex
