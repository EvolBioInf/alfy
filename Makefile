progs = alfy ms2nn
shellScripts = cross neiSim
awkScripts = accuracy
all: data
	test -d bin || mkdir bin
	for shellScript in $(shellScripts); do \
		make -C $$shellScript; \
		cp $$shellScript/$$shellScript.sh scripts; \
		cp $$shellScript/$$shellScript.sh playground;\
	done
	for awkScript in $(awkScripts); do \
		make -C $$awkScript; \
		cp $$awkScript/$$awkScript.awk scripts; \
		cp $$awkScript/$$awkScript.awk playground;\
	done
	cp scripts/quantifyGenotypes.awk playground
	make -C src
	cp src/alfy bin
	make -C ms2nn
	cp ms2nn/ms2nn bin
	make -C tut
	cp tut/testAlfy.sh playground
data: playground
	curl https://owncloud.gwdg.de/index.php/s/ch7WkkXD5GLEjJ7/download -o alfyData.tgz
	tar -xvzf alfyData.tgz
	rm alfyData.tgz
	cp data/A+DQ083238.fasta playground
	cp data/hiv42.fasta playground
playground:
	test -d playground || mkdir playground
.PHONY: test
test: data
	make test -C src/
clean:
	make clean -C src
	make clean -C ms2nn
	for shellScript in $(shellScripts); do \
		make clean -C $$shellScript; \
		rm -f scripts/$$shellScript.sh; \
	done
	for awkScript in $(awkScripts); do \
		make clean -C $$awkScript; \
		rm -f scripts/$$awkScript.awk; \
	done
	rm -fr bin/* playground
