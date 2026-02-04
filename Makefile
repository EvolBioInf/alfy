progs = alfy ms2nn
shellScripts = cross neiSim
awkScripts = alfy2bed senSpec
all:
	test -d bin || mkdir bin
	for shellScript in $(shellScripts); do \
		make -C $$shellScript; \
		cp $$shellScript/$$shellScript.sh scripts; \
	done
	for awkScript in $(awkScripts); do \
		make -C $$awkScripts; \
		cp $$awkScript/$$awkScript.awk scripts; \
	done
	make -C src
	cp src/alfy bin
	make -C ma2nn
	cp ms2nn/ms2nn bin
data:
	curl https://owncloud.gwdg.de/index.php/s/ch7WkkXD5GLEjJ7/download -o alfyData.tgz
	tar -xvzf alfyData.tgz
	rm alfyData.tgz
.PHONY: test
test: data
	make test -C src/
clean:
	make clean -C src
	make clean -C ms2nn
	for shellScript in $(shellScripts); do \
		make clean -C $$shellScript; \
		rm scripts/$$shellScript.sh; \
	done
	for awkScript in $(awkScripts); do \
		make clean -C $$awkScripts; \
		rm scripts/$$awkShript.awk; \
	done
	rm -f bin/*
