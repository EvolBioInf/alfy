all: alfy ms2nn

alfy: bin/alfy
bin/alfy:
	if [ ! -d bin ]; then \
		mkdir bin; \
	fi
	make -C src
	cp src/alfy bin
ms2nn: bin/ms2nn
bin/ms2nn:
	if [ ! -d bin ]; then \
		mkdir bin; \
	fi
	make -C ms2nn
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
	rm -f bin/*
