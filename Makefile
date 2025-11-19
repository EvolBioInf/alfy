all: alfy

alfy: bin/alfy
bin/alfy:
	if [ ! -d bin ]; then \
		mkdir bin; \
	fi
	make -C src
	cp src/alfy bin
data:
	curl https://owncloud.gwdg.de/index.php/s/ch7WkkXD5GLEjJ7/download -o alfyData.tgz
	tar -xvzf alfyData.tgz
	rm alfyData.tgz
.PHONY: test
test:
	make test -C src/
clean:
	make clean -C src
	rm bin/*
