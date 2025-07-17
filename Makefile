packs = util
progs = alfy sir

all:
	test -d bin || mkdir bin
	for pack in $(packs); do \
		make -C $$pack; \
	done
	for prog in $(progs); do \
		make -C $$prog; \
		cp $$prog/$$prog bin; \
	done

clean:
	rm -rf bin/*
	for pack in $(packs); do \
		make clean -C $$pack; \
	done
	for prog in $(progs) $(packs) doc; do \
		make clean -C $$prog; \
	done

