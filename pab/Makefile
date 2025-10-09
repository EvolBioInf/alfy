nw = $(shell which noweb)
all: pab.sh
pab.sh: pab.org
	if [ "$(nw)" != "" ]; then\
		bash ../scripts/org2nw pab.org | notangle -Rpab.sh > pab.sh;\
	fi
clean:
	rm -f pab.sh m.fa
