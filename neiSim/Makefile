nw = $(shell which noweb)
all: neiSim.sh
neiSim.sh: neiSim.org
	if [ "$(nw)" != "" ]; then \
		bash ../scripts/org2nw neiSim.org | notangle -RneiSim.sh > neiSim.sh; \
	fi
test: test.sh neiSim.sh
	bash test.sh
clean:
	rm -f neiSim.sh
