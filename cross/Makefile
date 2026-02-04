nw = $(shell which noweb)
all: cross.sh cli.sh
cross.sh: cross.org
	if [ "$(nw)" != "" ]; then \
		bash ../scripts/org2nw cross.org | notangle -Rcross.sh > cross.sh; \
	fi
cli.sh: cross.org
	if [ "$(nw)" != "" ]; then \
		bash ../scripts/org2nw cross.org | notangle -Rcli > cli.sh; \
	fi
test: cross.sh
	bash test.sh
clean:
	rm -f cli.sh cross.sh
