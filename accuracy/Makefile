nw = $(shell which noweb)
all: accuracy.awk
accuracy.awk: accuracy.org
	if [ "$(nw)" != "" ]; then \
		bash ../scripts/org2nw accuracy.org | notangle -Raccuracy.awk > accuracy.awk; \
	fi
test: accuracy.awk
	bash test.sh
clean:
	rm -f accuracy.awk

