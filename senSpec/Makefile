nw = $(shell which noweb)
all: senSpec.awk
senSpec.awk: senSpec.org
	if [ "$(nw)" != "" ]; then \
		bash ../scripts/org2nw senSpec.org | notangle -RsenSpec.awk > senSpec.awk; \
	fi
test: senSpec.awk
	bash test.sh
clean:
	rm -r senSpec.awk

