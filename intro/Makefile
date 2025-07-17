all: intro.sh

intro.sh: intro.org
	bash ../scripts/org2nw intro.org | notangle -Rintro.sh > intro.sh

clean:
	rm intro.sh
