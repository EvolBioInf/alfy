txt = tut
nw = $(shell which noweb)

all: sda.sh eda.sh testAlfy.sh
sda.sh: $(txt).org
	if [ "$(nw)" != "" ]; then\
		bash ../scripts/org2nw $(txt).org | notangle -Rsda.sh > sda.sh;\
	fi
eda.sh: $(txt).org
	if [ "$(nw)" != "" ]; then\
		bash ../scripts/org2nw $(txt).org | notangle -Reda.sh > eda.sh;\
	fi
testAlfy.sh: $(txt).org
	if [ "$(nw)" != "" ]; then\
		bash ../scripts/org2nw $(txt).org | notangle -RtestAlfy.sh > testAlfy.sh;\
	fi
clean:
	rm -f sda.sh eda.sh testAlfy.sh
