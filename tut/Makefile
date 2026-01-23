txt = tut
nw = $(shell which noweb)

all: cross.sh sda.sh eda.sh neiSim.sh alfy2bed.awk senspec.awk
cross.sh: $(txt).org
	if [ "$(nw)" != "" ]; then\
		bash ../scripts/org2nw $(txt).org | notangle -Rcross.sh > cross.sh;\
	fi
sda.sh: $(txt).org
	if [ "$(nw)" != "" ]; then\
		bash ../scripts/org2nw $(txt).org | notangle -Rsda.sh > sda.sh;\
	fi
eda.sh: $(txt).org
	if [ "$(nw)" != "" ]; then\
		bash ../scripts/org2nw $(txt).org | notangle -Reda.sh > eda.sh;\
	fi
neiSim.sh: $(txt).org
	if [ "$(nw)" != "" ]; then\
		bash ../scripts/org2nw $(txt).org | notangle -RneiSim.sh > neiSim.sh;\
	fi
alfy2bed.awk: $(txt).org
	if [ "$(nw)" != "" ]; then\
		bash ../scripts/org2nw $(txt).org | notangle -Ralfy2bed.awk > alfy2bed.awk;\
	fi
senspec.awk: $(txt).org
	if [ "$(nw)" != "" ]; then\
		bash ../scripts/org2nw $(txt).org | notangle -Rsenspec.awk > senspec.awk;\
	fi
clean:
	rm -f cross.sh sda.sh eda.sh neiSim.sh alfy2bed.awk senspec.awk
