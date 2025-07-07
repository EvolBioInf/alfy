EXECFILE = alfy
EXECFILE64 = alfy64

DIRECTORY = alfy
################################################

VERSION = 1.5
.PHONY : all
all : $(EXECFILE) $(EXECFILE64)

# 32 bit version
$(EXECFILE):
	cd src/alfy; make; cp alfy ../../
# 64 bit version
$(EXECFILE64):
	cd src/alfy; make; cp alfy64 ../../

###############################################

data:
	curl https://owncloud.gwdg.de/index.php/s/ch7WkkXD5GLEjJ7/download -o alfyData.tgz
	tar -xvzf alfyData.tgz
	rm alfyData.tgz

test:
	make test -C src/alfy

#
# Other Standard make rules
#
clean:
	rm $(EXECFILE) $(EXECFILE64)
	cd externSrc/deepShallow/src/; make clean;
	cd externSrc/deepShallow64/src; make clean;
	cd src/common/interval; make clean;
	cd src/common/sequence/; make clean;
	cd src/common/util/; make clean;
	cd src/alfy; make clean;
