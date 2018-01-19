# if you have an open-mp enabled compiler
CPP=g++  -fopenmp
# if you don't
#CPP=g++ -DNO_THREADS -g
REFCATPATH=$(shell pwd)/refcat/

LIBS=-lCCfits -L/home/rigel/dgruen/werc3/lib -I/home/rigel/dgruen/werc3/include -L/u/ki/dgruen/lib -I/u/ki/dgruen/include

all: betatree referencecats


betatree: src/betatree.cpp src/filter.cpp src/filter.h src/cfhtlib.h
	$(CPP) $(LIBS) -DREFCATPATH=\"$(REFCATPATH)\" -o betatree src/betatree.cpp src/filter.cpp

betatree_buzzard: src/betatree_buzzard.cpp src/filter.cpp src/filter.h src/cfhtlib.h
	$(CPP) $(LIBS) -DREFCATPATH=\"\" -o betatree_buzzard src/betatree_buzzard.cpp src/filter.cpp

referencecats:
	mkdir -p refcat
	if [ ! -s refcat/D.cat ]; then wget --no-check-certificate https://slac.stanford.edu/~dgruen/refcat/D.cat.gz; mv D.cat.gz refcat; gunzip refcat/D.cat.gz; fi
	if [ ! -s refcat/D.spectra.good.cat ]; then wget --no-check-certificate https://slac.stanford.edu/~dgruen/refcat/D.spectra.good.cat.gz; mv D.spectra.good.cat.gz refcat; gunzip refcat/D.spectra.good.cat.gz; fi
	if [ ! -s refcat/D.C2015.i.cat ]; then wget --no-check-certificate https://slac.stanford.edu/~dgruen/refcat/D.C2015.i.cat.gz; mv D.C2015.i.cat.gz refcat; gunzip refcat/D.C2015.i.cat.gz; fi

