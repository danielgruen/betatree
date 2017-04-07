

all: executable referencecats


executable:


referencecats:
	mkdir -p refcat
	if [ ! -s refcat/D.cat ]; then wget --no-check-certificate https://slac.stanford.edu/~dgruen/refcat/D.cat.gz; mv D.cat.gz refcat; gunzip refcat/D.cat.gz; fi

