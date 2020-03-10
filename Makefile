
MAKE = make
UTILS := spacegroup pp

include make.inc

all:
	cd src; $(MAKE)
	for UTIL in $(UTILS); do cd utilities/$$UTIL; $(MAKE); cd ../..; done

clean:
	cd src; $(MAKE) clean

test:
	cd tests; ./tests.sh

install:
	cp src/elk-cpu bin/elk
	for UTIL in $(UTILS); do cd utilities/$$UTIL; $(MAKE) install; cd ../..; done

docs:
	[[ -d docs ]] || mkdir docs
	cd src; $(MAKE) doc; cp elk.pdf ../docs/
	cd ../utilities/spacegroup; $(MAKE) doc; cp spacegroup.pdf ../../docs/
	cp src/addons/CRPA-Calculation.pdf docs/

lsvars:
	echo "export LIBS=\"$(LIBS)\"" > libs.sh; chmod +x libs.sh
	echo "export F90_OPTS=\"$(F90_OPTS)\"" > f90opts.sh; chmod +x f90opts.sh
