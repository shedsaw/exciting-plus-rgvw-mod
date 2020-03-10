include ../../make.inc

.SUFFIXES: .o .f90
.f90.o:
	$(F90SERIAL) $(F90_OPTS_SERIAL) -c -o $(*F).o $<

all: pdos bndchr

pdos: addorb.o pdos.o
	$(F90SERIAL) -o $@ addorb.o pdos.o

bndchr: addorb.o bndchr.o
	$(F90SERIAL) -o $@ addorb.o bndchr.o

install: pdos bndchr
	[[ -d ../../bin ]] || mkdir -p ../../bin
	cp ./pdos ../../bin/pdos-${COMPILER}
	cp ./bndchr ../../bin/bndchr-${COMPILER}

clean:
	rm pdos bndchr *.o
