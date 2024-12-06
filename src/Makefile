# Makefile for Qdyn program

# Compiler
FC = gfortran
FFLAGS =
# FFTW library needed
LIBS = -L/Users/janosj/Documents/Programs/fftw/fftw-3.3.10/lib -lfftw3 -lm -llapack
LDLIBS = -lfftw3 # This is necessary only on clusters, not on MacOS

# for Github tests, I use different fftw library path
GitHub_LIBS = -L/usr/local/lib -lfftw3 -lm -llapack

# Name of program
OUT = qdyn

#======================================================================

# List of all files (NAME ALL F90)
OBJS = vars.o fparser.o fftw.o utils.o init.o exactfactor.o propag.o qdyn.o

myprogram: $(OBJS) 
	   $(FC) -o $(OUT) ${LIBS} ${FFLAGS} $(OBJS) ${LDLIBS}

github-build: $(OBJS)
	   $(FC) -o $(OUT) ${GitHub_LIBS} ${FFLAGS} $(OBJS) ${LDLIBS}

clean:
	/bin/rm -f *.o *.mod $(OUT)

# run the full test suite
test:
	cd ../tests; /bin/bash run_test_suite.sh

# run the quick version of the test suite
qtest:
	cd ../tests; /bin/bash run_test_suite.sh 0

.PHONY: clean 

.SUFFIXES: .F90 .f90 .f03 

.F90.o:
	$(FC) ${FFLAGS} -c $<

.f90.o:
	$(FC) ${FFLAGS} -c $<

.f03.o:
	$(FC) ${FFLAGS} -c $<
