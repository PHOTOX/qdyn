# Makefile for Qdyn program

# Compiler
FC = gfortran
FFLAGS =
#FFTW library needed
# LIBS = -L/home/janos/Programs/fftw/fftw-3.3.10/lib -lfftw3 -lm 
LIBS = -L/Users/janosj/Documents/Programs/fftw/fftw-3.3.10/lib -lfftw3 -lm
LDLIBS = -lfftw3 # This is necessary only on clusters, not on MacOS
#FFLAGS = -I/home/suchanj/programs/fftw-3.3.6/include/ # This was uncommented in the original make file by JS

# Name of program
OUT = qdyn

#======================================================================

# List of all files (NAME ALL F90)
OBJS = vars.o fparser.o fftw.o utils.o init.o propag.o qdyn.o

myprogram: $(OBJS) 
	   $(FC) -o $(OUT) ${LIBS} ${FFLAGS} $(OBJS) ${LDLIBS}

clean:
	/bin/rm -f *.o *.mod $(OUT)

test:
	cd ../tests; /bin/bash run_test_suite.sh


.PHONY: clean 

.SUFFIXES: .F90 .f90 .f03 

.F90.o:
	$(FC) ${FFLAGS} -c $<

.f90.o:
	$(FC) ${FFLAGS} -c $<

.f03.o:
	$(FC) ${FFLAGS} -c $<
