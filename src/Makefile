# Makefile for Qdyn program

# Compiler
FC = gfortran
FFLAGS =
#FFTW library needed
# LIBS = -L/home/janos/Programs/fftw/fftw-3.3.6/lib -lfftw3 -lm 
LIBS = -L/home/janos/Programs/fftw/fftw-3.3.10/lib -lfftw3 -lm 
LDLIBS = -lfftw3
#FFLAGS = -I/home/suchanj/programs/fftw-3.3.6/include/

# Name of program
OUT = qdyn

#======================================================================

# List of all files (NAME ALL F90)
OBJS = vars.o fparser.o fftw.o utils.o init.o propag.o qdyn.o

myprogram: $(OBJS) 
	   $(FC) -o $(OUT) ${LIBS} ${FFLAGS} $(OBJS) ${LDLIBS}

clean :
	/bin/rm -f *.o *.mod $(OUT)


.PHONY: clean 

.SUFFIXES: .F90 .f90 .f03 

.F90.o:
	$(FC) ${FFLAGS} -c $<

.f90.o:
	$(FC) ${FFLAGS} -c $<

.f03.o:
	$(FC) ${FFLAGS} -c $<
