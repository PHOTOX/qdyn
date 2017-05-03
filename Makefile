# Makefile for Qdyn program

# Compiler
FC = gfortran
FFLAGS = 
# Name of program
OUT = qdyn

# List of all files (NAME ALL F90)
OBJS = init.o qdyn.o

myprogram: $(OBJS) 
	   $(FC) -o $(OUT) ${FFLAGS} $(OBJS) 

clean :
	/bin/rm -f *.o *.mod $(OUT)


.PHONY: clean 

.SUFFIXES: .F90 

.F90.o:
	$(FC) ${FFLAGS} -c $<

