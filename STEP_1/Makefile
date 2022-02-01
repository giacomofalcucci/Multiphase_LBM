#----------------------------------------------
DATE = $(shell date +%d%m%y)
#----------------------------------------------
FINE = "********** Fine Compilazione *********"
#----------------------------------------------
#
# generic option
FOPT = -O3
COPT = -O3
OPT = -O3
CUDA = cc70
OPENMP = -qopenmp
#
# generic compiler
FC = gfortran
CC = gcc
# 
#
OBJ = muphase.o \
      input.o \
      inithydro1.o \
      initpop.o \
      move.o \
      hydrovar1.o \
      hydrovar2.o \
      equili.o \
      collis.o \
      media.o \
      force.o \
      pbc.o \
      pbcdens.o \
      mbc.o \
      obstbc.o \
      output.o \
      savepop.o \
      resume.o \
      calcsurf.o \
      fourier.o \
      laplace.o

#----------------------------------------------

all: $(OBJ) 
	$(FC) $(FOPT) $(OBJ) $(FIX) -o muphase.x 

#----------------------------------------------

debug: FOPT = -O0 -g -C
debug: $(OBJ) 
	$(FC) $(FOPT) $(OBJ) $(FIX) -o muphase.x 

#----------------------------------------------
info: 
	@echo "Objects          =  "$(OBJ);
	@echo "Compiler         =  "$(FC);
	@echo "Compiler flags   =  "$(FOPT);
	@echo "Try      flags   =  "$(TRY);

#----------------------------------------------
clean:
	rm -rf *.o; 
	rm -rf muphase.x;

#----------------------------------------------

%.o %.mod: $(INC) %.f90
	$(FC) $(FOPT) $(FIX) -c $<

%.o: $(INC) %.F90
	$(FC) $(FOPT) $(FIX) -c $<

%.o: $(INC) %.f
	$(FC) $(FOPT) $(FIX) -c $<

%.o: $(INC) %.c
	$(CC) $(COPT)        -c $<


