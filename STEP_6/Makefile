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
OBJ = storage.o \
      muphase.o \
      input.o \
      alloca.o \
      inithydro1.o \
      initpop.o \
      move.o \
      hydrovarPRE.o \
      collis.o \
      pbc.o \
      pbcdens.o \
      mbc.o \
      obstbc.o \
      output.o \
      savepop.o \
      resume.o \
      calcsurf.o \
      laplace.o

#----------------------------------------------

all: $(OBJ) 
	$(FC) $(FOPT) $(FIX) $(OBJ) -o muphase.x 

#----------------------------------------------

debug: FOPT := -O0 -g -C -DDEBUG
debug: $(OBJ) 
	$(FC) $(FOPT) $(FIX) $(OBJ) -o muphase.x 

#----------------------------------------------
openmp: FOPT += $(OPENP)
openmp: $(OBJ)
	$(FC) $(FOPT) $(OPENMP) $(FIX) $(OBJ) -o muphase.x

#----------------------------------------------
info: 
	@echo "Objects          =  "$(OBJ);
	@echo "Compiler         =  "$(FC);
	@echo "Compiler flags   =  "$(FOPT);
	@echo "OpenMI   flag    =  "$(FOPT);
	@echo "Try      flag    =  "$(TRY);

#----------------------------------------------
clean:
	rm -rf *.o; 
	rm -rf *.mod; 
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


