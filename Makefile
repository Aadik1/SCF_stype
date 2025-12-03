#DEBUG   = -nbs -C -d2
#FFLAGC  = -g -m64
FFLAGC  = -O2 -fopenmp
opt     = 
#profile = -stack:fa00000
obj     = run.out
ildir   =  
LIBS    = -L/usr/local/opt/openblas/lib -lopenblas 
input   = 

CFLAGC  = 

#FCOMPL  = af90  $(DEBUG) 
FCOMPL  = gfortran
CCOMPL  = 
# list of other directories for source files
.PREFIXES: .

.SUFFIXES:
.SUFFIXES: .f90 .c .s .o .fil

.f90.o:
	$(FCOMPL) -c $(ildir) $(FFLAGC) $(opt) $(profile) $<

OBJECTS = lapack_mods.o algebra.o Hamiltonian.o mod_Gfunctions.o \
          stypesub.o StypeJunc.o 

INLINE  = 

APPLIC: $(INLINE) $(OBJECTS)  $(OBJECTS1)  $(OBJECTS2)  $(OBJECTS3)
	$(FCOMPL) $(FFLAGC) $(profile) -o $(obj)  $(OBJECTS)  $(LIBS)

test:
	@echo START TEST ON $(input) , opt = $(opt)
	@echo start test on $(input) , opt = $(opt) >> TIME.LOG
	@date >> TIME.LOG
	@( time $(obj) < $(input) > $(input).out ) 2>> TIME.LOG
	@echo - - - - - - - - - - - >> TIME.LOG

clean:
	@rm -f $(INLINE) $(OBJECTS) $(OBJECTS1) sf_* *.o *.mod *.out *~ iteration* fort* err* Volt_Current* Print.dat SpinJunction.* populations.*
	clear

# include file dependencies

lapack_mods.o : lapack_mods.f90
algebra.o     : algebra.f90
Hamiltonian.o : Hamiltonian.f90
mod_Gfunctions.o : mod_Gfunctions.f90 Hamiltonian.f90
StypeJunc.o : StypeJunc.f90 mod_Gfunctions.f90 Hamiltonian.f90
stypesub.o : stypesub.f90 mod_Gfunctions.f90 Hamiltonian.f90

