#FFLAGS  = -g -Og -fbacktrace
#FFLAGS = -m64
#FFLAGS = -O3 -g -I../src
FFLAGS = -O1 -g -fbounds-check
F77     = gfortran

%.o: %.F90 
	$(F77) $(FFLAGS) -c $<
VPATH = ../../src/common:../../src/ellipse

include dependencies.mk

OBJ= params.o types.o globalvars.o mod_stokes.o mod_SDLP.o mod_EHpqr.o mod_compcoeffdens.o mod_runge.o mod_products.o mod_geom.o mod_init.o mod_target.o runtime.o old.o

OBJtest= params.o types.o globalvars.o mod_EHpqr.o mod_geom.o mod_stokes.o mod_products.o mod_compcoeffdens.o mod_SDLP.o

#problems = shadow.o evolve.o fixedmax.o fixeddble.o

driver.out : $(OBJ) driver.o 
	$(F77) -o driver.out driver.o $(OBJ) 

testcorr.out : $(OBJtest) testcorr.o 
	$(F77) -o testcorr.out testcorr.o $(OBJtest) 

clean : 
	rm -f *.o
