
## Options To Choose:
## -DNIBM : Without immersed boundary model
## -DIBM  : With immersed boundary model
## -DSI   : Semi-Implicit algorithm
## -DEXP  : Explicit algorithm
## -DTwoD : Two dimensional code



# fortran compiler
#FC = mpif90
FC   = ftn
#FC = mpiifort

OTH :=  -r8 -fpconstant  -132  -DEXP -DIBM  -DTwoD
#OTH := -ffixed-line-length-none -fdefault-real-8 -DSI -DIBM #-DNIBM  #-DIBM #-DSI  #-DEXP -DTwoD
OPT := -O3 
BIG :=  -mcmodel=medium
##DBG := -g -fbacktrace -Wall  #Gnu
DB := -traceback  -g   #Intel
PROF :=# -pg
OMP := #-openmp
FFLAGS:= -cpp $(OTH) $(BIG) $(BIG) $(DBG) $(OMP) $(PROF) $(OMP)  
# 2decomp&fft
include libs/2decomp_fft/src/Makefile.inc
INCLUDE = -I libs/2decomp_fft/include
#
LIB = -L libs/2decomp_fft/lib -l2decomp_fft -L libs/fft/lib -lfft

TARGET = inho

SRC = param.f90 common.f90 IBM.f90  initmpi.f90 interface.f90 debug.f90  initsolver.f90 bound.f90 initIBM.f90 chkdiv.f90 chkdt.f90 loadd.f90 CahnHilliardSI.f90  mom.f90 rk3.f90 fillps.f90 zredistribute.f90 solver.f90 correc.f90 post.f90  output.f90   vtk_write.f90 init.f90  Impl_PredVel.f90  main.f90

OBJ = $(SRC:.f90=.o)

all: $(TARGET)

$(TARGET): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(LIB)

.PHONY: clean  
clean:
	rm -f *.o *.mod $(TARGET) .depend

veryclean: 
	rm -f *.o *.mod $(TARGET) .depend libs/fft/fft.o /libs/fft/lib/libfft.a; cd libs/2decomp_fft/src/; make clean; cd ../../../

libraries:
	cd libs/fft; $(FC) $(FFLAGS) $(BIG) $(DBG) $(PROF) $(OMP) -c fft.f; ar qc libfft.a fft.o; mv libfft.a lib; cd ../../;
	cd libs/2decomp_fft/src/; make; cd ../../../

%.o: %.f90
	$(FC) $(INCLUDE) -c -o $@ $(FFLAGS) $(BIG) $(DBG) $(PROF) $(OMP) $<

.depend dep:
	./.makedepo $(SRC) > .depend

include .depend
