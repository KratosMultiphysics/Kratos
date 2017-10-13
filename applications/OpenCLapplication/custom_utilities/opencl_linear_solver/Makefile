EXE = opencl_linear_solver
INC = -I ../../../../kratos/includes/  -I ../../../../kratos/ -I /usr/local/cuda/include -I /opt/AMDAPP/include -I ..

CC = gcc
CCFLAGS = -O3 -DNDEBUG -fopenmp

CXX = g++
CXXFLAGS = -O3 -DNDEBUG -fopenmp

LD = g++
LIB = -L /opt/AMDAPP/lib/x86_64 -lOpenCL -lrt -lgomp

SRC = main.cpp mmio.c
OBJ = main.o mmio.o

MatrixPath = ../../../../../matrices

# Support for CGTYPE, MAXIT and MAXERR

ifdef CGTYPE
	DEFS += -DCGTYPE=$(CGTYPE)
endif

ifdef MAXIT
	DEFS += -DMAXIT=$(MAXIT)
endif

ifdef MAXERR
	DEFS += -DMAXERR=$(MAXERR)
endif

all: $(EXE)

.c.o:
	$(CC) -c $(CCFLAGS) $(INC) $(DEFS) $< -o $@

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(INC) $(DEFS) $< -o $@

$(EXE): $(OBJ)
	$(LD) -o $(EXE) $(OBJ) $(LIB)

run1: $(EXE)
	./$(EXE) A_0.mm B_0.mm

run2: $(EXE)
	./$(EXE) mat65k.mm rhs65025.mm

FEMCantilever: $(EXE)
	./$(EXE) $(MatrixPath)/cant.mtx $(MatrixPath)/cant-rhs.mtx

FEMSpheres: $(EXE)
	./$(EXE) $(MatrixPath)/consph.mtx $(MatrixPath)/consph-rhs.mtx

FEMAccelerator: $(EXE)
	./$(EXE) $(MatrixPath)/cop20k_A.mtx $(MatrixPath)/cop20k_A-rhs.mtx

Dense: $(EXE)
	./$(EXE) $(MatrixPath)/dense2.mtx $(MatrixPath)/dense2-rhs.mtx

Economics: $(EXE)
	./$(EXE) $(MatrixPath)/mac_econ_fwd500.mtx $(MatrixPath)/mac_econ_fwd500-rhs.mtx

Epidemiology: $(EXE)
	./$(EXE) $(MatrixPath)/mc2depi.mtx $(MatrixPath)/mc2depi-rhs.mtx

Protein: $(EXE)
	./$(EXE) $(MatrixPath)/pdb1HYS.mtx $(MatrixPath)/pdb1HYS-rhs.mtx

WindTunnel: $(EXE)
	./$(EXE) $(MatrixPath)/pwtk.mtx $(MatrixPath)/pwtk-rhs.mtx

QCD: $(EXE)
	./$(EXE) $(MatrixPath)/qcd5_4.mtx $(MatrixPath)/qcd5_4-rhs.mtx

LP: $(EXE)
	./$(EXE) $(MatrixPath)/rail4284.mtx $(MatrixPath)/rail4284-rhs.mtx

FEMHarbor: $(EXE)
	./$(EXE) $(MatrixPath)/rma10.mtx $(MatrixPath)/rma10-rhs.mtx

Circuit: $(EXE)
	./$(EXE) $(MatrixPath)/scircuit.mtx $(MatrixPath)/scircuit-rhs.mtx

FEMShip: $(EXE)
	./$(EXE) $(MatrixPath)/shipsec1.mtx $(MatrixPath)/shipsec1-rhs.mtx

Webbase: $(EXE)
	./$(EXE) $(MatrixPath)/webbase-1M.mtx $(MatrixPath)/webbase-1M-rhs.mtx
