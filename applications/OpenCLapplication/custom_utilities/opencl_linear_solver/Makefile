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

# Only on NVIDIA systems
MC = cuda-memcheck

all: $(EXE)

.c.o:
	$(CC) -c $(CCFLAGS) $(INC) $< -o $@

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(INC) $< -o $@

$(EXE): $(OBJ)
	$(LD) -o $(EXE) $(OBJ) $(LIB)

run1: $(EXE)
	LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):/opt/AMDAPP/lib/x86_64 ./$(EXE) A_0.mm B_0.mm

run2: $(EXE)
	LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):/opt/AMDAPP/lib/x86_64 ./$(EXE) mat65k.mm rhs65025.mm

mc1: $(EXE)
	LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):/opt/AMDAPP/lib/x86_64 $(MC) ./$(EXE) A_0.mm B_0.mm

mc2: $(EXE)
	LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):/opt/AMDAPP/lib/x86_64 $(MC) ./$(EXE) mat65k.mm rhs65025.mm
