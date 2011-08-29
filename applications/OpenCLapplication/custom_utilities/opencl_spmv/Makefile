EXE = opencl_spmv
INC = -I ../../../../kratos/includes/  -I ../../../../kratos/ -I /home/mossaiby/NVIDIA_GPU_Computing_SDK/OpenCL/common/inc/ -I ..

CC = gcc
CCFLAGS = -O3 -DNDEBUG

CXX = g++
CXXFLAGS = -O3 -DNDEBUG

LD = g++
LIB = -L /usr/local/cuda/lib64/ -lOpenCL -lrt

SRC = main.cpp mmio.c
OBJ = main.o mmio.o

MC = cuda-memcheck

all: $(EXE)

.c.o:
	$(CC) -c $(CCFLAGS) $(INC) $< -o $@

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(INC) $< -o $@

$(EXE): $(OBJ)
	$(LD) -o $(EXE) $(OBJ) $(LIB)

run1: $(EXE)
	LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/lib64 ./$(EXE) A_0.mm B_0.mm

run2: $(EXE)
	LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/lib64 ./$(EXE) mat65k.mm rhs65025.mm

mc1: $(EXE)
	LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/lib64 $(MC) ./$(EXE) A_0.mm B_0.mm

mc2: $(EXE)
	LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/lib64 $(MC) ./$(EXE) mat65k.mm rhs65025.mm
