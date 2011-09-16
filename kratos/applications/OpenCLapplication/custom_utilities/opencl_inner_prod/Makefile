EXE = opencl_inner_prod
INC = -I ../../../../kratos/includes/  -I ../../../../kratos/ -I /usr/local/cuda/include -I ~/AMD-APP-SDK-v2.5-RC2-lnx64/include -I ..

CC = gcc
CCFLAGS = -O3 -DNDEBUG

CXX = g++
CXXFLAGS = -O3 -DNDEBUG

LD = g++
LIB = -L ~/AMD-APP-SDK-v2.5-RC2-lnx64/lib/x86_64 -lOpenCL -lrt

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
	LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/AMD-APP-SDK-v2.5-RC2-lnx64/lib/x86_64 ./$(EXE) B_0.mm B_0.mm

run2: $(EXE)
	LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/AMD-APP-SDK-v2.5-RC2-lnx64/lib/x86_64 ./$(EXE) rhs65025.mm rhs65025.mm

mc1: $(EXE)
	LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/AMD-APP-SDK-v2.5-RC2-lnx64/lib/x86_64 $(MC) ./$(EXE) B_0.mm B_0.mm

mc2: $(EXE)
	LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/AMD-APP-SDK-v2.5-RC2-lnx64/lib/x86_64 $(MC) ./$(EXE) rhs65025.mm rhs65025.mm
