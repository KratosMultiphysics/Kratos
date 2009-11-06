# nvcc compiler to use
NVCC = nvcc

# nvcc flags
NVCC_FLAGS = -Xcompiler -fPIC -arch compute_13 -code compute_13

# Behavior of gpu_sparse library
#   -DUSE_UMUL24: Use __umul24 in GlobalIdx() instead of simple multiplication
#   -DUSE_TEXTURE_CACHING: Use texture caching in matrix vector multiplication kernels
DEF = -DUSE_UMUL24 -DUSE_TEXTURE_CACHING

# Main source file
SRC = gpu_sparse.cu

# Dependencies
DEP = gpu_sparse.h gpu_sparse_kernels.h gpu_sparse_utils.h spai_preconditioner.h

# Libraries to link with
LNK = -lcublas

# Object file name
OBJ = gpu_sparse.o

# Library file name
LIB = libgpu_sparse.a

# Defaults: no emulation, no debug data
emu = 0
dbg = 0

# Device emulation mode
ifeq ($(emu), 1)
	NVCC_FLAGS += -deviceemu
endif

# Debug mode
ifeq ($(dbg), 1)
	NVCC_FLAGS += -g
else
	NVCC_FLAGS += -O3
	DEF += -DNDEBUG
endif

$(LIB): $(OBJ)
	ar cr $(LIB) $(OBJ)

$(OBJ): $(SRC) $(DEP)
	$(NVCC) $(SRC) -c -o $(OBJ) $(NVCC_FLAGS) $(DEF) $(LNK)

clean:
	rm *.o *.a *.linkinfo
	
