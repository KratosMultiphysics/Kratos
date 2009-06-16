/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Farshid Mossaiby
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
mossaiby@yahoo.com
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

//
// Utility functions and macros for gpu_sparse.cu

// Divide X by Y, adding one if neccessary
#define DIVIDE_INTO(X, Y)	(((X) + (Y) - 1) / (Y))

// Returns the global index of the current thread
#ifdef USE_UMUL24

	#define GlobalIdx()		((__umul24(blockDim.x, blockIdx.x + __umul24(blockIdx.y, gridDim.x)) + threadIdx.x))

#else

	#define GlobalIdx()		((blockDim.x * (blockIdx.x + blockIdx.y * gridDim.x) + threadIdx.x))

#endif

// Some needed definitions
// TODO: Get MAX_BLOCKS using API

// Max possible number of blocks
#define MAX_BLOCKS			65536

// Size of a block (threads per block)
#define BLOCK_SIZE			256

// Size of a warp (threads per warp)
#define WARP_SIZE			32

// Warp size mask
#define WARP_SIZE_MASK		31

// Warp size log_2
#define WARP_SIZE_BITS		5

// Size of a half warp (threads per half warp)
#define HALF_WARP_SIZE		16

// Half warp size mask
#define HALF_WARP_SIZE_MASK	15

// Half warp size log_2
#define HALF_WARP_SIZE_BITS	4

// Macro to enforce intrawarp sychronization during emulation
#ifdef __DEVICE_EMULATION__

	#define EMUSYNC __syncthreads()

#else

	#define EMUSYNC

#endif

// Includes, system

#include "cuda.h"

namespace Kratos
{

namespace GPUSparse
{

//
// Build_Grid
// Given total number of threads, returns an appropriate dim3 grid object

dim3 Build_Grid(const unsigned int Threads, const unsigned int BlockSize)
{
	const unsigned int Blocks = DIVIDE_INTO(Threads, BlockSize);
	
	if (Blocks < MAX_BLOCKS)
	{
		// Fits in a 1D grid
		return dim3(Blocks);
	}
	else
	{
		// 2D grid is required
		const unsigned int Side = (unsigned int) ceil(sqrt(static_cast <double> (Blocks)));
		
		return dim3(Side, Side);
	}
}

//
// CUDA_Success
// A macro to check the return value of CUDA calls

bool inline CUDA_Success(cudaError_t ErrorCode)
{
	//if (ErrorCode != cudaSuccess)
		//printf("*** CUDA error: %s ***\n", cudaGetErrorString(ErrorCode));
		
	return ErrorCode == cudaSuccess;
}

bool inline CUBLAS_Success(cublasStatus ErrorCode)
{
	return ErrorCode == CUBLAS_STATUS_SUCCESS;
}

}

}

//
// To use texture memory in GPU_MatrixVectorMultiply_CSR_Kernel; using two int2 to read a double
texture <int2, 1>  Texture_double;

bool Bind_X(const double *X)
{
	return Kratos::GPUSparse::CUDA_Success(cudaBindTexture(NULL, Texture_double, X));
}

bool Unbind_X()
{
	return Kratos::GPUSparse::CUDA_Success(cudaUnbindTexture(Texture_double));
}

__inline__ __device__ double Fetch_X(const double *X, const size_t &i)
{
	int2 V = tex1Dfetch(Texture_double, i);
	return __hiloint2double(V.y, V.x);
}
