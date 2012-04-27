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

#define EMUSYNC

#else

#define EMUSYNC

#endif

//defines for preconditioner
#define HWS 16
#define HWSB 4
#define HWSM 15
#define BS 256
#define WS 32
#define WSM 31
#define WSB 5

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

/** ADDED FUNCTIONS **/

void CUDA_CHECK(cudaError_t ErrorCode)
{
    if (ErrorCode != cudaSuccess)
        printf("*** CUDA error: %s ***\n", cudaGetErrorString(ErrorCode));
}

void sortElems_CPU(size_t* array, double* array2, size_t numElems)
{
    size_t i, temp, flag = 1, numLength = numElems;
    size_t d = numLength;
    double tempData;
    while( flag || (d > 1))      // boolean flag (true when not equal to 0)
    {
        flag = 0;           // reset flag to 0 to check for future swaps
        d = (d+1) / 2;
        for (i = 0; i < (numLength - d); i++)
        {
            if (array[i + d] < array[i])
            {
                temp = array[i + d];      // swap positions i+d and i
                array[i + d] = array[i];
                array[i] = temp;
                tempData = array2[i+d];
                array2[i+d] = array2[i];
                array2[i] = tempData;
                flag = 1;                  // tells swap has occurred
            }
        }
    }
}
void sortElems_CPU(size_t* array, size_t numElems)
{
    size_t i, temp, flag = 1, numLength = numElems;
    size_t d = numLength;
    while( flag || (d > 1))      // boolean flag (true when not equal to 0)
    {
        flag = 0;           // reset flag to 0 to check for future swaps
        d = (d+1) / 2;
        for (i = 0; i < (numLength - d); i++)
        {
            if (array[i + d] < array[i])
            {
                temp = array[i + d];      // swap positions i+d and i
                array[i + d] = array[i];
                array[i] = temp;
                flag = 1;                  // tells swap has occurred
            }
        }
    }
}
/** implementation of quicksort algorithm **/
void sortMatrix(size_t *CPU_RowIndices, size_t *CPU_Columns, double *CPU_Values, size_t Size1, size_t NNZ, bool sortingValues)
{
    switch(sortingValues)
    {
    case false:
        for(size_t i = 0; i < Size1; i++)
        {
            size_t firstElem = CPU_RowIndices[i];
            size_t numElems = CPU_RowIndices[i+1] - firstElem;
            sortElems_CPU(&CPU_Columns[firstElem], numElems);
        }
        break;
    case true:
        for(size_t i = 0; i < Size1; i++)
        {
            size_t firstElem = CPU_RowIndices[i];
            size_t numElems = CPU_RowIndices[i+1] - firstElem;
            sortElems_CPU(&CPU_Columns[firstElem], &CPU_Values[firstElem], numElems);
        }
        break;
    }

}

void printVector(const GPUVector& vec)
{
    //cout << "Num elems: " << vec.Size << endl;
    printf("Num elems: %u\n", vec.Size);
    for(size_t i = 0; i < vec.Size; i++)
    {
        //cout << vec.CPU_Values[i] << " ";
        printf("%lg ", vec.CPU_Values[i]);
    }
    printf("\n");
}

void printMatrix(const GPUCSRMatrix& mat)
{
    size_t index = 0;
    //cout << "Num rows: " << mat.Size1 << ", Num cols: " << mat.Size2 << ", numNNZ: " << mat.NNZ << endl;
    printf("Num rows: %u, Num cols: %u, numNNZ %u\n", mat.Size1, mat.Size2, mat.NNZ);
    for(size_t i = 0; i < mat.Size1; i++)
    {
        size_t currentCol = 0;
        size_t numPtr = mat.CPU_RowIndices[i+1] - mat.CPU_RowIndices[i];
        while(currentCol < mat.Size2 && numPtr > 0)
        {
            if(mat.CPU_Columns[index] == currentCol)
            {
//                cout << mat.CPU_Values[index] << " ";
                printf("%lg ", mat.CPU_Values[index]);
                index++;
                numPtr--;
            }
            else
                printf("0 ");
            currentCol++;
        }
        while(currentCol < mat.Size2)
        {
            //cout << "0" << " ";
            printf("0 ");
            currentCol++;
        }
        printf("\n");
    }
    printf("\n");
}

void printMatrix_csr(const GPUCSRMatrix& mat)
{
    for(size_t i = 0; i < mat.Size1; i++)
    {
        for(size_t j = mat.CPU_RowIndices[i]; j < mat.CPU_RowIndices[i+1]; j++)
        {
            //cout << "(" << i << ", " << mat.CPU_Columns[j] << ") -> " << mat.CPU_Values[j] << endl;
            printf("(%u, %u) -> %lg\n", i, mat.CPU_Columns[j], mat.CPU_Values[j]);
        }
    }
}

double computeSingleSize(GPUCSRMatrix& m)
{
    return (((double)m.NNZ * (sizeof(size_t) + sizeof(double)) + m.Size1 * sizeof(size_t))/1024)/1024;
}

double computeSingleSize(GPUVector& v)
{
    return (((double)v.Size * (sizeof(double))) / 1024) / 1024;
}

}
}

//
// To use texture memory in GPUGPUCSRMatrixVectorMultiply_CSR_Kernel; using two int2 to read a double
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



