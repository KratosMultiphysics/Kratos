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


#if !defined(KRATOS_GPU_SPARSE_H_INCLUDED)
#define  KRATOS_GPU_SPARSE_H_INCLUDED


#include <stddef.h>

#define GPU_CHECK(call)		if (!(call)) printf(#call " in file %s line %d failed.\n", __FILE__, __LINE__);

namespace Kratos
{

namespace GPUSparse
{

//
// CopyDirection
// Copy direction between CPU and GPU

enum CopyDirection
{
	CPU_GPU,
	GPU_CPU
};

//
// GPUVector
// GPUVector class


class GPUVector
{
public:
	size_t Size;
	double *CPU_Values;
	double *GPU_Values;
	bool Allocated;

	GPUVector(size_t _Size, double *_CPU_Values);
	GPUVector(size_t _Size);
	~GPUVector();

	bool GPU_Allocate();
	bool GPU_Free();

	bool Copy(CopyDirection Direction);
	bool CopyFromGPU(GPUVector &V);
};

//
// GPUCSRMatrix
// Sparse matrix class in CSR format

class GPUCSRMatrix
{
public:
	size_t NNZ, Size1, Size2;

	size_t *CPU_Columns, *CPU_RowIndices;
	double *CPU_Values;

	size_t *GPU_Columns, *GPU_RowIndices;
	double *GPU_Values;

	//variables for a dense representation
	bool haveDenseRepresentation;
	double* matAuxValues;
	size_t numValuesDenseRep;

	bool Allocated;

	GPUCSRMatrix(size_t _NNZ, size_t _Size1, size_t _Size2, size_t *_CPU_Columns, size_t *_CPU_RowIndices, double *_CPU_Values, bool _NZMultiple16 = true);
	GPUCSRMatrix(size_t _NNZ, size_t _Size1, size_t _Size2);
	~GPUCSRMatrix();

	bool GPU_Allocate();
	bool GPU_Free();

	//function for dense conversion
	bool GenerateDenseRepresentation(bool FortranRep = false);

	bool Copy(CopyDirection Direction, bool CopyValuesOnly);
	bool CopyFromGPU(GPUCSRMatrix &M, bool CopyStructure, bool CopyValues);
};

// Operations defined on GPUCSRMatrix and GPUVector

//
// CPUGPUCSRMatrixVectorMultiply
// Matrix-Vector multiply on CPU

bool CPUGPUCSRMatrixVectorMultiply(GPUCSRMatrix &A, GPUVector &X, GPUVector &Y);

//
// GPUGPUCSRMatrixVectorMultiply
// Matrix-Vector multiply on GPU

bool GPUGPUCSRMatrixVectorMultiply(GPUCSRMatrix &A, GPUVector &X, GPUVector &Y);

//
// CPUGPUCSRMatrixGetDiagonals
// Extract the diagonal elements of a matrix into a vector on CPU

bool CPUGPUCSRMatrixGetDiagonals(GPUCSRMatrix &A, GPUVector &X);

//
// GPUGPUCSRMatrixGetDiagonals
// Extract the diagonal elements of a matrix into a vector on GPU

bool GPUGPUCSRMatrixGetDiagonals(GPUCSRMatrix &A, GPUVector &X);

//
// CPUGPUCSRMatrixMatrixDiagonalMultiply
// Multiply a digonal matrix specified with a vector with a matrix on CPU

bool CPUGPUCSRMatrixMatrixDiagonalMultiply(GPUVector &X, GPUCSRMatrix &A);

//
// GPUGPUCSRMatrixMatrixDiagonalMultiply
// Multiply a digonal matrix specified with a vector with a matrix on GPU

bool GPUGPUCSRMatrixMatrixDiagonalMultiply(GPUVector &X, GPUCSRMatrix &A);

//
// CPUGPUVectorPrepareDiagonalPreconditionerValues
// Prepare diagonal values of the matrix for Diagonal Preconditioner on CPU

bool CPUGPUVectorPrepareDiagonalPreconditionerValues(GPUVector &X);

//
// GPUGPUVectorPrepareDiagonalPreconditionerValues
// Prepare diagonal values of the matrix for Diagonal Preconditioner on GPU

bool GPUGPUVectorPrepareDiagonalPreconditionerValues(GPUVector &X);

//
// GPU_PrepareSPAIPreconditioner
// Prepare SPAI preconditioner on GPU

bool GPU_PrepareSPAIPreconditioner(GPUCSRMatrix &A, GPUCSRMatrix &M);

//
// CPUGPUVectorVectorMultiply
// Vector-Vector multiply on CPU

bool CPUGPUVectorVectorMultiply(GPUVector &X, GPUVector &Y, double &Result);

//
// GPUGPUVectorVectorMultiply
// Vector-Vector multiply on GPU

bool GPUGPUVectorVectorMultiply(GPUVector &X, GPUVector &Y, double &Result);

//
// CPUGPUVectorVectorMultiplyElementWise
// Vector-Vector element-wise multiply on CPU

bool CPUGPUVectorVectorMultiplyElementWise(GPUVector &X, GPUVector &Y,  GPUVector &Z);

//
// GPUGPUVectorVectorMultiplyElementWise
// Vector-Vector element-wise multiply on GPU

bool GPUGPUVectorVectorMultiplyElementWise(GPUVector &X, GPUVector &Y, GPUVector &Z);

//
// CPUGPUVectorNorm2
// Vector norm 2 on CPU

bool CPUGPUVectorNorm2(GPUVector &X, double &Result);

//
// GPUGPUVectorNorm2
// Vector norm 2 on GPU

bool GPUGPUVectorNorm2(GPUVector &X, double &Result);

//
// CPUGPUVectorScaleAndAdd
// Vector scale-and-add on GPU

// Variant 1: Z = A * X + B * Y

bool CPUGPUVectorScaleAndAdd(double A, GPUVector &X, double B, GPUVector &Y, GPUVector &Z);

// Variant 2: Y = A * X + B * Y

bool CPUGPUVectorScaleAndAdd(double A, GPUVector &X, double B, GPUVector &Y);

//
// GPUGPUVectorScaleAndAdd
// Vector scale-and-add on GPU

// Variant 1: Z = A * X + B * Y

bool GPUGPUVectorScaleAndAdd(double A, GPUVector &X, double B, GPUVector &Y, GPUVector &Z);

// Variant 2: Y = A * X + B * Y

bool GPUGPUVectorScaleAndAdd(double A, GPUVector &X, double B, GPUVector &Y);


//Wrapper for cublas dotProduct in double precision
/*double GPU_dotProduct(size_t numElems, const double *firstVec, int incFirstVec,
    const double *secondVec, int incSecondVec);*/

/** ADDED FUNCTIONS AND STRUCTURES **/


//kernel for fill an array with zeros
void GPU_fillWithZeros(size_t numElems, double* gpuVec);

//kernel with VectorScaleAndAdd, adding to Z version
bool GPUGPUVectorScaleAndAdd_addingVersion(double A, GPUVector &X, double B, GPUVector &Y, GPUVector &Z);

//preconditioner classes
void generateResidual(const GPUCSRMatrix& R, const GPUVector& b, const GPUCSRMatrix& A, const GPUVector& u, GPUVector& r);
void calculateInstantVector(GPUVector& u, const GPUVector& b, const GPUCSRMatrix& A, const GPUCSRMatrix& G);
void generateResidualGPUVectorized(const GPUCSRMatrix& R, const GPUVector& b, const GPUCSRMatrix& A, const GPUVector& u, GPUVector& r);
void calculateInstantVectorGPUVectorized(GPUVector& u, const GPUVector& b, const GPUCSRMatrix& A, const GPUCSRMatrix& G);
void multilevel(GPUCSRMatrix**& A, GPUCSRMatrix**& P, GPUCSRMatrix**& R, GPUCSRMatrix**& G, GPUVector& b, GPUVector& u,
			unsigned short lvl, unsigned short maxLevels, size_t* preSweeps, size_t* postSweeps, bool assumeZeros);
size_t generateHierarchy(GPUCSRMatrix**& Matrices, GPUCSRMatrix**& Pmat, GPUCSRMatrix**& Qmat,
        GPUCSRMatrix**& Gmat, double W, size_t numLevelsRoh, size_t max_levels, size_t min_system_size);
void generateP_vCPU(const GPUCSRMatrix& A, const GPUVector& diag, double W, size_t numLevelsRoh, size_t &pNNZ, size_t &pSize1, size_t &pSize2, size_t *&pindices, size_t *&pptr, double *&pvalues);
void generateQ(const GPUCSRMatrix &P, size_t &NNZ, size_t &Size1, size_t &Size2, size_t *&indices, size_t *&ptr, double *&values);

void createPTent(const GPUCSRMatrix &A, size_t &NNZ_ptent, size_t &Size1_ptent, size_t &Size2_ptent, size_t *&ptr_ptent, size_t *&indices_ptent, double *&values_ptent);
void createDiagonal(const GPUCSRMatrix &A, size_t &Size, double *&values);
void createDiagonal_vCPU(const GPUCSRMatrix& A, size_t &NNZ, size_t &Size1, size_t &Size2, size_t *&indices, size_t *&ptr, double *&values, double *&vecDiag);
void computeDenseMatrix(const GPUCSRMatrix& A, double *& vec);
void subIdentityMatrix_cpu(const GPUCSRMatrix& A, size_t* CPU_RowIndices, size_t* CPU_Columns, double* CPU_Values);
void csr_tocsc(const size_t n_row,
	           const size_t n_col,
	           const size_t Ap[],
	           const size_t Aj[],
	           const double Ax[],
	                 size_t*& Bp,
	                 size_t* Bi,
	                 double* Bx);
double checkResidual(const GPUVector& u, const GPUVector& b, const GPUCSRMatrix& A);
double roh(const GPUCSRMatrix& A, size_t iter);
bool checkConvergence(const GPUVector& u, const GPUVector& b, const GPUCSRMatrix& A, const double lastResidual, const double threshold);
void mallocAndCopyMem(double*& CPU, double*& GPU, size_t size);
void mallocAndCopyMem(size_t*& CPU, size_t*& GPU, size_t size);
void malloc_(double*& GPU, size_t size);
void malloc_(size_t*& GPU, size_t size);
template <class Q>
void copyMem(Q*& source, Q*& destiny, size_t size, unsigned short way);

void deletingStuff(size_t* stuff);
void deletingStuff(double* stuff);

void GPUGPUVectorMultiply(double* sourceVec, double* destinyVec, size_t N);

}
}

#endif //KRATOS_GPU_SPARSE_H_INCLUDED
