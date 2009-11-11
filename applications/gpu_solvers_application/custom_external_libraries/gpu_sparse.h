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

	bool Allocated;

	GPUCSRMatrix(size_t _NNZ, size_t _Size1, size_t _Size2, size_t *_CPU_Columns, size_t *_CPU_RowIndices, double *_CPU_Values, bool _NZMultiple16 = true);
	GPUCSRMatrix(size_t _NNZ, size_t _Size1, size_t _Size2);
	~GPUCSRMatrix();

	bool GPU_Allocate();
	bool GPU_Free();

	bool Copy(CopyDirection Direction, bool CopyValuesOnly);
	bool CopyFromGPU(GPUCSRMatrix &M, bool CopyStructure, bool CopyValues);
};

// Operations defined on GPUCSRMatrix and GPUVector

//
// CPU_MatrixVectorMultiply
// Matrix-Vector multiply on CPU

bool CPU_MatrixVectorMultiply(GPUCSRMatrix &A, GPUVector &X, GPUVector &Y);

//
// GPU_MatrixVectorMultiply
// Matrix-Vector multiply on GPU

bool GPU_MatrixVectorMultiply(GPUCSRMatrix &A, GPUVector &X, GPUVector &Y);

//
// CPU_MatrixGetDiagonals
// Extract the diagonal elements of a matrix into a vector on CPU

bool CPU_MatrixGetDiagonals(GPUCSRMatrix &A, GPUVector &X);

//
// GPU_MatrixGetDiagonals
// Extract the diagonal elements of a matrix into a vector on GPU

bool GPU_MatrixGetDiagonals(GPUCSRMatrix &A, GPUVector &X);

//
// CPU_MatrixMatrixDiagonalMultiply
// Multiply a digonal matrix specified with a vector with a matrix on CPU

bool CPU_MatrixMatrixDiagonalMultiply(GPUVector &X, GPUCSRMatrix &A);

//
// GPU_MatrixMatrixDiagonalMultiply
// Multiply a digonal matrix specified with a vector with a matrix on GPU

bool GPU_MatrixMatrixDiagonalMultiply(GPUVector &X, GPUCSRMatrix &A);

//
// CPU_VectorPrepareDiagonalPreconditionerValues
// Prepare diagonal values of the matrix for Diagonal Preconditioner on CPU

bool CPU_VectorPrepareDiagonalPreconditionerValues(GPUVector &X);

//
// GPU_VectorPrepareDiagonalPreconditionerValues
// Prepare diagonal values of the matrix for Diagonal Preconditioner on GPU

bool GPU_VectorPrepareDiagonalPreconditionerValues(GPUVector &X);

//
// GPU_PrepareSPAIPreconditioner
// Prepare SPAI preconditioner on GPU

bool GPU_PrepareSPAIPreconditioner(GPUCSRMatrix &A, GPUCSRMatrix &M);

//
// CPU_VectorVectorMultiply
// Vector-Vector multiply on CPU

bool CPU_VectorVectorMultiply(GPUVector &X, GPUVector &Y, double &Result);

//
// GPU_VectorVectorMultiply
// Vector-Vector multiply on GPU

bool GPU_VectorVectorMultiply(GPUVector &X, GPUVector &Y, double &Result);

//
// CPU_VectorVectorMultiplyElementWise
// Vector-Vector element-wise multiply on CPU

bool CPU_VectorVectorMultiplyElementWise(GPUVector &X, GPUVector &Y,  GPUVector &Z);

//
// GPU_VectorVectorMultiplyElementWise
// Vector-Vector element-wise multiply on GPU

bool GPU_VectorVectorMultiplyElementWise(GPUVector &X, GPUVector &Y, GPUVector &Z);

//
// CPU_VectorNorm2
// Vector norm 2 on CPU

bool CPU_VectorNorm2(GPUVector &X, double &Result);

//
// GPU_VectorNorm2
// Vector norm 2 on GPU

bool GPU_VectorNorm2(GPUVector &X, double &Result);

//
// CPU_VectorScaleAndAdd
// Vector scale-and-add on GPU

// Variant 1: Z = A * X + B * Y

bool CPU_VectorScaleAndAdd(double A, GPUVector &X, double B, GPUVector &Y, GPUVector &Z);

// Variant 2: Y = A * X + B * Y

bool CPU_VectorScaleAndAdd(double A, GPUVector &X, double B, GPUVector &Y);

//
// GPU_VectorScaleAndAdd
// Vector scale-and-add on GPU

// Variant 1: Z = A * X + B * Y

bool GPU_VectorScaleAndAdd(double A, GPUVector &X, double B, GPUVector &Y, GPUVector &Z);

// Variant 2: Y = A * X + B * Y

bool GPU_VectorScaleAndAdd(double A, GPUVector &X, double B, GPUVector &Y);


//Wrapper for cublas dotProduct in double precision
/*double GPU_dotProduct(size_t numElems, const double *firstVec, int incFirstVec,
    const double *secondVec, int incSecondVec);*/

/** ADDED FUNCTIONS AND STRUCTURES **/

typedef struct{
    size_t numRows;
    size_t numCols;
    size_t numNNZ;
    size_t* ptr_cpu;
    size_t* ptr_gpu;
    size_t* indices_cpu;
    size_t* indices_gpu;

    double* values_cpu;
    double* values_gpu;

    double* matAuxValues;
    size_t numValuesDenseRep;
}_Matrix;

typedef struct{
    size_t numElems;
    double* values_cpu;
    double* values_gpu;
}_Vector;


//kernel for fill an array with zeros
void GPU_fillWithZeros(size_t numElems, double* gpuVec);

//kernel with VectorScaleAndAdd, adding to Z version
bool GPU_VectorScaleAndAdd_addingVersion(double A, GPUVector &X, double B, GPUVector &Y, GPUVector &Z);

//preconditioner classes
void generateResidual(const _Matrix& R, const _Vector& b, const _Matrix& A, const _Vector& u, _Vector& r);
void calculateInstantVector(_Vector& u, const _Vector& b, const _Matrix& A, const _Matrix& G);
void generateResidual_vectorized(const _Matrix& R, const _Vector& b, const _Matrix& A, const _Vector& u, _Vector& r);
void calculateInstantVector_vectorized(_Vector& u, const _Vector& b, const _Matrix& A, const _Matrix& G);
void multilevel(_Matrix*& A, _Matrix*& P, _Matrix*& R, _Matrix*& G, _Vector& b, _Vector& u,
			unsigned short lvl, unsigned short maxLevels, size_t* preSweeps, size_t* postSweeps, bool assumeZeros);
size_t generateHierarchy(_Matrix*& Matrices, _Matrix*& Pmat, _Matrix*& Qmat,
        _Matrix*& Gmat, double W, size_t numLevelsRoh, size_t max_levels, size_t min_system_size);
_Matrix generateP_vCPU(const _Matrix& A, const _Vector& diag, double W, size_t numLevelsRoh);
_Matrix generateQ(const _Matrix& P);

void createPTent(const _Matrix& A, _Matrix& P);
void createDiagonal(const _Matrix& A, _Vector& G);
void createDiagonal_vCPU(const _Matrix& A, _Matrix& G, _Vector& res);
void computeDenseMatrix(const _Matrix& A, double *& vec);
void subIdentityMatrix_cpu(const _Matrix& A, _Matrix& sub);
void csr_tocsc(const size_t n_row,
	           const size_t n_col,
	           const size_t Ap[],
	           const size_t Aj[],
	           const double Ax[],
	                 size_t*& Bp,
	                 size_t* Bi,
	                 double* Bx);
double checkResidual(const _Vector& u, const _Vector& b, const _Matrix& A);
double roh(const _Matrix& A, size_t iter);
bool checkConvergence(const _Vector& u, const _Vector& b, const _Matrix& A, const double lastResidual, const double threshold);
void mallocAndCopyMem(double*& CPU, double*& GPU, size_t size);
void mallocAndCopyMem(size_t*& CPU, size_t*& GPU, size_t size);
void malloc_(double*& GPU, size_t size);
void malloc_(size_t*& GPU, size_t size);
template <class Q>
void copyMem(Q*& source, Q*& destiny, size_t size, unsigned short way);

void deletingStuff(size_t* stuff);
void deletingStuff(double* stuff);

void GPU_VectorMultiply(double* sourceVec, double* destinyVec, size_t N);

}
}

#endif //KRATOS_GPU_SPARSE_H_INCLUDED
