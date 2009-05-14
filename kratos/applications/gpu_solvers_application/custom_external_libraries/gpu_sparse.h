#include <stddef.h>

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
	bool CopyGPUValuesFrom(GPUVector &V);
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

	GPUCSRMatrix(size_t _NNZ, size_t _Size1, size_t _Size2, size_t *_CPU_Columns, size_t *_CPU_RowIndices, double *_CPU_Values);
	~GPUCSRMatrix();

	bool GPU_Allocate();
	bool GPU_Free();
	
	bool Copy(CopyDirection Direction, bool CopyValuesOnly);
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
// CPU_VectorVectorMultiply
// Vector-Vector multiply on CPU

bool CPU_VectorVectorMultiply(GPUVector &X, GPUVector &Y, double &Result);

//
// GPU_VectorVectorMultiply
// Vector-Vector multiply on GPU

bool GPU_VectorVectorMultiply(GPUVector &X, GPUVector &Y, double &Result);

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

}

}
