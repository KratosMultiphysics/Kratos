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
// SPAI preconditioner by Borja Servan

#include "spai_preconditioner.h"

//
// Kernels used in gpu_sparse.cu

//
// GPU_MatrixVectorMultiply_CSR_Kernel
// CSRMatrix Vector multiply kernel

__global__ void GPU_MatrixVectorMultiply_CSR_Kernel(const size_t Rows, const size_t *A_Columns, const size_t *A_RowIndices, const double *A_Values, const double *X_Values, double *Y_Values)
{
	size_t Idx = GlobalIdx();
	
	if (Idx < Rows)
	{
		double YI = static_cast <double> (0);
		
		for (size_t j = A_RowIndices[Idx]; j < A_RowIndices[Idx + 1]; j++)
		
#ifdef USE_TEXTURE_CACHING

			// With caching
			YI += A_Values[j] * Fetch_X(X_Values, A_Columns[j]);

#else

			// Without caching
			YI += A_Values[j] * X_Values[A_Columns[j]];

#endif
			
		Y_Values[Idx] = YI;
	}
}

//
// GPU_MatrixVectorMultiply_CSR_Vectorized_Kernel
// CSRMatrix Vector multiply vectorized kernel



__global__ void GPU_MatrixVectorMultiply_CSR_Vectorized_Kernel(const size_t Rows, const size_t *A_Columns, const size_t *A_RowIndices, const double *A_Values, const double *X_Values, double *Y_Values)
{
	const size_t Idx = GlobalIdx();									// Global thread index
	const size_t Lane = Idx & (HALF_WARP_SIZE - 1);					// Thread index within a half warp
	const size_t HWIdx = Idx >> HALF_WARP_SIZE_BITS;				// Half warp index
	const size_t HWLane = HWIdx & (HALF_WARP_SIZE - 1);				// Half warp lane
	
	__shared__ double Buffer[BLOCK_SIZE];							// Reduction buffer
	__shared__ size_t Limits[BLOCK_SIZE / HALF_WARP_SIZE][2];		// Fetch buffer for upper and lower limits of a row
	
	// We are doing one row per half warp, so the row number is the same as HWIdx
	
	if (HWIdx < Rows)
	{
		Buffer[threadIdx.x] = 0;									// Each thread zeros its own element in the Buffer
		
		if (Lane < 2)
			Limits[HWLane][Lane] = A_RowIndices[HWIdx + Lane];

		const size_t Start = Limits[HWLane][0];
		const size_t End = Limits[HWLane][1];

		for (size_t i = Start + Lane; i < End; i += HALF_WARP_SIZE)

#ifdef USE_TEXTURE_CACHING

			// With caching

			Buffer[threadIdx.x] += A_Values[i] * Fetch_X(X_Values, A_Columns[i]);

#else

			// Without caching
			Buffer[threadIdx.x] += A_Values[i] * X_Values[A_Columns[i]];

#endif

			// Reduce the results in the Buffer; loops are unrolled!
			// There is 16 = 2 ^ 4 threads in a half warp
			
			if (Lane < 8) 
			{
				Buffer[threadIdx.x] += Buffer[threadIdx.x + 8];
				EMUSYNC;
			}
				
			if (Lane < 4)
			{
				Buffer[threadIdx.x] += Buffer[threadIdx.x + 4];
				EMUSYNC;
			}
				
			if (Lane < 2)
			{
				Buffer[threadIdx.x] += Buffer[threadIdx.x + 2];
				EMUSYNC;
			}
			
			if (Lane < 1)
			{
				Buffer[threadIdx.x] += Buffer[threadIdx.x + 1];
				EMUSYNC;
			}
			
			// The first thread in warp has the answer; write it back
			if (Lane == 0)
				 Y_Values[HWIdx] = Buffer[threadIdx.x];
	}
}

//
// GPU_MatrixGetDiagonals_Kernel
// Extract the diagonal elements of a matrix into a vector kernel

__global__ void GPU_MatrixGetDiagonals_Kernel(const size_t Rows, const size_t *A_Columns, const size_t *A_RowIndices, const double *A_Values, double *X_Values)
{
	const size_t Idx = GlobalIdx();
	
	if (Idx < Rows)
	{
		X_Values[Idx] = static_cast <double> (0);
		
		for (size_t j = A_RowIndices[Idx]; j < A_RowIndices[Idx + 1]; j++)
			if (A_Columns[j] == Idx)
				X_Values[Idx] = A_Values[j];
	}
}

//
// GPU_MatrixGetDiagonals_Kernel
// Extract the diagonal elements of a matrix into a vector kernel

__global__ void GPU_MatrixMatrixDiagonalMultiply_Kernel(const size_t Size, const double *X_Values, const size_t *A_Columns, const size_t *A_RowIndices, double *A_Values)
{
	const size_t Idx = GlobalIdx();
	
	if (Idx < Size)
	{
		double t = X_Values[Idx];
	
		for (size_t j = A_RowIndices[Idx]; j < A_RowIndices[Idx + 1]; j++)
			A_Values[j] *= t;
	}
}

//
// GPU_VectorPrepareDiagonalPreconditionerValues_Kernel
// Prepare diagonal values of the matrix for Diagonal Preconditioner kernel

__global__ void GPU_VectorPrepareDiagonalPreconditionerValues_Kernel(const size_t Size, double *X_Values)
{
	const size_t Idx = GlobalIdx();
	
	if (Idx < Size)
		if (X_Values[Idx] == 0.00)
			X_Values[Idx] = 1.00;
			
		else
			X_Values[Idx] = 1.00 / X_Values[Idx];
//			X_Values[Idx] = 1.00 / sqrt(abs(X_Values[Idx]));
}

//
// GPU_VectorVectorMultiplyElementWise_Kernel
// Vector-Vector element-wise multiply kernel

__global__ void GPU_VectorVectorMultiplyElementWise_Kernel(const size_t N, const double *X, const double *Y, double *Z)
{
	const size_t Idx = GlobalIdx();

	if (Idx < N)
		Z[Idx] = X[Idx] * Y[Idx];
}

//
// GPU_VectorScaleAndAdd_1_Kernel
// VectorScaleAndAdd kernel

// Variant 1: Z = A * X + B * Y

__global__ void GPU_VectorScaleAndAdd_1_Kernel(const size_t N, const double A, const double *X, const double B, const double *Y, double *Z)
{
	const size_t Idx = GlobalIdx();

	if (Idx < N)
		Z[Idx] = A * X[Idx] + B * Y[Idx];
}

// Variant 1-A: Z = A * X + B * Y (A = 1.00; B = 1.00)

__global__ void GPU_VectorScaleAndAdd_1_A_Kernel(const size_t N, const double A, const double *X, const double B, const double *Y, double *Z)
{
	const size_t Idx = GlobalIdx();

	if (Idx < N)
		Z[Idx] = X[Idx] + Y[Idx];
}

// Variant 1-B: Z = A * X + B * Y (A = 1.00; B = -1.00)

__global__ void GPU_VectorScaleAndAdd_1_B_Kernel(const size_t N, const double A, const double *X, const double B, const double *Y, double *Z)
{
	const size_t Idx = GlobalIdx();

	if (Idx < N)
		Z[Idx] = X[Idx] - Y[Idx];
}

// Variant 1-C: Z = A * X + B * Y (A = -1.00; B = 1.00)

__global__ void GPU_VectorScaleAndAdd_1_C_Kernel(const size_t N, const double A, const double *X, const double B, const double *Y, double *Z)
{
	const size_t Idx = GlobalIdx();

	if (Idx < N)
		Z[Idx] = -X[Idx] + Y[Idx];
}

// Variant 1-D: Z = A * X + B * Y (A = -1.00; B = -1.00)

__global__ void GPU_VectorScaleAndAdd_1_D_Kernel(const size_t N, const double A, const double *X, const double B, const double *Y, double *Z)
{
	const size_t Idx = GlobalIdx();

	if (Idx < N)
		Z[Idx] = -X[Idx] - Y[Idx];
}

// Variant 1-E: Z = A * X + B * Y (A = 1.00; B = *)

__global__ void GPU_VectorScaleAndAdd_1_E_Kernel(const size_t N, const double A, const double *X, const double B, const double *Y, double *Z)
{
	const size_t Idx = GlobalIdx();

	if (Idx < N)
		Z[Idx] = X[Idx] + B * Y[Idx];
}

// Variant 1-F: Z = A * X + B * Y (A = -1.00; B = *)

__global__ void GPU_VectorScaleAndAdd_1_F_Kernel(const size_t N, const double A, const double *X, const double B, const double *Y, double *Z)
{
	const size_t Idx = GlobalIdx();

	if (Idx < N)
		Z[Idx] = -X[Idx] + B * Y[Idx];
}

// Variant 1-G: Z = A * X + B * Y (A = *; B = 1.00)

__global__ void GPU_VectorScaleAndAdd_1_G_Kernel(const size_t N, const double A, const double *X, const double B, const double *Y, double *Z)
{
	const size_t Idx = GlobalIdx();

	if (Idx < N)
		Z[Idx] = A * X[Idx] + Y[Idx];
}

// Variant 1-H: Z = A * X + B * Y (A = *; B = -1.00)

__global__ void GPU_VectorScaleAndAdd_1_H_Kernel(const size_t N, const double A, const double *X, const double B, const double *Y, double *Z)
{
	const size_t Idx = GlobalIdx();

	if (Idx < N)
		Z[Idx] = A * X[Idx] - Y[Idx];
}

// Variant 2: Y = A * X + B * Y

__global__ void GPU_VectorScaleAndAdd_2_Kernel(const size_t N, const double A, const double *X, const double B, double *Y)
{
	const size_t Idx = GlobalIdx();

	if (Idx < N)
		Y[Idx] = A * X[Idx] + B * Y[Idx];
}

// Variant 2-A: Y = A * X + B * Y (A = 1.00, B = 1.00)

__global__ void GPU_VectorScaleAndAdd_2_A_Kernel(const size_t N, const double A, const double *X, const double B, double *Y)
{
	const size_t Idx = GlobalIdx();

	if (Idx < N)
		Y[Idx] += X[Idx];
}

// Variant 2-B: Y = A * X + B * Y (A = 1.00, B = -1.00)

__global__ void GPU_VectorScaleAndAdd_2_B_Kernel(const size_t N, const double A, const double *X, const double B, double *Y)
{
	const size_t Idx = GlobalIdx();

	if (Idx < N)
		Y[Idx] = X[Idx] - Y[Idx];
}

// Variant 2-C: Y = A * X + B * Y (A = -1.00, B = 1.00)

__global__ void GPU_VectorScaleAndAdd_2_C_Kernel(const size_t N, const double A, const double *X, const double B, double *Y)
{
	const size_t Idx = GlobalIdx();

	if (Idx < N)
		Y[Idx] += -X[Idx];
}

// Variant 2-D: Y = A * X + B * Y (A = -1.00, B = -1.00)

__global__ void GPU_VectorScaleAndAdd_2_D_Kernel(const size_t N, const double A, const double *X, const double B, double *Y)
{
	const size_t Idx = GlobalIdx();

	if (Idx < N)
		Y[Idx] = -X[Idx] - Y[Idx];
}

// Variant 2-E: Y = A * X + B * Y (A = 1.00, B = *)

__global__ void GPU_VectorScaleAndAdd_2_E_Kernel(const size_t N, const double A, const double *X, const double B, double *Y)
{
	const size_t Idx = GlobalIdx();
	if (Idx < N){
		Y[Idx] *= B;
		Y[Idx] += X[Idx];
	}
//	if (Idx < N)
//		Y[Idx] = X[Idx] + B * Y[Idx];
}

// Variant 2-F: Y = A * X + B * Y (A = -1.00, B = *)

__global__ void GPU_VectorScaleAndAdd_2_F_Kernel(const size_t N, const double A, const double *X, const double B, double *Y)
{
	const size_t Idx = GlobalIdx();

	if (Idx < N)
		Y[Idx] = -X[Idx] + B * Y[Idx];
}

// Variant 2-G: Y = A * X + B * Y (A = *, B = 1.00)

__global__ void GPU_VectorScaleAndAdd_2_G_Kernel(const size_t N, const double A, const double *X, const double B, double *Y)
{
	const size_t Idx = GlobalIdx();

	if (Idx < N)
		Y[Idx] += A * X[Idx];
}

// Variant 2-H: Y = A * X + B * Y (A = *, B = -1.00)

__global__ void GPU_VectorScaleAndAdd_2_H_Kernel(const size_t N, const double A, const double *X, const double B, double *Y)
{
	const size_t Idx = GlobalIdx();

	if (Idx < N)
		Y[Idx] = A * X[Idx] - Y[Idx];
}

	/** ADDED KERNELS  **/
/** FOR GPU_SPARSE **/

//
// GPU_VectorScaleAndAdd_1_Kernel ADDING VERSIONS
// VectorScaleAndAdd kernel

// Variant 1: Z = A * X + B * Y

__global__ void GPU_VectorScaleAndAdd_1_Kernel_addingVersion(const size_t N, const double A, const double *X, const double B, const double *Y, double *Z)
{
	const size_t Idx = GlobalIdx();

	if (Idx < N)
		Z[Idx] += A * X[Idx] + B * Y[Idx];
}

// Variant 1-A: Z = A * X + B * Y (A = 1.00; B = 1.00)

__global__ void GPU_VectorScaleAndAdd_1_A_Kernel_addingVersion(const size_t N, const double A, const double *X, const double B, const double *Y, double *Z)
{
	const size_t Idx = GlobalIdx();

	if (Idx < N)
		Z[Idx] += X[Idx] + Y[Idx];
}

// Variant 1-B: Z = A * X + B * Y (A = 1.00; B = -1.00)

__global__ void GPU_VectorScaleAndAdd_1_B_Kernel_addingVersion(const size_t N, const double A, const double *X, const double B, const double *Y, double *Z)
{
	const size_t Idx = GlobalIdx();

	if (Idx < N)
		Z[Idx] += X[Idx] - Y[Idx];
}

// Variant 1-C: Z = A * X + B * Y (A = -1.00; B = 1.00)

__global__ void GPU_VectorScaleAndAdd_1_C_Kernel_addingVersion(const size_t N, const double A, const double *X, const double B, const double *Y, double *Z)
{
	const size_t Idx = GlobalIdx();

	if (Idx < N)
		Z[Idx] += -X[Idx] + Y[Idx];
}

// Variant 1-D: Z = A * X + B * Y (A = -1.00; B = -1.00)

__global__ void GPU_VectorScaleAndAdd_1_D_Kernel_addingVersion(const size_t N, const double A, const double *X, const double B, const double *Y, double *Z)
{
	const size_t Idx = GlobalIdx();

	if (Idx < N)
		Z[Idx] += -X[Idx] - Y[Idx];
}

// Variant 1-E: Z = A * X + B * Y (A = 1.00; B = *)

__global__ void GPU_VectorScaleAndAdd_1_E_Kernel_addingVersion(const size_t N, const double A, const double *X, const double B, const double *Y, double *Z)
{
	const size_t Idx = GlobalIdx();

	if (Idx < N)
		Z[Idx] += X[Idx] + B * Y[Idx];
}

// Variant 1-F: Z = A * X + B * Y (A = -1.00; B = *)

__global__ void GPU_VectorScaleAndAdd_1_F_Kernel_addingVersion(const size_t N, const double A, const double *X, const double B, const double *Y, double *Z)
{
	const size_t Idx = GlobalIdx();

	if (Idx < N)
		Z[Idx] += -X[Idx] + B * Y[Idx];
}

// Variant 1-G: Z = A * X + B * Y (A = *; B = 1.00)

__global__ void GPU_VectorScaleAndAdd_1_G_Kernel_addingVersion(const size_t N, const double A, const double *X, const double B, const double *Y, double *Z)
{
	const size_t Idx = GlobalIdx();

	if (Idx < N)
		Z[Idx] += A * X[Idx] + Y[Idx];
}

// Variant 1-H: Z = A * X + B * Y (A = *; B = -1.00)

__global__ void GPU_VectorScaleAndAdd_1_H_Kernel_addingVersion(const size_t N, const double A, const double *X, const double B, const double *Y, double *Z)
{
	const size_t Idx = GlobalIdx();

	if (Idx < N)
		Z[Idx] += A * X[Idx] - Y[Idx];
}

/** FOR PRECONDITIONER **/
/** Adding version of MatrixVector multiply **/
__global__ void GPU_MatrixVectorMultiply_CSR_Kernel_addingVersion(const size_t Rows, const size_t *A_Columns, const size_t *A_RowIndices, const double *A_Values, const double *X_Values, double *Y_Values)
{
	size_t Idx = GlobalIdx();

	if (Idx < Rows)
	{
		double YI = static_cast <double> (0);

		for (size_t j = A_RowIndices[Idx]; j < A_RowIndices[Idx + 1]; j++)

#ifdef USE_TEXTURE_CACHING

			// With caching
			YI += A_Values[j] * Fetch_X(X_Values, A_Columns[j]);

#else

			// Without caching
			YI += A_Values[j] * X_Values[A_Columns[j]];

#endif

		Y_Values[Idx] += YI;
	}
}

/** Kernel for fill a vector with zeros **/
__global__ void fillWithZeros (double* values, size_t rows){
    size_t idx = GlobalIdx();
    if(idx < rows)
        values[idx] = 0.0;
}

/** Kernel to add 2 vectors
 *  version that let the result on first vector **/
__global__ void sumVectorVector(double* source, double* dest, size_t size){
    size_t idx = GlobalIdx();
    if(idx < size){
        dest[idx] = dest[idx] + source[idx];
    }
}

        /** NEW KERNELS **/
/** kernel to subtract source of dest **/
template < class G >
__global__ void subVectorVector(G* source, G* dest, size_t size){
    size_t idx = GlobalIdx();
    if(idx < size){
        dest[idx] -= source[idx];
    }
}
/** kernel to subtract source2 from source1 and save result on dest **/
template < class G >
__global__ void subVectorVector(G* source1, G* source2, G* dest, size_t size){
    size_t idx = GlobalIdx();
    if(idx < size){
        dest[idx] = source1[idx] - source2[idx];
    }
}

/** kernel to subtract a value from each elem of G **/
template < class G >
__global__ void subVectorConstantValue (G* subVector, G num, G* sourceVector, size_t numElems){
    size_t idx = GlobalIdx();
    if(idx < numElems){
        subVector[idx] -= num * sourceVector[idx];
    }
}


/** kernel to divide a value from each elem of G **/
template < class G >
__global__ void divideVectorConstantValue (G* vector, G elem, size_t numElems){
    size_t idx = GlobalIdx();
    if(idx < numElems){
        vector[idx] /= elem;
    }
}


