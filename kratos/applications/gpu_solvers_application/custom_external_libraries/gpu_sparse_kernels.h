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

// WARNING! THIS HAS BUGS! DO NOT USE!

__global__ void GPU_MatrixVectorMultiply_CSR_Vectorized_Kernel(const size_t Rows, const size_t *A_Columns, const size_t *A_RowIndices, const double *A_Values, const double *X_Values, double *Y_Values)
{
	size_t Idx = GlobalIdx();				// Global thread index
	size_t WIdx = Idx / WARP_SIZE;			// There are 32 threads in a warp
	size_t Lane = Idx & (WARP_SIZE - 1);	// Thread index within a warp
	
	__shared__ double Tree[BLOCK_SIZE];
	
	// We are doing one row per warp, so the row number is the same as WarpIdx
	
	if (WIdx < Rows)
	{
		Tree[threadIdx.x] = 0;		// Each thread zeros its own element in the tree
		
		size_t Start = A_RowIndices[WIdx];
		size_t End = A_RowIndices[WIdx + 1];
		
		for (size_t i = Start + Lane; i < End; i += WARP_SIZE)

#ifdef USE_TEXTURE_CACHING

			// With caching

			Tree[threadIdx.x] += A_Values[i] * Fetch_X(X_Values, A_Columns[i]);

#else

			// Without caching
			Tree[threadIdx.x] += A_Values[i] * X_Values[A_Columns[i]];

#endif

			// Reduce the results in the Tree; loops are unrolled!
			// Note 1: There is 32 = 2 ^ 5 threads in a warp
			// Note 2: This code relies on implicit synchronization of threads in a warp
			
			if (Lane < 16)
				Tree[threadIdx.x] += Tree[threadIdx.x + 16];
				
			if (Lane < 8) 
				Tree[threadIdx.x] += Tree[threadIdx.x + 8];
				
			if (Lane < 4)
				Tree[threadIdx.x] += Tree[threadIdx.x + 4];
				
			if (Lane < 2)
				Tree[threadIdx.x] += Tree[threadIdx.x + 2];
			
			if (Lane < 1)
				Tree[threadIdx.x] += Tree[threadIdx.x + 1];
			
			// The first thread in warp has the answer; write it back
			if (Lane == 0)
				 Y_Values[WIdx] = Tree[threadIdx.x];
	}
}

//
// GPU_VectorScaleAndAdd_1_Kernel
// VectorScaleAndAdd kernel

// Variant 1: Z = A * X + B * Y

__global__ void GPU_VectorScaleAndAdd_1_Kernel(const size_t N, const double A, const double *X, const double B, const double *Y, double *Z)
{
	size_t Idx = GlobalIdx();

	if (Idx < N)
		Z[Idx] = A * X[Idx] + B * Y[Idx];
}

// Variant 1-A: Z = A * X + B * Y (A = 1.00; B = 1.00)

__global__ void GPU_VectorScaleAndAdd_1_A_Kernel(const size_t N, const double A, const double *X, const double B, const double *Y, double *Z)
{
	size_t Idx = GlobalIdx();

	if (Idx < N)
		Z[Idx] = X[Idx] + Y[Idx];
}

// Variant 1-B: Z = A * X + B * Y (A = 1.00; B = -1.00)

__global__ void GPU_VectorScaleAndAdd_1_B_Kernel(const size_t N, const double A, const double *X, const double B, const double *Y, double *Z)
{
	size_t Idx = GlobalIdx();

	if (Idx < N)
		Z[Idx] = X[Idx] - Y[Idx];
}

// Variant 1-C: Z = A * X + B * Y (A = -1.00; B = 1.00)

__global__ void GPU_VectorScaleAndAdd_1_C_Kernel(const size_t N, const double A, const double *X, const double B, const double *Y, double *Z)
{
	size_t Idx = GlobalIdx();

	if (Idx < N)
		Z[Idx] = -X[Idx] + Y[Idx];
}

// Variant 1-D: Z = A * X + B * Y (A = -1.00; B = -1.00)

__global__ void GPU_VectorScaleAndAdd_1_D_Kernel(const size_t N, const double A, const double *X, const double B, const double *Y, double *Z)
{
	size_t Idx = GlobalIdx();

	if (Idx < N)
		Z[Idx] = -X[Idx] - Y[Idx];
}

// Variant 1-E: Z = A * X + B * Y (A = 1.00; B = *)

__global__ void GPU_VectorScaleAndAdd_1_E_Kernel(const size_t N, const double A, const double *X, const double B, const double *Y, double *Z)
{
	size_t Idx = GlobalIdx();

	if (Idx < N)
		Z[Idx] = X[Idx] + B * Y[Idx];
}

// Variant 1-F: Z = A * X + B * Y (A = -1.00; B = *)

__global__ void GPU_VectorScaleAndAdd_1_F_Kernel(const size_t N, const double A, const double *X, const double B, const double *Y, double *Z)
{
	size_t Idx = GlobalIdx();

	if (Idx < N)
		Z[Idx] = -X[Idx] + B * Y[Idx];
}

// Variant 1-G: Z = A * X + B * Y (A = *; B = 1.00)

__global__ void GPU_VectorScaleAndAdd_1_G_Kernel(const size_t N, const double A, const double *X, const double B, const double *Y, double *Z)
{
	size_t Idx = GlobalIdx();

	if (Idx < N)
		Z[Idx] = A * X[Idx] + Y[Idx];
}

// Variant 1-H: Z = A * X + B * Y (A = *; B = -1.00)

__global__ void GPU_VectorScaleAndAdd_1_H_Kernel(const size_t N, const double A, const double *X, const double B, const double *Y, double *Z)
{
	size_t Idx = GlobalIdx();

	if (Idx < N)
		Z[Idx] = A * X[Idx] - Y[Idx];
}

// Variant 2: Y = A * X + B * Y

__global__ void GPU_VectorScaleAndAdd_2_Kernel(const size_t N, const double A, const double *X, const double B, double *Y)
{
	size_t Idx = GlobalIdx();

	if (Idx < N)
		Y[Idx] = A * X[Idx] + B * Y[Idx];
}

// Variant 2-A: Y = A * X + B * Y (A = 1.00, B = 1.00)

__global__ void GPU_VectorScaleAndAdd_2_A_Kernel(const size_t N, const double A, const double *X, const double B, double *Y)
{
	size_t Idx = GlobalIdx();

	if (Idx < N)
		Y[Idx] += X[Idx];
}

// Variant 2-B: Y = A * X + B * Y (A = 1.00, B = -1.00)

__global__ void GPU_VectorScaleAndAdd_2_B_Kernel(const size_t N, const double A, const double *X, const double B, double *Y)
{
	size_t Idx = GlobalIdx();

	if (Idx < N)
		Y[Idx] = X[Idx] - Y[Idx];
}

// Variant 2-C: Y = A * X + B * Y (A = -1.00, B = 1.00)

__global__ void GPU_VectorScaleAndAdd_2_C_Kernel(const size_t N, const double A, const double *X, const double B, double *Y)
{
	size_t Idx = GlobalIdx();

	if (Idx < N)
		Y[Idx] += -X[Idx];
}

// Variant 2-D: Y = A * X + B * Y (A = -1.00, B = -1.00)

__global__ void GPU_VectorScaleAndAdd_2_D_Kernel(const size_t N, const double A, const double *X, const double B, double *Y)
{
	size_t Idx = GlobalIdx();

	if (Idx < N)
		Y[Idx] = -X[Idx] - Y[Idx];
}

// Variant 2-E: Y = A * X + B * Y (A = 1.00, B = *)

__global__ void GPU_VectorScaleAndAdd_2_E_Kernel(const size_t N, const double A, const double *X, const double B, double *Y)
{
	size_t Idx = GlobalIdx();

	if (Idx < N)
		Y[Idx] = X[Idx] + B * Y[Idx];
}

// Variant 2-F: Y = A * X + B * Y (A = -1.00, B = *)

__global__ void GPU_VectorScaleAndAdd_2_F_Kernel(const size_t N, const double A, const double *X, const double B, double *Y)
{
	size_t Idx = GlobalIdx();

	if (Idx < N)
		Y[Idx] = -X[Idx] + B * Y[Idx];
}

// Variant 2-G: Y = A * X + B * Y (A = *, B = 1.00)

__global__ void GPU_VectorScaleAndAdd_2_G_Kernel(const size_t N, const double A, const double *X, const double B, double *Y)
{
	size_t Idx = GlobalIdx();

	if (Idx < N)
		Y[Idx] += A * X[Idx];
}

// Variant 2-H: Y = A * X + B * Y (A = *, B = -1.00)

__global__ void GPU_VectorScaleAndAdd_2_H_Kernel(const size_t N, const double A, const double *X, const double B, double *Y)
{
	size_t Idx = GlobalIdx();

	if (Idx < N)
		Y[Idx] = A * X[Idx] - Y[Idx];
}

