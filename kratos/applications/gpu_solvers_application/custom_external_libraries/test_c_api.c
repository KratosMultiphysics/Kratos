
#include "cwrapper.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

void selectB(size_t *ptr, size_t* indices, double* values, double* bValues, size_t rowChoice, size_t numRows){
	size_t currentPointer = 0;
	int i;
	for(i = (int)ptr[rowChoice]; i < (int)ptr[rowChoice+1]; i++){
		if(indices[i] <= currentPointer){
			bValues[currentPointer] = values[i];
		}else{
			bValues[currentPointer] = 0.0;
			i--;
		}	
		currentPointer++;
	}
	while(currentPointer < numRows){
		bValues[currentPointer] = 0.0;
		currentPointer++;
	}
}

void setZeros(double* x, size_t numRows){
	size_t i = 0;
	for(i = 0; i < numRows; i++){
		x[i] = 0.0;
	}
}


int main(){
	size_t numRows, NNZ, numCols;
	size_t *ptr, *indices;
	double *values;
	double mBNorm, mResidualNorm;
	size_t mIterationsNumber;
	//FillMatrix_CPU(ptr, indices, values, NNZ, numRows);
//	char* name = "mat8.mm";
//	char* name_vec = "vecb8.mm";

//	char* name = "examples/aa.mm";
//	char* name_vec = "examples/bb.mm";

//	char* name = "data/mat90.mm";
//	char* name_vec = "data/vecb90.mm";

	char* name = "examples/mat40.mm";
	char* name_vec = "examples/vecb40.mm";
	

	double tolerance = 1e-9;
	size_t maxIters = 5000;

	double* x;
	double* b;
	size_t b_rows;
	size_t b_cols;

	clock_t in, out;

	loadDataFromFile_vector(name_vec, &b, &b_rows, &b_cols);	
	loadDataFromFile(name, &indices, &ptr, &values, &numRows, &numCols, &NNZ);
	
	x = allocate_double(numRows);
	

	//Solving strategies	
		//BICGSTAB without preconditioner
	setZeros(x, numRows);
	printf("	STARTING TESTS\n");
	in = clock();
	Solve(numRows, numCols, NNZ, values, indices, ptr, numRows, x, b, tolerance, maxIters, &mBNorm, &mResidualNorm, &mIterationsNumber, GPU_BICGSTAB, NOPRECOND);
	out = clock();
	printf("Finished BICGSTAB without preconditioner with norm(b) = %lg, norm(residual) = %lg, numIterationsPerformed = %u\n", mBNorm, mResidualNorm, mIterationsNumber);
	printf("Timing for this function = %lg\n\n", (((double)out-in)/CLOCKS_PER_SEC));
	

		//BICGSTAB with diagonal preconditioner
	setZeros(x, numRows);
	in = clock();
	Solve(numRows, numCols, NNZ, values, indices, ptr, numRows, x, b, tolerance, maxIters, &mBNorm, &mResidualNorm, &mIterationsNumber, GPU_BICGSTAB, DIAGONAL);
	out = clock();
	printf("Finished BICGSTAB with diagonal preconditioner with norm(b) = %lg, norm(residual) = %lg, numIterationsPerformed = %u\n", mBNorm, mResidualNorm, mIterationsNumber);
	printf("Timing for this function = %lg\n\n", (((double)out-in)/CLOCKS_PER_SEC));
	
		//CG without preconditioner
	setZeros(x, numRows);
	in = clock();
	Solve(numRows, numCols, NNZ, values, indices, ptr, numRows, x, b, tolerance, maxIters, &mBNorm, &mResidualNorm, &mIterationsNumber, GPU_CG, NOPRECOND);
	out = clock();
	printf("Finished CG without preconditioner with norm(b) = %lg, norm(residual) = %lg, numIterationsPerformed = %u\n", mBNorm, mResidualNorm, mIterationsNumber);
	printf("Timing for this function = %lg\n\n", (((double)out-in)/CLOCKS_PER_SEC));
		
		//CG with diagonal preconditioner
	setZeros(x, numRows);
	in = clock();
	Solve(numRows, numCols, NNZ, values, indices, ptr, numRows, x, b, tolerance, maxIters, &mBNorm, &mResidualNorm, &mIterationsNumber, GPU_CG, DIAGONAL);
	out = clock();
	printf("Finished CG with diagonal preconditioner with norm(b) = %lg, norm(residual) = %lg, numIterationsPerformed = %u\n", mBNorm, mResidualNorm, mIterationsNumber);
	printf("Timing for this function = %lg\n\n", (((double)out-in)/CLOCKS_PER_SEC));
		
		//CG with amg preconditioner
	setZeros(x, numRows);
	in = clock();
	Solve(numRows, numCols, NNZ, values, indices, ptr, numRows, x, b, tolerance, maxIters, &mBNorm, &mResidualNorm, &mIterationsNumber, GPU_CG, AMG);
	out = clock();
	printf("Finished CG with amg preconditioner with norm(b) = %lg, norm(residual) = %lg, numIterationsPerformed = %u\n", mBNorm, mResidualNorm, mIterationsNumber);
	printf("Timing for this function = %lg\n\n", (((double)out-in)/CLOCKS_PER_SEC));

	//Cleaning resources	
	
	deallocate_double(b);
	deallocate_double(x);
	deallocate_size_t(ptr);
	deallocate_size_t(indices);
	deallocate_double(values);
	return 0;
}
