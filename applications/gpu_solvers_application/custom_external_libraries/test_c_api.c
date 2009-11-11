
#include "cwrapper.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <stdlib.h>
#include <stdio.h>


void selectB(size_t *ptr, size_t* indices, double* values, double* bValues, size_t rowChoice, size_t numRows){
	size_t currentPointer = 0;
	int i;
	for(i = ptr[rowChoice]; i < ptr[rowChoice+1]; i++){
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
	char* name = "mat8.mm";
	char* name_vec = "vecb8.mm";
	

	loadDataFromFile(name, &indices, &ptr, &values, &numRows, &numCols, &NNZ);
	

	double* x;
	double* b;
	size_t b_rows;
	size_t b_cols;
//	b = allocate_double(numRows);
	x = allocate_double(numRows);
	loadDataFromFile_vector(name_vec, &b,
		&b_rows, &b_cols);
	//b = A[rowChoice]; b is copied from the row number rowChoice of A
/*	size_t rowChoice = 2;
	selectB(ptr, indices, values, b, rowChoice, numRows);	

	size_t cn;
	for(cn = 0; cn < numRows; cn++){
		printf("%lg, ", b[cn]);
	}
	printf("\n\n");*/
	
	//Solving strategies	
		//BICGSTAB without preconditioner


	setZeros(x, numRows);
	Solve(numRows, numCols, NNZ, values, indices, ptr, numRows, x, b, 1e-9, 100, &mBNorm, &mResidualNorm, &mIterationsNumber, GPU_BICGSTAB, NOPRECOND);
	printf("Finished BICGSTAB without preconditioner with norm(b) = %lg, norm(residual) = %lg, numIterationsPerformed = %u\n\n", mBNorm, mResidualNorm, mIterationsNumber);

		//BICGSTAB with diagonal preconditioner
	setZeros(x, numRows);
	Solve(numRows, numCols, NNZ, values, indices, ptr, numRows, x, b, 1e-9, 100, &mBNorm, &mResidualNorm, &mIterationsNumber, GPU_BICGSTAB, DIAGONAL);
	printf("Finished BICGSTAB with diagonal preconditioner with norm(b) = %lg, norm(residual) = %lg, numIterationsPerformed = %u\n\n", mBNorm, mResidualNorm, mIterationsNumber);
	
		//CG without preconditioner
	setZeros(x, numRows);
	Solve(numRows, numCols, NNZ, values, indices, ptr, numRows, x, b, 1e-9, 100, &mBNorm, &mResidualNorm, &mIterationsNumber, GPU_CG, NOPRECOND);
	printf("Finished CG without preconditioner with norm(b) = %lg, norm(residual) = %lg, numIterationsPerformed = %u\n\n", mBNorm, mResidualNorm, mIterationsNumber);
		
		//CG with diagonal preconditioner
	setZeros(x, numRows);
	Solve(numRows, numCols, NNZ, values, indices, ptr, numRows, x, b, 1e-9, 100, &mBNorm, &mResidualNorm, &mIterationsNumber, GPU_CG, DIAGONAL);
	printf("Finished CG with diagonal preconditioner with norm(b) = %lg, norm(residual) = %lg, numIterationsPerformed = %u\n\n", mBNorm, mResidualNorm, mIterationsNumber);
		
		//CG with amg preconditioner
	setZeros(x, numRows);
	Solve(numRows, numCols, NNZ, values, indices, ptr, numRows, x, b, 1e-9, 100, &mBNorm, &mResidualNorm, &mIterationsNumber, GPU_CG, AMG);
	printf("Finished CG with amg preconditioner with norm(b) = %lg, norm(residual) = %lg, numIterationsPerformed = %u\n\n", mBNorm, mResidualNorm, mIterationsNumber);

	//Cleaning resources	
	
	deallocate_double(b);
	deallocate_double(x);
	deallocate_size_t(ptr);
	deallocate_size_t(indices);
	deallocate_double(values);
	return 0;
}
