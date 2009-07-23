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


#if !defined(KRATOS_MATRIX_MARKET_INTERFACE_H_INCLUDED )
#define  KRATOS_MATRIX_MARKET_INTERFACE_H_INCLUDED


// System includes
#include <stdio.h>


// External includes

// To avoid linking problems
extern "C"
{ 
	#include "includes/mmio.h"
}

#include <boost/numeric/ublas/matrix_sparse.hpp>


// Project includes
#include "includes/ublas_interface.h"


namespace Kratos {

	// Matrix I/O routines
	
	static bool ReadMatrixMarketMatrix(char *FileName, CompressedMatrix &M)
	{
		// Open MM file for reading
		FILE *f = fopen(FileName, "r");
	
		if (f == NULL)
		{
			printf("ReadMatrixMarketMatrix(): unable to open %s.\n", FileName);
			return false;
		}

		// Process MM file header
		MM_typecode mm_code;
	
		if (mm_read_banner(f, &mm_code) != 0)
		{
			printf("ReadMatrixMarketMatrix(): unable to read MatrixMarket banner.\n");
			fclose(f);
			return false;
		}
	
		if (!mm_is_valid(mm_code))
		{
			printf("ReadMatrixMarketMatrix(): invalid MatrixMarket banner.\n");
			fclose(f);
			return false;
		}
	
		// Check for supported types of MM file
		if (!((mm_is_real(mm_code) || mm_is_integer(mm_code) || mm_is_pattern(mm_code)) && mm_is_coordinate(mm_code) && mm_is_sparse(mm_code)))
		{
			printf("ReadMatrixMarketMatrix(): invalid MatrixMarket type, \"%s\".\n",  mm_typecode_to_str(mm_code));
			fclose(f);
			return false;
		}
	
		// Read MM dimensions and NNZ
		int size1, size2, nnz;
	
		if (mm_read_mtx_crd_size(f, &size1, &size2, &nnz) != 0)
		{
			printf("ReadMatrixMarketMatrix(): cannot read dimensions and NNZ.\n");
			fclose(f);
			return false;
		}
	
		// Allocate temporary arrays
		int *I = new int[nnz];
		int *J = new int[nnz];
		double *V = new double[nnz];
	
		// Read MM file
	
		// Pattern file, only non-zero structure
		if (mm_is_pattern(mm_code))
			for (int i = 0; i < nnz; i++)
			{
				if (fscanf(f, "%d %d", &I[i], &J[i]) != 2)
				{
					printf("ReadMatrixMarketMatrix(): invalid data.\n");
					fclose(f);					
				
					delete[] I;
					delete[] J;
					delete[] V;
				
					return false;
				}
			
				// Adjust to 0-based
				I[i]--;
				J[i]--;
			
				// Set all values to 1.00
				V[i] = 1.00;
			}
		else		
			for (int i = 0; i < nnz; i++)
			{
				if (fscanf(f, "%d %d %lg", &I[i], &J[i], &V[i]) != 3)
				{
					printf("ReadMatrixMarketMatrix(): invalid data.\n");
					fclose(f);
				
					delete[] I;
					delete[] J;
					delete[] V;
				
					return false;
				}
			
				// Adjust to 0-based
				I[i]--;
				J[i]--;
			}
	
		fclose(f);
	
		// Second stage
		int nz[size1];
	
		for (int i = 0; i < size1; i++)
			nz[i] = 0;
		
		// Count non-zeros on each line
		if (mm_is_symmetric(mm_code))
			for (int i = 0; i < nnz; i++)
			{
				if (I[i] == J[i])
					nz[I[i]]++;
				else
				{
					nz[I[i]]++;
					nz[J[i]]++;
				}
			}
		else
			for (int i = 0; i < nnz; i++)
				nz[I[i]]++;

		// Find out total number of non-zeros			
		int nnz2;
	
		if (mm_is_symmetric(mm_code))
		{
			int diagonals = 0;
		
			for (int i = 0; i < nnz; i++)
				if (I[i] == J[i])
					diagonals++;
				
			nnz2 = diagonals + 2 * (nnz - diagonals);
		}
		else
			nnz2 = nnz;

		// Fill in an almost-CSR data structure
		int *filled = new int[size1];
		int *indices = new int[size1];
		int *columns = new int[nnz2];
		double *values = new double[nnz2];
	
		indices[0] = 0;
		for (int i = 1; i < size1; i++)
			indices[i] = indices[i - 1] + nz[i - 1];
		
		for (int i = 0; i < size1; i++)
			filled[i] = 0;

		if (mm_is_symmetric(mm_code))
			for (int i = 0; i < nnz; i++)
				if (I[i] == J[i])
				{
					int index;
				
					index = indices[I[i]] + filled[I[i]];
					columns[index] = J[i];
					values[index] = V[i];
					filled[I[i]]++;
				}
				else
				{
					int index;
				
					index = indices[I[i]] + filled[I[i]];
					columns[index] = J[i];
					values[index] = V[i];
					filled[I[i]]++;
				
					index = indices[J[i]] + filled[J[i]];
					columns[index] = I[i];
					values[index] = V[i];
					filled[J[i]]++;
				}
		else
			for (int i = 0; i < nnz; i++)
				{
					int index;
				
					index = indices[I[i]] + filled[I[i]];
					columns[index] = J[i];
					values[index] = V[i];
					filled[I[i]]++;
				}
	
		// Create the matrix
		CompressedMatrix *m = new CompressedMatrix(size1, size2, nnz2);
	
		int k = 0;
	
		for (int i = 0; i < size1; i++)
			for (int j = 0; j < nz[i]; j++)
				(*m)(i, columns[indices[i] + j]) = values[k++];
		
		M = *m;
		
		delete[] I;
		delete[] J;
		delete[] V;

		delete[] filled;
		delete[] indices;
		delete[] columns;
		delete[] values;

		delete m;

		return true;
	}

	static bool WriteMatrixMarketMatrix(char *FileName, CompressedMatrix &M, bool Symmetric)
	{
		// Open MM file for writing
		FILE *f = fopen(FileName, "w");
	
		if (f == NULL)
		{
			printf("WriteMatrixMarketMatrix(): unable to open %s.\n", FileName);
			return false;
		}

		// Write MM file header
		MM_typecode mm_code;
	
		mm_initialize_typecode(&mm_code);	
	
		mm_set_matrix(&mm_code);
		mm_set_coordinate(&mm_code);
		mm_set_real(&mm_code);
	
		if (Symmetric)
			mm_set_symmetric(&mm_code);
		else
			mm_set_general(&mm_code);

		mm_write_banner(f, mm_code);
		
		// Find out the actual number of non-zeros in case of a symmetric matrix
		int nnz;
		
		if (Symmetric)
		{
			nnz = 0;
			
			CompressedMatrix::iterator1 a_iterator = M.begin1();

			for (unsigned int i = 0; i < M.size1(); i++)
			{
				#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
				for (CompressedMatrix::iterator2 row_iterator = a_iterator.begin(); row_iterator != a_iterator.end(); ++row_iterator) 
				#else
				for (CompressedMatrix::iterator2 row_iterator = begin(a_iterator, iterator1_tag()); row_iterator != end(a_iterator, iterator1_tag()); ++row_iterator)
				#endif
				{
					if (a_iterator.index1() >= row_iterator.index2())
						nnz++;
				}

			   a_iterator++;
			}
		}
		else
			nnz = M.nnz();
		
		// Write MM file sizes
		mm_write_mtx_crd_size(f, M.size1(), M.size2(), nnz);

		if (Symmetric)
		{
			CompressedMatrix::iterator1 a_iterator = M.begin1();

			for (unsigned int i = 0; i < M.size1(); i++)
			{
				#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
				for (CompressedMatrix::iterator2 row_iterator = a_iterator.begin(); row_iterator != a_iterator.end(); ++row_iterator) 
				#else
				for (CompressedMatrix::iterator2 row_iterator = begin(a_iterator, iterator1_tag()); row_iterator != end(a_iterator, iterator1_tag()); ++row_iterator)
				#endif
				{
					int I = a_iterator.index1(), J = row_iterator.index2();
				
					if (I >= J)
						if (fprintf(f, "%d %d %g\n", I + 1, J + 1, *row_iterator) < 0)
						{
							printf("WriteMatrixMarketMatrix(): unable to write data.\n");
							fclose(f);
							return false;
						}
				}

			   a_iterator++;
			}
		}
		else
		{
			CompressedMatrix::iterator1 a_iterator = M.begin1();

			for (unsigned int i = 0; i < M.size1(); i++)
			{
				#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
				for (CompressedMatrix::iterator2 row_iterator = a_iterator.begin(); row_iterator != a_iterator.end(); ++row_iterator) 
				#else
				for (CompressedMatrix::iterator2 row_iterator = begin(a_iterator, iterator1_tag()); row_iterator != end(a_iterator, iterator1_tag()); ++row_iterator)
				#endif
				{
					int I = a_iterator.index1(), J = row_iterator.index2();
				
					if (fprintf(f, "%d %d %g\n", I + 1, J + 1, *row_iterator) < 0)
					{
						printf("WriteMatrixMarketMatrix(): unable to write data.\n");
						fclose(f);
						return false;
					}
				}

			   a_iterator++;
			}
		}
	
		fclose(f);
	
		return true;
	}
	
	// Vector I/O routines
	
	static bool ReadMatrixMarketVector(char *FileName, Vector &V)
	{
		// Open MM file for reading
		FILE *f = fopen(FileName, "r");
	
		if (f == NULL)
		{
			printf("ReadMatrixMarketVector(): unable to open %s.\n", FileName);
			return false;
		}

		// Process MM file header
		MM_typecode mm_code;
	
		if (mm_read_banner(f, &mm_code) != 0)
		{
			printf("ReadMatrixMarketVector(): unable to read MatrixMarket banner.\n");
			fclose(f);
			return false;
		}
	
		if (!mm_is_valid(mm_code))
		{
			printf("ReadMatrixMarketVector(): invalid MatrixMarket banner.\n");
			fclose(f);
			return false;
		}
	
		// Check for supported types of MM file
		if (!((mm_is_real(mm_code) || mm_is_integer(mm_code)) && mm_is_array(mm_code)))
		{
			printf("ReadMatrixMarketVector(): invalid MatrixMarket type, \"%s\".\n",  mm_typecode_to_str(mm_code));
			fclose(f);
			return false;
		}
	
		// Read MM dimensions
		int size1, size2;
	
		if (mm_read_mtx_array_size(f, &size1, &size2) != 0)
		{
			printf("ReadMatrixMarketVector(): cannot read dimensions.\n");
			fclose(f);
			return false;
		}
		
		// Check MM dimensions
		if (size2 != 1)
		{
			printf("ReadMatrixMarketVector(): not a N x 1 array.\n");
			fclose(f);
			return false;
		}
	
		Vector *v = new Vector(size1);
		double T;
	
		// Read MM file
	
		for (int i = 0; i < size1; i++)
		{
			if (fscanf(f, "%lg", &T) != 1)
			{
				printf("ReadMatrixMarketVector(): invalid data.\n");
				fclose(f);

				return false;
			}
			
			(*v)(i) = T;		
		}
	
		fclose(f);
		
		V = *v;
		
		delete v;

		return true;
	}

	static bool WriteMatrixMarketVector(char *FileName, Vector &V)
	{
		// Open MM file for writing
		FILE *f = fopen(FileName, "w");
	
		if (f == NULL)
		{
			printf("WriteMatrixMarketVector(): unable to open %s.\n", FileName);
			return false;
		}

		// Write MM file header
		MM_typecode mm_code;
	
		mm_initialize_typecode(&mm_code);	
	
		mm_set_matrix(&mm_code);
		mm_set_array(&mm_code);
		mm_set_real(&mm_code);
	
		mm_write_banner(f, mm_code);
		
		// Write MM file sizes
		mm_write_mtx_array_size(f, V.size(), 1);

		for (unsigned int i = 0; i < V.size(); i++)
			if (fprintf(f, "%g\n", V(i)) < 0)
			{
				printf("WriteMatrixMarketVector(): unable to write data.\n");
				fclose(f);
				return false;
			}
	
		fclose(f);
	
		return true;
	}

} // namespace Kratos

#endif // KRATOS_MATRIX_MARKET_INTERFACE_H_INCLUDED  defined 
