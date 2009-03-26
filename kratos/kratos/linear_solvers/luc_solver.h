/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
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
//   Project Name:        Kratos       
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-03-06 10:30:33 $
//   Revision:            $Revision: 1.3 $
//
//


#if !defined(KRATOS_LUC_SOLVER_H_INCLUDED )
#define  KRATOS_LUC_SOLVER_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <fstream>


// External includes 
#include <boost/timer.hpp>


// Project includes
#include "includes/define.h"
#include "linear_solvers/iterative_solver.h"


namespace Kratos
{

  ///@name Kratos Globals
  ///@{ 
  
  ///@} 
  ///@name Type Definitions
  ///@{ 
  
  ///@} 
  ///@name  Enum's
  ///@{
      
  ///@}
  ///@name  Functions 
  ///@{
      
  ///@}
  ///@name Kratos Classes
  ///@{
  
  /// Short class definition.
  /** Detail class definition.
  */
  template<class TSparseSpaceType, class TDenseSpaceType, 
    class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
    class LUCSolver : public DirectSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>
    {
      public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of LUCSolver
      KRATOS_CLASS_POINTER_DEFINITION(LUCSolver);

      typedef DirectSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType; 
  
      typedef typename TSparseSpaceType::MatrixType SparseMatrixType;
  
      typedef typename TSparseSpaceType::VectorType VectorType;
  
      typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

	  typedef std::size_t IndexType;
  
	  typedef std::vector<IndexType> IndicesVectorType;

	  ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      LUCSolver(){}

      /// Constructor with specific reorderer.
	  LUCSolver(typename TReordererType::Pointer pNewReorderer) : BaseType(pNewReorderer){}

	  /// Copy constructor.
      LUCSolver(const LUCSolver& Other) : BaseType(Other) {}


      /// Destructor.
      virtual ~LUCSolver(){}
      

      ///@}
      ///@name Operators 
      ///@{
      
      /// Assignment operator.
      LUCSolver& operator=(const LUCSolver& Other)
      {
        BaseType::operator=(Other);
		return *this;
      }
      
      ///@}
      ///@name Operations
      ///@{
      
      /** Normal solve method.
	  Solves the linear system Ax=b and puts the result on SystemVector& rX. 
	  @param rA. System matrix
	  @param rX. Solution vector. it's also the initial 
	  guess for iterative linear solvers.
 	  @param rB. Right hand side vector.
      */
	  bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
	  {
		  GetReorderer()->Initialize(rA, rX, rB);
		  std::cout << "	Decomposing Matrix A..." << std::endl;
		  boost::timer decomposing_timer;
		  if(LUCDecompose(rA) != 0)
			  return false;
		  std::cout << "	Decomposing time : " << decomposing_timer.elapsed() << std::endl;

		  std::cout << "	Solving L U X = B..." << std::endl;
		  boost::timer lu_solve_timer;
		  if(LUSolve(rA, rX, rB) == false)
			  return false;
		  std::cout << "	Solving L U X = B time : " << lu_solve_timer.elapsed() << std::endl;

		  return true;
	  }

      
      /** Multi solve method for solving a set of linear systems with same coefficient matrix.
	  Solves the linear system Ax=b and puts the result on SystemVector& rX. 
	  @param rA. System matrix
	  @param rX. Solution vector. it's also the initial 
	  guess for iterative linear solvers.
 	  @param rB. Right hand side vector.
      */
      bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB)
	{

	  return false;
	}
      
      ///@}
      ///@name Access
      ///@{ 
      
      
      ///@}
      ///@name Inquiry
      ///@{
      
      
      ///@}      
      ///@name Input and output
      ///@{

      /// Turn back information as a string.
      virtual std::string Info() const
	{
	  std::stringstream buffer;
	  buffer << "LUC linear solver";
	  return  buffer.str();
	}
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
	{
	  rOStream << Info();
	}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const
	{
	  BaseType::PrintData(rOStream);
	}
      
            
      ///@}      
      ///@name Friends
      ///@{
      
            
      ///@}
      
    protected:
      ///@name Protected static Member Variables 
      ///@{ 
        
        
      ///@} 
      ///@name Protected member Variables 
      ///@{ 
        
        
      ///@} 
      ///@name Protected Operators
      ///@{ 
        
        
      ///@} 
      ///@name Protected Operations
      ///@{ 
        
        
      ///@} 
      ///@name Protected  Access 
      ///@{ 
        
        
      ///@}      
      ///@name Protected Inquiry 
      ///@{ 
        
        
      ///@}    
      ///@name Protected LifeCycle 
      ///@{ 
      
            
      ///@}
      
    private:
      ///@name Static Member Variables 
      ///@{ 
        
        
      ///@} 
      ///@name Member Variables 
      ///@{ 
        
	  IndicesVectorType mIL, mJL, mIU, mJU;
	  std::vector<double> mL, mU;
        
      ///@} 
      ///@name Private Operators
      ///@{ 
        
        
      ///@} 
      ///@name Private Operations
      ///@{ 

	  void ClearLU()
	  {
			mIL.clear();
			mJL.clear();
			mIU.clear();
			mJU.clear();
			mL.clear();
			mU.clear();
	  }
        
 void WriteMatrixForGid(std::string MatrixFileName, std::string MatrixName)
{
	KRATOS_TRY
		



		//GiD_FILE matrix_file = GiD_fOpenPostMeshFile("matricesA.post.msh", GiD_PostAscii);
		std::ofstream matrix_file(MatrixFileName.c_str());

		//GiD_fBeginMesh(matrix_file, "Kratos Matrix",GiD_2D,GiD_Point,1);
		matrix_file << "MESH \"" << MatrixName << "\" dimension 2 ElemType Point Nnode 1" << std::endl;
        //GiD_fBeginCoordinates(matrix_file);
		matrix_file << "Coordinates" << std::endl;

		std::vector<double> results;
		int index = 1;
		double SolutionTag=0;
  
		int size_1 = mIU.size() - 1;

		// Creating a node for each nonzero of matrix U
		for (unsigned int i = 0 ; i < size_1 ; i++) 
			for (unsigned int j = mIU[i] ; j < mIU[i+1] ; j++)
			{
				//GiD_fWriteCoordinates(matrix_file, index++, i2.index1(), i2.index2(), 0);
				matrix_file << index++ << "  " << i << "  " << size_1 - mJU[j] << "  0" << std::endl;
				results.push_back(mU[j]);
			}
		// Creating a node for each nonzero of matrix L
		//for (unsigned int j = 0 ; j < mJL.size() - 1 ; j++) 
		//	for (unsigned int i = mJL[j] ; j < mJL[j+1] ; i++)
		//	{
		//		//GiD_fWriteCoordinates(matrix_file, index++, i2.index1(), i2.index2(), 0);
		//		matrix_file << index++ << "  " << mIL[i] << "  " << size_1 - j << "  0" << std::endl;
		//		results.push_back(mL[j]);
		//	}
  //              
		//GiD_fEndCoordinates(matrix_file);
		matrix_file << "End Coordinates" << std::endl;

		// Creating an element for each nonzero of matrix
		int nodes_id[1];
        //GiD_fBeginElements(matrix_file);
		matrix_file << "Elements" << std::endl;
        for(  int i = 1 ; i < index ; i++)
		{
                    nodes_id[0] = i;
                    //GiD_fWriteElement(matrix_file, i,nodes_id);
					matrix_file << i << "  " << i << std::endl;
		}
                
		//GiD_fEndElements(matrix_file);
		matrix_file << "End Elements" << std::endl;
        //GiD_fEndMesh(matrix_file);
		//GiD_fClosePostMeshFile(matrix_file);


            
		//GiD_OpenPostResultFile("matrices.post.bin", GiD_PostBinary);
		//GiD_BeginResult( "Matrix", "Kratos", 
  //                       SolutionTag, GiD_Scalar, 
  //                       GiD_OnNodes, NULL, NULL, 0, NULL );

  //   
		//for( int i = 1 ; i < index ; i++)
		//{
		//	GiD_WriteScalar( i, results[i-1]);
		//}
  //              
		//GiD_EndResult();
		//GiD_ClosePostResultFile();
	
	KRATOS_CATCH("")

}


	  /// Decomposes the Matrix A to L and U and store it in mL and mU

	  std::size_t LUCDecompose(SparseMatrixType& rA)
	  {
		  typedef typename TSparseSpaceType::DataType DataType;

		  const SizeType size = TSparseSpaceType::Size1(rA);

		  VectorType z(size); // working vector for row k
		  VectorType w(size); // working vector for column k
		  VectorType temp(size);

		  IndicesVectorType l_first(size, 0); // The pointer to the first element in each column of L with row index >= k (L is stored by columns)
		  IndicesVectorType u_first(size, 0); // The pointer to the first element in each row of U with column index >= k (U is stored by rows)
		  IndicesVectorType l_list(size, size+1); // The link list of elements in row k the L
		  IndicesVectorType u_list(size, size+1); // The link list of elements in column k the U
		  std::vector<typename SparseMatrixType::iterator2> a_first; // The pointer to the first element in each row of rA with column index >= k (U is stored by rows)
		  typename SparseMatrixType::iterator1 a_iterator = rA.begin1();

		  for(SizeType i = 0 ; i < size ; i++)
			  a_first.push_back((a_iterator++).begin()); // initializing the row_pointers of A

		  ClearLU();

		  mJL.push_back(0);
		  mIU.push_back(0);

		  SizeType output_index = 0;
		  for(SizeType k = 0 ; k < size ; k++)
		  {
			  //initializing z: z_1:k-1 = 0, z_k:n = rA_k,kn
			  TSparseSpaceType::GetRow(k, rA, z); // NOTE: I don't need to set the first part to zero while I'm not using it!

			  for(SizeType i = l_list[k] ; i < size ; i = l_list[i]) // loop over nonzeros in row k of L
				  for(SizeType j = u_first[i] ; j < mIU[i+1] ; j++) // for all nonzeros in row i of U
					  z[mJU[j]] -= mL[l_first[i]] * mU[j]; // z_j = Z_j - L_ki * U_ij

			  //initializing w: w_1:k = 0, w_k:n = rA_k+1:n,k
			  for(SizeType i = k ; i < size ; i++)
			  {
				  if(a_first[i].index2() < k)
					  (a_first[i])++;

				  if(a_first[i].index2() == k)
					  w[i] = *(a_first[i]);
				  else
					  w[i] = 0;
			  }


			  for(SizeType i = u_list[k] ; i < size ; i = u_list[i]) // loop over nonzeros in column k of U
				  for(SizeType j = l_first[i] ; j < mJL[i+1] ; j++) // for all nonzeros in column i of L
					  w[mIL[j]] -= mU[u_first[i]] * mL[j]; // w_j = w_j - L_ki * U_ij


			  // adding nonzeros of z to the U
			  for(SizeType i = k ; i < size ; i++)
			  {
				  if(z[i] != 0.00)
				  {
					  mJU.push_back(i);
					  mU.push_back(z[i]);
				  }
			  }

			  mIU.push_back(mU.size());

			  double u_kk = z[k];

			  if(u_kk == 0.00)
				  //{
				  //	KRATOS_WATCH(k);
				  //	for(SizeType k1 = k + 1 ; k1 < size ; k1++)
				  //		if(rA(r_index_permutation[k1], r_index_permutation[k1]) != 0.00)
				  //		{
				  //			KRATOS_WATCH(k1);
				  //			SizeType temp_index = r_index_permutation[k1];
				  //			r_index_permutation[k1] = r_index_permutation[k];
				  //			r_index_permutation[k] = temp_index;
				  //			break;
				  //		}
				  //		k--;
				  //		continue;
				  //}
				  KRATOS_ERROR(std::runtime_error, "Zero pivot found in row ",k);


			  // adding nonzeros of w to the L
			  for(SizeType i = k + 1 ; i < size ; i++)
			  {
				  if(w[i] != 0.00)
				  {
					  mIL.push_back(i);
					  mL.push_back(w[i]/u_kk);
				  }
			  }	  
				
			  mJL.push_back(mL.size());

			  // updating the l_first, l_list for added column
			  l_first[k] = mJL[k];

			  if(l_first[k]  < mIL.size())
				IndexPushBack(l_list, k, mIL[l_first[k]]);

			  // updating l_first and l_list in this columns
			  SizeType j = k;
			  if(k + 1 < size) // not for the last column!
			  {
				for(SizeType i = l_list[k] ; i < size ; i = l_list[i]) // loop over nonzeros in row k of L
				  {
					  l_list[j] = size + 1; // reseting the list
					  l_first[i]++;
					  if(l_first[i] < mJL[i+1])
					  {	
							IndexPushBack(l_list, i, mIL[l_first[i]]);
					  }
					  
					  j = i;
				  }
			  }

			  // updating the u_first, u_list for added row
				  u_first[k] = mIU[k];

			  if(u_first[k]  < mJU.size())
				IndexPushBack(u_list, k, mJU[u_first[k]]);

			  // updating u_first and u_list in this columns
			  j = k;
			  if(k + 1 < size) // not for the last row!
			  {
				for(SizeType i = u_list[k] ; i < size ; i = u_list[i]) // loop over nonzeros in column k of U
				  {
					  u_list[j] = size + 1; // reseting the list
					  u_first[i]++;
					  if(u_first[i] < mIU[i+1])
					  {	
							IndexPushBack(u_list, i, mJU[u_first[i]]);
					  }
					  
					  j = i;
				  }
			  }

				  if(output_index++ >= (0.1 * size))
				  {
					  std::cout << "            " << int(double(k) / double(size) * 100.00) << " % : L nonzeros = " << mL.size() << " and U nonzeros = " << mU.size() << std::endl;
					  output_index = 0;
				  }

				  //std::cout << "Finishing loop #" << k << " of " << size << std::endl;
		  }

		  std::cout << "            " << 100 << " % : L nonzeros = " << mL.size() << " and U nonzeros = " << mU.size() << std::endl;

		  // //Writing U
		  //for(SizeType i = 0 ; i < size ; i++)
		  //{
		  //	int k = mIU[i];
		  //	for(SizeType j = 0 ; j < size ; j++)
		  //	{
		  //		double value = 0.00;
		  //		if(mJU[k] == j)
		  //		{
		  //			value = mU[k];
		  //			k++;
		  //		}
		  //		std::cout << value << ", ";
		  //	}
		  //	std::cout << std::endl;
		  //}

		  //// Writing L
		  //for(SizeType i = 0 ; i < size ; i++)
		  //{
		  //	for(SizeType j = 0 ; j < i ; j++)
		  //	{
		  //		double value = 0.00;
		  //		for(SizeType k = mJL[j] ; k < mJL[j+1] ; k++)
		  //			if(mIL[k] == i)
		  //				value = mL[k];
		  //		std::cout << value << ", ";
		  //	}
		  //	std::cout << std::endl;
		  //}


		  //WriteMatrixForGid("matrixLU.post.msh", "lu_matrix");

		  return 0;      
	  }

	  void IndexPushBack(IndicesVectorType& ThisLinkList, IndexType ThisIndex, IndexType ThisRow)
	  {
		  SizeType i = ThisRow;
		   // Loop to the end of the link list of row ThisRow
		  for(; ThisLinkList[i] < ThisLinkList.size() ; i = ThisLinkList[i]);

		  ThisLinkList[i] = ThisIndex;
	  }
			

 
	  std::size_t LUCPermuteAndDecompose(SparseMatrixType& rA)
	  {
		    typename TReordererType::IndexVectorType& r_index_permutation = GetReorderer()->GetIndexPermutation();

			typedef typename TSparseSpaceType::DataType DataType;

			const SizeType size = TSparseSpaceType::Size1(rA);

			DataType* z = new DataType[size]; // working vector for row k
			VectorType w(size); // working vector for column k
			VectorType temp(size);

			IndicesVectorType l_first(size, 0); // The pointer to the first element in each column of L with row index >= k (L is stored by columns)
			IndicesVectorType u_first(size, 0); // The pointer to the first element in each row of U with column index >= k (U is stored by rows)

			std::vector<IndicesVectorType> l_rows(size,IndicesVectorType());
			std::vector<IndicesVectorType> u_columns(size,IndicesVectorType());

			mIL.clear();
			mJL.clear();
			mIU.clear();
			mJU.clear();
			mL.clear();
			mU.clear();

			mJL.push_back(0);
			mIU.push_back(0);

			// find inverse permutation
			typename TReordererType::IndexVectorType inverse_permutation(size);
			for (SizeType i=0; i<size; i++) 
				inverse_permutation[r_index_permutation[i]]=i;

			SizeType output_index = 0;
			for(SizeType k = 0 ; k < size ; k++)
			{
				//KRATOS_WATCH(k);
				//initializing z: z_1:k-1 = 0, z_k:n = rA_k,kn
				// this has to be changed to be without temporary vector
				TSparseSpaceType::GetRow(r_index_permutation[k], rA, temp);

				for(SizeType i = 0 ; i < size ; i++)
					z[i] = temp[r_index_permutation[i]];
				for(SizeType i = 0 ; i < k ; i++)
					z[i] = 0;


				for(SizeType h = 0 ; h < l_rows[k].size() ; h++) // iterating over nonzeros of row k of L
				{
					SizeType i = l_rows[k][h];
					while((u_first[i] < mJU.size()) && (mJU[u_first[i]] < k)) // updating u_first
						u_first[i]++;
					while((l_first[i] < mIL.size()) && (mIL[l_first[i]] < k)) // updating l_first 
						l_first[i]++;

					for(SizeType j = u_first[i] ; j < mIU[i+1] ; j++) // for all nonzeros in row i of U
						z[mJU[j]] -= mL[l_first[i]] * mU[j]; // z_j = Z_j - L_ki * U_ij
				}

				//initializing w: w_1:k = 0, w_k:n = rA_k+1:n,k
				// this has to be changed to be without temporary vector
				TSparseSpaceType::GetColumn(r_index_permutation[k], rA, temp);
				for(SizeType i = 0 ; i < size ; i++)
					w[i] = temp[r_index_permutation[i]];
				for(SizeType i = 0 ; i < k ; i++)
					w[i] = 0;

				for(SizeType h = 0 ; h < u_columns[k].size() ; h++) // iterating over nonzeros of column k of U
				{
					SizeType i = u_columns[k][h];
					while((u_first[i] < mJU.size()) && (mJU[u_first[i]] < k)) // updating u_first
						u_first[i]++;
					while((l_first[i] < mIL.size()) && (mIL[l_first[i]] < k)) // updating l_first 
						l_first[i]++;

					for(SizeType j = l_first[i] ; j < mJL[i+1] ; j++) // for all nonzeros in column i of L
						w[mIL[j]] -= mU[u_first[i]] * mL[j]; // w_j = w_j - L_ki * U_ij
				}

				// adding nonzeros of z to the U
				for(SizeType i = k ; i < size ; i++)
					if(z[i] != 0.00)
					{
						u_columns[i].push_back(k);
						mJU.push_back(i);
						mU.push_back(z[i]);
					}
				mIU.push_back(mU.size());

				double u_kk = z[k];

				if(u_kk == 0.00)
					KRATOS_ERROR(std::runtime_error, "Zero pivot found in row ",k);


				// adding nonzeros of w to the L
				for(SizeType i = k + 1 ; i < size ; i++)
					if(w[i] != 0.00)
					{
						l_rows[i].push_back(k);
						mIL.push_back(i);
						mL.push_back(w[i]/u_kk);
					}
				mJL.push_back(mL.size());

				l_first[k] = mJL[k];
				u_first[k] = mIU[k];

				if(output_index++ >= (0.1 * size))
				{
					std::cout << "            " << int(double(k) / double(size) * 100.00) << " % : L nonzeros = " << mL.size() << " and U nonzeros = " << mU.size() << std::endl;
					output_index = 0;
				}
			}
					
			std::cout << "            " << 100 << " % : L nonzeros = " << mL.size() << " and U nonzeros = " << mU.size() << std::endl;
			
			
			delete [] z;

	  WriteMatrixForGid("matrixLU.post.msh", "lu_matrix");

			return 0;      
	  }
        
	template<class T>
	  void WriteVector(std::string Name, std::vector<T>& Data)
	  {
		  std::cout << Name << " : ";
			for(int i = 0 ; i < Data.size() ; i++)
				std::cout << Data[i] << ", ";
			std::cout << std::endl;
	  }

	  bool LUSolve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
	  {
		  VectorType y = ZeroVector(rX.size());
		  std::cout << "		Solving L Y = B..." << std::endl;
		  boost::timer u_solve_timer;
		  if(LSolve(rA, y, rB) == false)
			  return false;
		  std::cout << "		Solving L Y = B time : " << u_solve_timer.elapsed() << std::endl;

		  //KRATOS_WATCH(y);

		  std::cout << "		Solving U X = Y..." << std::endl;
		  boost::timer l_solve_timer;
		  if(USolve(rA, rX, y) == false)
			  return false;
		  std::cout << "		Solving U X = Y time : " << l_solve_timer.elapsed() << std::endl;

		  return true;
	  }

	  bool LSolve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
	  {
		   typename TReordererType::IndexVectorType& r_index_permutation = GetReorderer()->GetIndexPermutation();
		  SizeType size = rX.size();
		
		  for(SizeType i = 0 ; i < size ; i++)
		  {
				  rX[i] = rB[r_index_permutation[i]];
		  }

		  for(SizeType i = 0 ; i < size ; i++)
		  {
			  for(SizeType j = mJL[i] ; j < mJL[i+1] ; j++)
				  rX[mIL[j]] -= mL[j] * rX[i];
		  }

		  return true;
	  }
 
	  bool USolve(SparseMatrixType& rA, VectorType& rX, VectorType& rY)
	  {
		   typename TReordererType::IndexVectorType& r_index_permutation = GetReorderer()->GetIndexPermutation();

		  SizeType size = rX.size();
		  VectorType x(size);

		  for(SizeType i = 0 ; i < size ; i++) x[i] = 0.00;

		  for(int i = size - 1 ; i >= 0 ; i--)
		  {
			  double temp = rY[i];
			  for(SizeType j = mIU[i] + 1 ; j < mIU[i+1] ; j++)
				  temp -= mU[j] * x[mJU[j]];
			  x[i] = temp / mU[mIU[i]];
		  }

		  for(SizeType i = 0 ; i < size ; i++) rX[r_index_permutation[i]] = x[i];
		  return true;
	  }
 
      ///@} 
      ///@name Private  Access 
      ///@{ 
        
        
      ///@}    
      ///@name Private Inquiry 
      ///@{ 
        
        
      ///@}    
      ///@name Un accessible methods 
      ///@{ 
      
        
      ///@}    
        
    }; // Class LUCSolver 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  template<class TSparseSpaceType, class TDenseSpaceType, 
    class TPreconditionerType, 
    class TReordererType>
  inline std::istream& operator >> (std::istream& IStream, 
				      LUCSolver<TSparseSpaceType, TDenseSpaceType, 
				      TReordererType>& rThis)
    {
		return IStream;
    }

  /// output stream function
  template<class TSparseSpaceType, class TDenseSpaceType, 
    class TPreconditionerType, 
    class TReordererType>
  inline std::ostream& operator << (std::ostream& OStream, 
				    const LUCSolver<TSparseSpaceType, TDenseSpaceType, 
				      TReordererType>& rThis)
    {
      rThis.PrintInfo(OStream);
      OStream << std::endl;
      rThis.PrintData(OStream);

      return OStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_LUC_SOLVER_H_INCLUDED  defined 


