//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2004/02/02 16:51:26 $
//   Revision:            $Revision: 1.3 $
//
//


#if !defined(KRATOS_CUTHILL_MCKEE_REORDERER_H_INCLUDED )
#define  KRATOS_CUTHILL_MCKEE_REORDERER_H_INCLUDED



// System includes 
#include <vector>

// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"


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
  
  /// 
  /** .
      two template parameter: 
      - TSparseSpaceType which specify type
        of the unknowns, coefficients, sparse matrix, vector of
	unknowns, right hand side vector and their respective operators.
      - TDenseMatrixType which specify type of the
        matrices used as temporary matrices or multi solve unknowns and
	right hand sides and their operators.  
  */
  template<class TSparseSpaceType, class TDenseSpaceType>
    class CuthillMcKeeReorderer : public Reorderer<TSparseSpaceType, TDenseSpaceType>
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Counted pointer of CuthillMcKeeReorderer
      typedef boost::shared_ptr<CuthillMcKeeReorderer> Pointer;

	  typedef Reorderer<TSparseSpaceType, TDenseSpaceType> BaseType;

      typedef typename TSparseSpaceType::MatrixType SparseMatrixType;
  
      typedef typename TSparseSpaceType::VectorType VectorType;
  
      typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

      typedef BaseType::IndexType IndexType;

	  typedef typename BaseType::IndexVectorType IndexVectorType;
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      CuthillMcKeeReorderer(){}

      /// Destructor.
      virtual ~CuthillMcKeeReorderer(){}
      

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{

      virtual void Initialize(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
	{
	  CalculateIndexPermutation(rA);
	}

      virtual void Reorder(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
	{
	}
      
      virtual void Reorder(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB)
	{
	}

 template<class TMatrixType>
 void WriteMatrixForGid(TMatrixType& rM, std::string MatrixFileName, std::string MatrixName)
{
	KRATOS_TRY
		
		unsigned int size_1 = rM.size1();
		unsigned int size_2 = rM.size2();



		//GiD_FILE matrix_file = GiD_fOpenPostMeshFile("matricesA.post.msh", GiD_PostAscii);
		std::ofstream matrix_file(MatrixFileName.c_str());

		//GiD_fBeginMesh(matrix_file, "Kratos Matrix",GiD_2D,GiD_Point,1);
		matrix_file << "MESH \"" << MatrixName << "\" dimension 2 ElemType Point Nnode 1" << std::endl;
        //GiD_fBeginCoordinates(matrix_file);
		matrix_file << "Coordinates" << std::endl;

		std::vector<double> results;
		int index = 1;
		double SolutionTag=0;
                
  
		// Creating a node for each nonzero of matrix
		for (TMatrixType::iterator1 i1 = rM.begin1(); i1 != rM.end1(); ++i1) 
			for (TMatrixType::iterator2 i2 = i1.begin(); i2 != i1.end(); ++i2)
			{
				//GiD_fWriteCoordinates(matrix_file, index++, i2.index1(), i2.index2(), 0);
				matrix_file << index++ << "  " << i2.index1() << "  " << size_1 - i2.index2() << "  0" << std::endl;
				results.push_back(*i2);
			}
                
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
 template<class TMatrixType>
 void WritePermutedMatrixForGid(TMatrixType& rM, std::string MatrixFileName, std::string MatrixName)
{
	KRATOS_TRY
		
		unsigned int size_1 = rM.size1();
		unsigned int size_2 = rM.size2();

		IndexVectorType& r_index_permutation = GetIndexPermutation();
		IndexVectorType invperm(size_1);

		//GiD_FILE matrix_file = GiD_fOpenPostMeshFile("matricesA.post.msh", GiD_PostAscii);
		std::ofstream matrix_file(MatrixFileName.c_str());

		//GiD_fBeginMesh(matrix_file, "Kratos Matrix",GiD_2D,GiD_Point,1);
		matrix_file << "MESH \"" << MatrixName << "\" dimension 2 ElemType Point Nnode 1" << std::endl;
        //GiD_fBeginCoordinates(matrix_file);
		matrix_file << "Coordinates" << std::endl;

		std::vector<double> results;
		int index = 1;
		double SolutionTag=0;
		
		for(unsigned int i = 0 ; i < size_1 ; i++) invperm[r_index_permutation[i]]=i;



		// Creating a node for each nonzero of matrix
		for (TMatrixType::iterator1 i1 = rM.begin1(); i1 != rM.end1(); ++i1) 
			for (TMatrixType::iterator2 i2 = i1.begin(); i2 != i1.end(); ++i2)
			{
				//GiD_fWriteCoordinates(matrix_file, index++, i2.index1(), i2.index2(), 0);
				//matrix_file << index++ << "  " << r_index_permutation[i2.index1()] << "  " << size_1 - r_index_permutation[i2.index2()] << "  0" << std::endl;
				matrix_file << index++ << "  " << invperm[i2.index1()] << "  " << size_1 - invperm[i2.index2()] << "  0" << std::endl;
				results.push_back(*i2);
			}
                
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
     virtual IndexVectorType& CalculateIndexPermutation(SparseMatrixType& rA, IndexType InitialIndex = IndexType())
	{

		typedef std::multimap<SizeType,IndexType, std::greater<SizeType> > level_set_type;

	  const unsigned int size = TSparseSpaceType::Size1(rA);

	  IndexVectorType& r_index_permutation = GetIndexPermutation();

	  r_index_permutation.resize(size);

	  std::vector<bool> is_marked(size, false);
	  level_set_type level_set;
	  IndexVectorType connectivity;
	  //IndexVectorType next_level_set;
	  level_set_type next_level_set;
	  

	  r_index_permutation[0] = InitialIndex;
	  is_marked[InitialIndex] = true;

	  //WriteMatrixForGid(rA, "matrixA.post.msh", "unordered_matrix");
	  
	  unsigned int next = 1 ;
	  
	  level_set.insert(level_set_type::value_type(TSparseSpaceType::GraphDegree(InitialIndex, rA), InitialIndex));
	  while(next < size)
	  {
		  //std::cout << "LevelSet" << level_set << std::endl;
		  

		  //while(level_set.size() != 0)
		  //{
			 // TSparseSpaceType::GraphNeighbors(level_set.begin()->second, rA, connectivity);
			 // 
			 // //std::cout << "Connectivity" << connectivity << std::endl;

			 // for(IndexVectorType::iterator j = connectivity.begin() ; j != connectivity.end() ; ++j)
				//  if(is_marked[*j] == false)
				//  {
				//	  r_index_permutation[next++] = *j;
				//	  is_marked[*j] = true;
				//		
				//	  SizeType degree = TSparseSpaceType::GraphDegree(*j, rA);
				//	  if(degree != 0)
				//		level_set.insert(level_set_type::value_type(degree,*j));
				//  }

			 // level_set.erase(level_set.begin());
		  //}
		  for(level_set_type::iterator i = level_set.begin() ; i != level_set.end() ; ++i)
		  {
			  TSparseSpaceType::GraphNeighbors(i->second, rA, connectivity);
			  
			  //std::cout << "Connectivity" << connectivity << std::endl;

			  for(IndexVectorType::iterator j = connectivity.begin() ; j != connectivity.end() ; ++j)
				  if(is_marked[*j] == false)
				  {
					  r_index_permutation[next++] = *j;
					  is_marked[*j] = true;

					  SizeType degree = TSparseSpaceType::GraphDegree(*j, rA);
					  if(degree != 0)
						next_level_set.insert(level_set_type::value_type(degree,*j));
				  }
				  
		  }

		  //std::cout << "NextLevelSet" << next_level_set << std::endl;

		  level_set.swap(next_level_set);
		  next_level_set.clear();
		   	      if((level_set.size() == 0) && (next < size)) // No connected graph
		   		{
		   		  SizeType k = 0;
		   		  while(is_marked[k])
		   		    k++;
		   		  level_set.insert(level_set_type::value_type(TSparseSpaceType::GraphDegree(k, rA), k));
		   		  r_index_permutation[next++] = k;
		   		  is_marked[k] = true;
				  //next_level_set.push_back(*j);
		  
		   		}
	  }
	  // for testing reverse cuthill mckee
	  //IndexVectorType temp(size);
	  //for(int i = 0 ; i < size ; i++) temp[size - i - 1] = r_index_permutation[i];
	  //r_index_permutation = temp;


	  //std::cout << "r_index_permutation : ";
	  //for(int i = 0 ; i < r_index_permutation.size() ; i++)
		 // std::cout << r_index_permutation[i] << ",";
	  //std::cout << std::endl;

  	  WritePermutedMatrixForGid(rA, "matrixApermuted.post.msh", "unordered_matrix");

	  return r_index_permutation;
	}

      ///@}
      ///@name Access
      ///@{ 
      
      
      ///@}
      ///@name Inquiry
      ///@{
      
      
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

        
        
      ///@} 
      ///@name Private Operators
      ///@{ 
        
        
      ///@} 
      ///@name Private Operations
      ///@{ 
        
        
      ///@} 
      ///@name Private  Access 
      ///@{ 
        
        
      ///@}    
      ///@name Private Inquiry 
      ///@{ 
        
        
      ///@}    
      ///@name Un accessible methods 
      ///@{ 
      
      /// Assignment operator.
      CuthillMcKeeReorderer& operator=(const CuthillMcKeeReorderer& Other);

      /// Copy constructor.
        CuthillMcKeeReorderer(const CuthillMcKeeReorderer& Other);

        
      ///@}    
        
    }; // Class CuthillMcKeeReorderer 

  ///@} 
  
  ///@name Type Definitions       
  ///@{


  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
//   inline std::istream& operator >> (std::istream& IStream, 
// 				    CuthillMcKeeReorderer& rThis);
  
  /// output stream function
//   inline std::ostream& operator << (std::ostream& OStream, 
// 				    const CuthillMcKeeReorderer& rThis);
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_CUTHILL_MCKEE_REORDERER_H_INCLUDED  defined 


