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
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2008-06-20 17:01:18 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_PETSC_SPACE_H_INCLUDED )
#define  KRATOS_PETSC_SPACE_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <cstddef>


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"


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
  template<class TMatrixType, class TVectorType>
  class PetscSpace
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of PetscSpace
      KRATOS_CLASS_POINTER_DEFINITION(PetscSpace);
  
      typedef double DataType;

      typedef TMatrixType MatrixType;

      typedef TVectorType VectorType;

      typedef std::size_t IndexType;
  
      typedef std::size_t SizeType;
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      PetscSpace(){}

      /// Destructor.
      virtual ~PetscSpace(){}
      

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{

      /// return size of vector rV
      static IndexType Size(VectorType const& rV)
      {
	PetscInt size;
	PetscInt ierr = VecGetSize(rV, &size); CHKERRQ(ierr);
	return size;
      } 
      
      /// return number of rows of rM 
      static IndexType Size1(MatrixType const& rM)
      {
	PetscInt size1;
	PetscInt size2;
	PetscInt ierr = MatGetSize(rM, &size1, &size2); CHKERRQ(ierr);
	return size1;
      }

      /// return number of columns of rM
      static IndexType Size2(MatrixType const& rM)
      {
	PetscInt size1;
	PetscInt size2;
	PetscInt ierr = MatGetSize(rM, &size1, &size2); CHKERRQ(ierr);
	return size2;
      }

      /// rXi = rMij
      static void GetColumn(unsigned int j, MatrixType& rM, VectorType& rX)
      {
	MatGetColumnVector(rM,rX,static_cast<PetscInt>(j));
      } 


///////////////////////////////// TODO: Take a close look to this method!!!!!!!!!!!!!!!!!!!!!!!!!
      /// rMij = rXi
      //      static void SetColumn(unsigned int j, MatrixType& rM, VectorType& rX){rX = row(rM, j);}  

      /// rY = rX
      static void Copy(MatrixType const& rX, MatrixType& rY)
      {
	MatDuplicate(rX,MAT_COPY_VALUES,rY);
      } 

      /// rY = rX
      static void Copy(VectorType const& rX, VectorType& rY)
      {
	VecCopy(rX, rY);
      } 

      /// rX * rY
      static double Dot(VectorType& rX, VectorType& rY)
      {
	PetscScalar value;
	PetscInt ierr = VecDot(rX,rY,&value); CHKERRQ(ierr);
	return value;
      } 

      /// ||rX||2
      static double TwoNorm(VectorType const& rX) 
      {
	PetscScalar value;
	PetscInt ierr = VecNorm(rX,NORM_2, &value); CHKERRQ(ierr);
	return value;
      } 


      static void Mult(MatrixType& rA, VectorType& rX, VectorType& rY)
	{
	  MatMult(rA, rX, rY); //CHKERRQ(ierr);
	}

//       static void TransposeMult(MatrixType& rA, VectorType& rX, VectorType& rY)
// 	{
// 	  axpy_prod(rX, rA, rY, true);
// 	} // rY = rAT * rX
 
//         static inline SizeType GraphDegree( IndexType i, TMatrixType& A)
//         {
//             typename MatrixType::iterator1 a_iterator = A.begin1();
//             std::advance(a_iterator,i);
//             #ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
//             return( std::distance( a_iterator.begin(), a_iterator.end() ) );
//             #else
//             return( std::distance( begin(a_iterator, boost::numeric::ublas::iterator1_tag()),
//                           end(a_iterator, boost::numeric::ublas::iterator1_tag()) ) );
//             #endif
//         }
        
//         static inline void GraphNeighbors( IndexType i, TMatrixType& A, std::vector<IndexType>& neighbors)
//         {
//             neighbors.clear();
//             typename MatrixType::iterator1 a_iterator = A.begin1();
//             std::advance(a_iterator,i);
//             #ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
//             for (typename MatrixType::iterator2 row_iterator = a_iterator.begin() ;
//             row_iterator != a_iterator.end() ; ++row_iterator) 
//             {
//             #else
//             for ( typename MatrixType::iterator2 row_iterator = begin(a_iterator,
//                   boost::numeric::ublas::iterator1_tag()); 
//             row_iterator != end(a_iterator, 
//                                 boost::numeric::ublas::iterator1_tag()); ++row_iterator )
//             {
//             #endif
//                 neighbors.push_back( row_iterator.index2() );
//             }
//         }


	//********************************************************************
	//checks if a multiplication is needed and tries to do otherwise
	static void InplaceMult(VectorType& rX, const double A)
	{
	  if( A != 1.00) 
		  VecScale (rX, A);
	}

	//********************************************************************
	//checks if a multiplication is needed and tries to do otherwise
	//ATTENTION it is assumed no aliasing between rX and rY
	// X = A*y;
	static void Assign(VectorType& rX, const double A, const VectorType& rY)
	{
	  
	  VecCopy(rY, rX);
	  if( A != 1.00) 
		  VecScale (rX, A);
	}

	//********************************************************************
	//checks if a multiplication is needed and tries to do otherwise
	//ATTENTION it is assumed no aliasing between rX and rY
	// X += A*y;
	static void UnaliasedAdd(VectorType& rX, const double A, const VectorType& rY)
	{
	  VecAXPY(rX,A,rY);
	}

 	//********************************************************************
     static void ScaleAndAdd(const double A, const VectorType& rX, const double B, const VectorType& rY, VectorType& rZ)  // rZ = (A * rX) + (B * rY)
	{
		Assign(rZ,A,rX); //rZ = A*rX
		UnaliasedAdd(rZ,B,rY); //rZ += B*rY

	} 

      static void ScaleAndAdd(const double A,const  VectorType& rX, const double B, VectorType& rY) // rY = (A * rX) + (B * rY) 
	{
		InplaceMult(rY,B);
		UnaliasedAdd(rY,A,rX);
	} 

      
      /// rA[i] * rX
//       static double RowDot(unsigned int i, MatrixType& rA, VectorType& rX)
// 	{
// 	  return inner_prod(row(rA, i), rX);
// 	} 

      
      /// rX = A
      static void Set(VectorType& rX, DataType A)
      {
	VecSet(rX, A);
      } 

      static void Resize(MatrixType& rA, SizeType m, SizeType n)
      {
	MatSetSizes(rA, PETSC_DECIDE, PETSC_DECIDE, m, n);
      }

      static void Resize(VectorType& rX, SizeType n) 
      {
	VecSetSizes(rX, PETSC_DECIDE, n);
      }
      
      static void Clear(MatrixType& rA)
      {
	Resize(rA,0,0);
      }

      static void Clear(VectorType& rX) 
      {
	Resize(rX,0);
      }
      
       template<class TOtherMatrixType>
	 inline static void ClearData(TOtherMatrixType& rA){SetToZero(rA);}

//       inline static void ClearData(compressed_matrix<TDataType>& rA)
//       {
// 	rA.clear();
// //    	rA.value_data() = unbounded_array<TDataType>();
//       	//if(rA.non_zeros() != 0) rA.value_data() = unbounded_array<TDataType>();
//       }

//       inline static void ClearData(VectorType& rX) {rX = VectorType();}
      
//       template<class TOtherMatrixType>
//       inline static void ResizeData(TOtherMatrixType& rA, SizeType m){rA.resize(m,false);
//       std::fill(rA.begin(), rA.end(), TDataType());}

//       inline static void ResizeData(compressed_matrix<TDataType>& rA, SizeType m) {rA.value_data().resize(m);
//       std::fill(rA.value_data().begin(), rA.value_data().end(), TDataType());}

//       inline static void ResizeData(VectorType& rX, SizeType m) {rX.resize(m);
//       std::fill(rX.begin(), rX.end(), TDataType());}
      
      inline static void SetToZero(MatrixType& rA)
      {
	MatZeroEntries(rA);
      }


      inline static void SetToZero(VectorType& rX) 
      {
	VecZeroEntries(rX);
      }

      template<class TOtherMatrixType, class TEquationIdVectorType>
      inline static void AssembleLHS(
	  MatrixType& A,
	   TOtherMatrixType& LHS_Contribution,
	  TEquationIdVectorType& EquationId
	  )
	  {
	    unsigned int system_size = Size1(A);
	    unsigned int local_size = LHS_Contribution.size();
	    PetscInt ierr;
	      
	      for (unsigned int i_local=0; i_local<local_size; i_local++)
	      {
		  unsigned int i_global=EquationId[i_local];
		  if ( i_global < system_size )
		  {
		      for (unsigned int j_local=0; j_local<local_size; j_local++)
		      {
			  unsigned int j_global=EquationId[j_local];
			  if ( j_global < system_size )
			      MatSetValue(A, i_global, j_global, LHS_Contribution(i_local,j_local), ADD_VALUES);
		      }
		  }
	      }
	      ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	      ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
	  }

/// TODO: creating the the calculating reaction version
      template<class TOtherVectorType, class TEquationIdVectorType>
	void AssembleRHS(
			VectorType& b,
			TOtherVectorType& RHS_Contribution,
			Element::EquationIdVectorType& EquationId
			)
		{
		  unsigned int system_size = Size(b);
		  unsigned int local_size = RHS_Contribution.size();
		  PetscInt ierr;

// 			if (BaseType::mCalculateReactionsFlag==false) //if we don't need to calculate reactions
// 			{
				for (unsigned int i_local=0; i_local<local_size; i_local++)
				{
					unsigned int i_global=EquationId[i_local];
					if ( i_global < system_size ) //on "free" DOFs
					{	// ASSEMBLING THE SYSTEM VECTOR
					  VecSetValue(b, i_global, RHS_Contribution[i_local], ADD_VALUES);
						  b[i_global] += RHS_Contribution[i_local];
						
					}
				}
// 			}
// 			else //when the calculation of reactions is needed
// 			{
// 				for (unsigned int i_local=0; i_local<local_size; i_local++)
// 				{
// 					unsigned int i_global=EquationId[i_local];
// 					if ( i_global < BaseType::mEquationSystemSize ) //on "free" DOFs
// 					{	// ASSEMBLING THE SYSTEM VECTOR
// 						b[i_global] += RHS_Contribution[i_local];
// 					}
// 					else //on "fixed" DOFs
// 					{	// Assembling the Vector of REACTIONS
// 						BaseType::mReactionsVector[i_global-BaseType::mEquationSystemSize] -= RHS_Contribution[i_local];
// 					}
// 				}
// 			}
				ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
				ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
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
	return "UBlasSpace";
	}
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
	{
	rOStream << "UBlasSpace";
	}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const
	{
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
      PetscSpace& operator=(PetscSpace const& rOther);

      /// Copy constructor.
      PetscSpace(PetscSpace const& rOther);

        
      ///@}    
        
    }; // Class PetscSpace 



  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
//   inline std::istream& operator >> (std::istream& rIStream, 
// 				    PetscSpace& rThis);

//   /// output stream function
//   inline std::ostream& operator << (std::ostream& rOStream, 
// 				    const PetscSpace& rThis)
//     {
//       rThis.PrintInfo(rOStream);
//       rOStream << std::endl;
//       rThis.PrintData(rOStream);

//       return rOStream;
//     }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_PETSC_SPACE_H_INCLUDED  defined 


