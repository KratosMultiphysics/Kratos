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
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_DIAGONAL_PRECONDITIONER_H_INCLUDED )
#define  KRATOS_DIAGONAL_PRECONDITIONER_H_INCLUDED



// System includes 


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
  
  ///@name  Preconditioners 
  ///@{ 
  
  /// DiagonalPreconditioner class. 
  /** DiagonalPreconditioner for linesr system solvers.  
   */
  template<class TSparseSpaceType, class TDenseSpaceType>
    class DiagonalPreconditioner : public Preconditioner<TSparseSpaceType, TDenseSpaceType>
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Counted pointer of DiagonalPreconditioner
      typedef boost::shared_ptr<DiagonalPreconditioner> Pointer;

      typedef  Preconditioner<TSparseSpaceType, TDenseSpaceType> BaseType;
  
      typedef typename TSparseSpaceType::DataType DataType;
  
      typedef typename TSparseSpaceType::MatrixType SparseMatrixType;
  
      typedef typename TSparseSpaceType::VectorType VectorType;
  
      typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      DiagonalPreconditioner(){}

      /// Copy constructor.
      DiagonalPreconditioner(const DiagonalPreconditioner& Other) 
	: BaseType(Other), mDiagonal(Other.mDiagonal), mTemp(Other.mTemp) {} 

      /// Destructor.
      virtual ~DiagonalPreconditioner(){}
      

      ///@}
      ///@name Operators 
      ///@{
      
      /// Assignment operator.
      DiagonalPreconditioner& operator=(const DiagonalPreconditioner& Other)
	{
	  BaseType::operator=(Other);
	  mDiagonal = Other.mDiagonal;
	  mTemp = Other.mTemp;
	  return *this;
	}

      
      ///@}
      ///@name Operations
      ///@{
      

      /** DiagonalPreconditioner Initialize
	  Initialize preconditioner for linear system rA*rX=rB
	  @param rA  system matrix.
	  @param rX Unknows vector
	  @param rB Right side linear system of equations.
      */
      void Initialize(SparseMatrixType& rA, VectorType& rX, VectorType& rB) 
	{
	  mDiagonal.resize(TSparseSpaceType::Size(rX));
	  mTemp.resize(TSparseSpaceType::Size(rX));

	  unsigned int i;

 	  const DataType zero = DataType();
// 	  const DataType one = DataType(1.00);

	  #pragma omp parallel for private(i)
	  for(i = 0 ; i < rA.size1() ; ++i)
	  {
		double diag_Aii = rA(i,i);
	    if(diag_Aii != zero)
	      mDiagonal[i] = 1.00 / sqrt(fabs(diag_Aii));
	    else
	      KRATOS_ERROR(std::logic_error,"zero found in the diagonal. Diagonal preconditioner can not be used","");
	  }
// 	      mDiagonal[i] = one; 

/* 	  std::cout << "mDiagonal : " << mDiagonal << std::endl; */
	  
	  /* for(i = 0 ; i < rA.RowsNumber() ; ++i)
	     for(j = 0 ; j < rA.RowsNumber() ; ++j) */
	      
	  
	}

      void Initialize(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB) 
	{
	  BaseType::Initialize(rA, rX, rB);
	}

      void Mult(SparseMatrixType& rA, VectorType& rX, VectorType& rY)
      {
	unsigned int i;
	#pragma omp parallel for private(i)
	for(i = 0 ; i < TSparseSpaceType::Size(rX) ; ++i)
	    mTemp[i] = rX[i] * mDiagonal[i];
	TSparseSpaceType::Mult(rA,mTemp, rY);
	ApplyLeft(rY);
      }
      
      void TransposeMult(SparseMatrixType& rA, VectorType& rX, VectorType& rY)
      {
	unsigned int i;
	#pragma omp parallel for private(i)
	for(i = 0 ; i < TSparseSpaceType::Size(rX) ; ++i)
	    mTemp[i] = rX[i] * mDiagonal[i];
	TSparseSpaceType::TransposeMult(rA,mTemp, rY);
	ApplyRight(rY);
      }
      
      VectorType& ApplyLeft(VectorType& rX)
	{
		unsigned int i;
		#pragma omp parallel for private(i)
	  for(i = 0 ; i < TSparseSpaceType::Size(rX) ; ++i)
	    rX[i] *= mDiagonal[i];
	  
	  return rX;
	}
      
      VectorType& ApplyRight(VectorType& rX)
	{
unsigned int i;
	#pragma omp parallel for private(i)
	  for(i = 0 ; i < TSparseSpaceType::Size(rX) ; ++i)
	    rX[i] *= mDiagonal[i];
	  
	  return rX;
        }
      
      /** DiagonalPreconditioner transpose solver.
	  Solving tranpose preconditioner system M^T*x=y, where m^T means transpose.
	  @param rMatrix   DiagonalPreconditioner system matrix.
	  @param rXVector  Unknows of preconditioner suystem
	  @param rYVector  Right side of preconditioner system.
      */    
      VectorType& ApplyTransposeLeft(VectorType& rX)
	{
unsigned int i;
	#pragma omp parallel for private(i)
	  for(i = 0 ; i < TSparseSpaceType::Size(rX) ; ++i)
	    rX[i] *= mDiagonal[i];
	  
	  return rX;
	}
      
      VectorType& ApplyTransposeRight(VectorType& rX)
	{
unsigned int i;
	#pragma omp parallel for private(i)
	  for(i = 0 ; i < TSparseSpaceType::Size(rX) ; ++i)
	    rX[i] *= mDiagonal[i];
	  
	  return rX;
	}
      
      VectorType& ApplyInverseRight(VectorType& rX)
	{
// 	  const DataType zero = DataType();

unsigned int i;
	#pragma omp parallel for private(i)
	  for(i = 0 ; i < TSparseSpaceType::Size(rX) ; ++i)
/*	    if(mDiagonal[i] != zero)*/
	      rX[i] /= mDiagonal[i];
	  
	  return rX;
	}
      
      VectorType& Finalize(VectorType& rX)
	{
	unsigned int i;
	#pragma omp parallel for private(i)
	  for(i = 0 ; i < TSparseSpaceType::Size(rX) ; ++i)
	    rX[i] *= mDiagonal[i];
	  
	  return rX;
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
      
      /// Return information about this object.
       virtual std::string  Info() const
	{
	  return "Diagonal preconditioner";
	}
      
      /// Print information about this object.
      virtual void  PrintInfo(std::ostream& OStream) const
	{
	  OStream << "Diagonal preconditioner";
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

      VectorType mDiagonal;
        
      VectorType mTemp;
        
        
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
      
        
      ///@}    
        
    }; // Class DiagonalPreconditioner 
  
  ///@} 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  template<class TSparseSpaceType, class TDenseSpaceType>
   inline std::istream& operator >> (std::istream& IStream, 
 				    DiagonalPreconditioner<TSparseSpaceType, TDenseSpaceType>& rThis)
    {
		return IStream;
    }

  /// output stream function
  template<class TSparseSpaceType, class TDenseSpaceType>
   inline std::ostream& operator << (std::ostream& OStream, 
 				    const DiagonalPreconditioner<TSparseSpaceType, TDenseSpaceType>& rThis)
   {
      rThis.PrintInfo(OStream);
      OStream << std::endl;
      rThis.PrintData(OStream);

      return OStream;
   }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_DIAGONAL_PRECONDITIONER_H_INCLUDED  defined 

