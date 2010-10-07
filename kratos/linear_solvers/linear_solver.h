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


#if !defined(KRATOS_LINEAR_SOLVER_H_INCLUDED )
#define  KRATOS_LINEAR_SOLVER_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <sstream>
#include <cstddef>


// External includes 


// Project includes
#include "includes/define.h"
#include "reorderer.h"


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
  
  /// Base class for all the linear solvers in Kratos.
  /** This class define the general interface for the linear solvers in Kratos.
      There is three template parameter: 
      - TSparseSpaceType which specify type
        of the unknowns, coefficients, sparse matrix, vector of
	unknowns, right hand side vector and their respective operators.
      - TDenseMatrixType which specify type of the
        matrices used as temporary matrices or multi solve unknowns and
	right hand sides and their operators.  
      - TReordererType which specify type of the Orderer that performs the reordering of matrix to optimize the solution.
  */
  template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
    class LinearSolver
    {
      public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of LinearSolver
      KRATOS_CLASS_POINTER_DEFINITION(LinearSolver);
  
      typedef typename TSparseSpaceType::MatrixType SparseMatrixType;
  
      typedef typename TSparseSpaceType::VectorType VectorType;
  
      typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
  
      typedef typename TDenseSpaceType::VectorType DenseVectorType;

      typedef std::size_t  SizeType;
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      LinearSolver() : mpReorderer(new TReordererType()){}

      /// Constructor with specific reorderer.
      LinearSolver(TReordererType NewReorderer) : mpReorderer(NewReorderer){}

      /// Copy constructor.
      LinearSolver(const LinearSolver& Other) : mpReorderer(Other.mpReorderer){}

      /// Destructor.
      virtual ~LinearSolver(){}
      

      ///@}
      ///@name Operators 
      ///@{
      
      /// Assignment operator.
      LinearSolver& operator=(const LinearSolver& Other)
	{
	  mpReorderer = Other.mpReorderer;

	  return *this;
	}

      
      ///@}
      ///@name Operations
      ///@{
      
      virtual void Initialize(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
	{
	  mpReorderer->Initialize(rA, rX, rB);
	}
      
      /** Normal solve method.
	  Solves the linear system Ax=b and puts the result on SystemVector& rX. 
	  rVectorx is also th initial guess for iterative methods.
	  @param rA. System matrix
	  @param rX. Solution vector. it's also the initial 
	  guess for iterative linear solvers.
 	  @param rB. Right hand side vector.
      */
      virtual bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
	{
	  return false;
	}
      
      
      /** Multi solve method for solving a set of linear systems with same coefficient matrix.
	  Solves the linear system Ax=b and puts the result on SystemVector& rX. 
	  rVectorx is also th initial guess for iterative methods.
	  @param rA. System matrix
	  @param rX. Solution vector. it's also the initial 
	  guess for iterative linear solvers.
 	  @param rB. Right hand side vector.
      */
      virtual bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB)
	{
	  return false;
	}

      /** Eigenvalue and eigenvector solve method for derived eigensolvers */
      virtual  void Solve(SparseMatrixType& K, 
			  SparseMatrixType& M,
			  DenseVectorType& Eigenvalues,
			  DenseMatrixType& Eigenvectors){}
      
      ///@}
      ///@name Access
      ///@{ 
      
      virtual typename TReordererType::Pointer GetReorderer(void)
	{
	  return mpReorderer;
	}

      virtual void SetReorderer(typename TReordererType::Pointer pNewReorderer)
	{
	  mpReorderer = pNewReorderer;
	}

      virtual void SetTolerance(double NewTolerance)
	{
          std::cout << "WARNING: Accessed base function Kratos::LinearSolver::SetTolerance(double). This does nothing !" << std::endl;
	}

      virtual double GetTolerance()
	{
          std::cout << "WARNING: Accessed base function Kratos::LinearSolver::GetTolerance(). No tolerance defined, returning 0 !" << std::endl ;
	  return 0;
	}
      
      
      ///@}
      ///@name Inquiry
      ///@{
      
      virtual bool IsConsistent(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
	{
	  const SizeType size = TSparseSpaceType::Size1(rA);

	  return ((size ==  TSparseSpaceType::Size2(rA)) && 
		  (size ==  TSparseSpaceType::Size(rX)) &&
		  (size ==  TSparseSpaceType::Size(rB))); 
	}
      
      virtual bool IsConsistent(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB)
	{
	  const SizeType size = TSparseSpaceType::Size1(rA);

	  return ((size ==  TSparseSpaceType::Size2(rA)) && 
		  (size ==  TDenseSpaceType::Size1(rX)) &&
		  (size ==  TDenseSpaceType::Size1(rB)) && 
		  (TDenseSpaceType::Size2(rX) == TDenseSpaceType::Size2(rB)));
	}
      
      
      virtual bool IsNotConsistent(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
	{
	  return (!IsConsistent(rA, rX, rB));
	}

      virtual bool IsNotConsistent(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB)
	{
	  return (!IsConsistent(rA, rX, rB));
	}
      
      ///@}      
      ///@name Input and output
      ///@{

      /// Turn back information as a string.
      virtual std::string Info() const
	{
	  return "Linear solver";
	}
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
	{
	   rOStream << "Linear solver";
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
        
      /// A counted pointer to the reorderer object.
      typename TReordererType::Pointer mpReorderer;

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
        
    }; // Class LinearSolver 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
    inline std::istream& operator >> (std::istream& IStream, 
				      LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis)
    {
      return IStream;
    }

    /// output stream function
  template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
    inline std::ostream& operator << (std::ostream& rOStream, 
				      const LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_LINEAR_SOLVER_H_INCLUDED  defined 




