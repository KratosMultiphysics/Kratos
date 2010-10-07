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


#if !defined(KRATOS_CG_SOLVER_H_INCLUDED )
#define  KRATOS_CG_SOLVER_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


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
    class TPreconditionerType = Preconditioner<TSparseSpaceType, TDenseSpaceType>, 
    class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
    class CGSolver : public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType>
    {
      public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of CGSolver
      KRATOS_CLASS_POINTER_DEFINITION(CGSolver);

      typedef IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType> BaseType; 
  
      typedef typename TSparseSpaceType::MatrixType SparseMatrixType;
  
      typedef typename TSparseSpaceType::VectorType VectorType;
  
      typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      CGSolver(){}

      CGSolver(double NewMaxTolerance) : BaseType(NewMaxTolerance){}

      CGSolver(double NewMaxTolerance, unsigned int NewMaxIterationsNumber) : BaseType(NewMaxTolerance, NewMaxIterationsNumber){}

      CGSolver(double NewMaxTolerance, unsigned int NewMaxIterationsNumber, typename TPreconditionerType::Pointer pNewPreconditioner) : 
      BaseType(NewMaxTolerance, NewMaxIterationsNumber, pNewPreconditioner){}

      /// Copy constructor.
      CGSolver(const CGSolver& Other) : BaseType(Other) {}


      /// Destructor.
      virtual ~CGSolver(){}
      

      ///@}
      ///@name Operators 
      ///@{
      
      /// Assignment operator.
      CGSolver& operator=(const CGSolver& Other)
      {
        BaseType::operator=(Other);
	return *this;
      }
      
      ///@}
      ///@name Operations
      ///@{
      
      /** Normal solve method.
	  Solves the linear system Ax=b and puts the result on SystemVector& rX. 
	  rX is also th initial guess for iterative methods.
	  @param rA. System matrix
	  @param rX. Solution vector. it's also the initial 
	  guess for iterative linear solvers.
 	  @param rB. Right hand side vector.
      */
      bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
	{
	  if(this->IsNotConsistent(rA, rX, rB))
	    return false;
	  
// 	  GetTimeTable()->Start(Info());

	  BaseType::GetPreconditioner()->Initialize(rA,rX,rB);
 	  BaseType::GetPreconditioner()->ApplyInverseRight(rX);
	  BaseType::GetPreconditioner()->ApplyLeft(rB);

	  bool is_solved = IterativeSolve(rA,rX,rB);

 	  BaseType::GetPreconditioner()->Finalize(rX);

// 	  GetTimeTable()->Stop(Info());

	  return is_solved;
	}
      
      
      /** Multi solve method for solving a set of linear systems with same coefficient matrix.
	  Solves the linear system Ax=b and puts the result on SystemVector& rX. 
	  rX is also th initial guess for iterative methods.
	  @param rA. System matrix
	  @param rX. Solution vector. it's also the initial 
	  guess for iterative linear solvers.
 	  @param rB. Right hand side vector.
      */
      bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB)
	{
// 	  GetTimeTable()->Start(Info());

   	  BaseType::GetPreconditioner()->Initialize(rA,rX,rB);

 	  bool is_solved = true;
	  VectorType x(TDenseSpaceType::Size1(rX));
	  VectorType b(TDenseSpaceType::Size1(rB));
	  for(unsigned int i = 0 ; i < TDenseSpaceType::Size2(rX) ; i++)
	    {
	      TDenseSpaceType::GetColumn(i,rX, x);
	      TDenseSpaceType::GetColumn(i,rB, b);
	      
	      BaseType::GetPreconditioner()->ApplyInverseRight(x);
	      BaseType::GetPreconditioner()->ApplyLeft(b);

	      is_solved &= IterativeSolve(rA,x,b);

	      BaseType::GetPreconditioner()->Finalize(x);
	    }

// 	  GetTimeTable()->Stop(Info());

	  return is_solved;
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
	  buffer << "Conjugate gradient linear solver with " << BaseType::GetPreconditioner()->Info();
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
        
        
      ///@} 
      ///@name Private Operators
      ///@{ 
        
        
      ///@} 
      ///@name Private Operations
      ///@{ 
        
      bool IterativeSolve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
      {
	const int size = TSparseSpaceType::Size(rX);
	
	BaseType::mIterationsNumber = 0;
    
	VectorType r(size);
	
	this->PreconditionedMult(rA,rX,r);
	TSparseSpaceType::ScaleAndAdd(1.00, rB, -1.00, r);

	BaseType::mBNorm = TSparseSpaceType::TwoNorm(rB);
    
	VectorType p(r);
	VectorType q(size);
         
	double roh0 = TSparseSpaceType::Dot(r, r);
	double roh1 = roh0;
	double beta = 0;

	if(fabs(roh0) < 1.0e-30) //modification by Riccardo
//	if(roh0 == 0.00)
	  return false;
	    
	do
	  {
	    this->PreconditionedMult(rA,p,q);

	    double pq = TSparseSpaceType::Dot(p,q);

	    //if(pq == 0.00)
	    if(fabs(pq) <= 1.0e-30)
	      break;

	    double alpha = roh0 / pq;
        
	    TSparseSpaceType::ScaleAndAdd(alpha, p, 1.00, rX);
	    TSparseSpaceType::ScaleAndAdd(-alpha, q, 1.00, r);
         
	    roh1 = TSparseSpaceType::Dot(r,r);

	    beta = (roh1 / roh0);
	    TSparseSpaceType::ScaleAndAdd(1.00, r, beta, p);
	      
	    roh0 = roh1;

	    BaseType::mResidualNorm = sqrt(roh1);
	    BaseType::mIterationsNumber++;
	  } while(BaseType::IterationNeeded() && (fabs(roh0) > 1.0e-30)/*(roh0 != 0.00)*/);
	  
	return BaseType::IsConverged();
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
        
    }; // Class CGSolver 

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
				      CGSolver<TSparseSpaceType, TDenseSpaceType, 
				      TPreconditionerType, TReordererType>& rThis)
    {
		return IStream;
    }

  /// output stream function
  template<class TSparseSpaceType, class TDenseSpaceType, 
    class TPreconditionerType, 
    class TReordererType>
  inline std::ostream& operator << (std::ostream& OStream, 
				    const CGSolver<TSparseSpaceType, TDenseSpaceType, 
				      TPreconditionerType, TReordererType>& rThis)
    {
      rThis.PrintInfo(OStream);
      OStream << std::endl;
      rThis.PrintData(OStream);

      return OStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_CG_SOLVER_H_INCLUDED  defined 


