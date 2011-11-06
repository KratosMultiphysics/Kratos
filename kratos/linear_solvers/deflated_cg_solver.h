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
//   Last Modified by:    $Author: mossaiby $
//   Date:                $Date: 2008-12-22 14:46:36 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_DEFLATED_CG_SOLVER_H_INCLUDED )
#define  KRATOS_DEFLATED_CG_SOLVER_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <vector>
#include <set>

// External includes 


// Project includes
#include "includes/define.h"
#include "linear_solvers/iterative_solver.h"
#include "linear_solvers/skyline_lu_factorization_solver.h"
#include "utilities/deflation_utils.h"

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
    class DeflatedCGSolver : public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType>
    {
      public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of DeflatedCGSolver
      KRATOS_CLASS_POINTER_DEFINITION(DeflatedCGSolver);

      typedef IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType> BaseType;

//      typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> LinearSolverType;
  
      typedef typename TSparseSpaceType::MatrixType SparseMatrixType;
  
      typedef typename TSparseSpaceType::VectorType SparseVectorType;
  
      typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
  
      typedef typename TDenseSpaceType::VectorType DenseVectorType;
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      DeflatedCGSolver() {}

      DeflatedCGSolver(double NewMaxTolerance, bool assume_constant_structure, int max_reduced_size) : 
	  BaseType(NewMaxTolerance)
	  ,massume_constant_structure(assume_constant_structure)
	  ,mmax_reduced_size(max_reduced_size) {}

      DeflatedCGSolver(double NewMaxTolerance, unsigned int NewMaxIterationsNumber, bool assume_constant_structure, int max_reduced_size) : BaseType(NewMaxTolerance, NewMaxIterationsNumber),massume_constant_structure(assume_constant_structure),mmax_reduced_size(max_reduced_size) {}

      DeflatedCGSolver(double NewMaxTolerance, unsigned int NewMaxIterationsNumber, 
		       typename TPreconditionerType::Pointer pNewPreconditioner, bool assume_constant_structure, int max_reduced_size) : 
      BaseType(NewMaxTolerance, NewMaxIterationsNumber, pNewPreconditioner)
      ,massume_constant_structure(assume_constant_structure)
,mmax_reduced_size(max_reduced_size)
      {}

      /// Copy constructor.
      DeflatedCGSolver(const DeflatedCGSolver& Other) : BaseType(Other) {}


      /// Destructor.
      virtual ~DeflatedCGSolver(){}
      

      ///@}
      ///@name Operators 
      ///@{
      
      /// Assignment operator.
      DeflatedCGSolver& operator=(const DeflatedCGSolver& Other)
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
      bool Solve(SparseMatrixType& rA, SparseVectorType& rX, SparseVectorType& rB)
	{
	  KRATOS_WATCH("entered in solver")
	  if(this->IsNotConsistent(rA, rX, rB))
	    return false;
	  
// 	  GetTimeTable()->Start(Info());

// 	  BaseType::GetPreconditioner()->Initialize(rA,rX,rB);
//  	  BaseType::GetPreconditioner()->ApplyInverseRight(rX);
// 	  BaseType::GetPreconditioner()->ApplyLeft(rB);
KRATOS_WATCH("ln173")
	  bool is_solved = IterativeSolve(rA,rX,rB);

//  	  BaseType::GetPreconditioner()->Finalize(rX);
KRATOS_WATCH("ln177")
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

	  std::cout << "************ DeflatedCGSolver::Solve(SparseMatrixType&, DenseMatrixType&, DenseMatrixType&) not defined! ************" << std::endl;

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
	  buffer << "Deflated Conjugate gradient linear solver with " << BaseType::GetPreconditioner()->Info();
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
      int mmax_reduced_size;
      bool massume_constant_structure;
      std::vector<int> mw;
      SparseMatrixType mAdeflated;
      
      //typename LinearSolverType::Pointer  mpLinearSolver;
        
      ///@} 
      ///@name Private Operators
      ///@{ 
        
        
      ///@} 
      ///@name Private Operations
      ///@{ 
	
      
      bool IterativeSolve(SparseMatrixType& rA, SparseVectorType& rX, SparseVectorType& rB)
      {
	const int full_size = TSparseSpaceType::Size(rX);

	//construct "coloring" structure and fill reduced matrix structure 
	//note that this has to be done only once if the matrix structure is preserved
	if(massume_constant_structure == false || mw.size()==0)
	{
	  std::cout << "constructing the W matrix and the reduced size one" << std::endl;
	    DeflationUtils::ConstructW(mmax_reduced_size,rA, mw, mAdeflated);
	}

	//fill reduced matrix mmAdeflated
	DeflationUtils::FillDeflatedMatrix(rA,mw,mAdeflated);

	std::size_t reduced_size = mAdeflated.size1();
	
	std::cout << "within solver: full size=" << full_size << " reduced_size="<<reduced_size << " deflation factor = " << double(full_size)/double(reduced_size) << std::endl;
	
	// To save some time, we do the factorization once, and do the solve several times.
	// When this is available through the LinearSolver interface, replace this.

	LUSkylineFactorization<TSparseSpaceType, TDenseSpaceType> Factorization;

	//mpLinearSolver = LinearSolverType::Pointer(new SkylineLUFactorizationSolver<TSparseSpaceType, TDenseSpaceType>);

	Factorization.copyFromCSRMatrix(mAdeflated);
        Factorization.factorize();

	std::cout << "********** Factorization done!" << std::endl;

	SparseVectorType r(full_size), t(full_size), d(full_size), p(full_size), q(full_size);
	SparseVectorType th(reduced_size), dh(reduced_size);

	// r = rA * rX
	TSparseSpaceType::Mult(rA, rX, r);

	// r = rB - r
	TSparseSpaceType::ScaleAndAdd(1.00, rB, -1.00, r);

	std::cout << "********** ||r|| = " << TSparseSpaceType::TwoNorm(r) << std::endl;

	// th = W^T * r -> form reduced problem
	DeflationUtils::ApplyWtranspose(mw,r,th);
// 	TSparseSpaceType::TransposeMult(W, r, th);

	// Solve mAdeflated * th = dh
	Factorization.backForwardSolve(reduced_size, th, dh);

	// t = W * dh -> transfer reduced problem to large scale one
	DeflationUtils::ApplyW(mw,dh,t);
// 	TSparseSpaceType::Mult(W, dh, t);

	// rX = rX + t
	TSparseSpaceType::ScaleAndAdd(1.00, t, 1.00, rX);

	//r = rA * rX
	TSparseSpaceType::Mult(rA, rX, r);

	// r = B - r
	TSparseSpaceType::ScaleAndAdd(1.00, rB, -1.00, r);

	// t = A * r
	//TSparseSpaceType::Mult(rA, r, t);
	this->PreconditionedMult(rA, r, t);

	// th = W^T * t
	DeflationUtils::ApplyWtranspose(mw,t,th);
	//	TSparseSpaceType::TransposeMult(W, t, th);
	
	// Solve mAdeflated * th = dh
	Factorization.backForwardSolve(reduced_size, th, dh);

	// p = W * dh
	DeflationUtils::ApplyW(mw,dh,p);
// 	TSparseSpaceType::Mult(W, dh, p);

	// p = r - p
	TSparseSpaceType::ScaleAndAdd(1.00, r, -1.00, p);
	
	// Iteration counter
	BaseType::mIterationsNumber = 0;
	
	BaseType::mBNorm = TSparseSpaceType::TwoNorm(rB);

	double roh0 = TSparseSpaceType::Dot(r, r);
	double roh1 = roh0;
	double beta = 0;

	if(fabs(roh0) < 1.0e-30) //modification by Riccardo
//	if(roh0 == 0.00)

	  return false;
	    
	do
	  {
	    TSparseSpaceType::Mult(rA, p, q);

	    double pq = TSparseSpaceType::Dot(p, q);

	    //std::cout << "********** pq = " << pq << std::endl;

	    //if(pq == 0.00)
	    if(fabs(pq) <= 1.0e-30)
	      break;

	    double alpha = roh0 / pq;

	    TSparseSpaceType::ScaleAndAdd(alpha, p, 1.00, rX);
	    TSparseSpaceType::ScaleAndAdd(-alpha, q, 1.00, r);

	    roh1 = TSparseSpaceType::Dot(r, r);

	    beta = (roh1 / roh0);

	    // t = A * r
	    //TSparseSpaceType::Mult(rA, r, t);
	    TSparseSpaceType::Mult(rA, r, t);

	    // th = W^T * t
	    DeflationUtils::ApplyWtranspose(mw,t,th);
// 	    TSparseSpaceType::TransposeMult(W, t, th);

	    // Solve mAdeflated * th = dh
	    Factorization.backForwardSolve(reduced_size, th, dh);

	    // t = W * dh
	    DeflationUtils::ApplyW(mw,dh,t);
	    // TSparseSpaceType::Mult(W, dh, t);

	    // t = r - t
	    TSparseSpaceType::ScaleAndAdd(1.00, r, -1.00, t);

	    // p = beta * p + t
	    TSparseSpaceType::ScaleAndAdd(1.00, t, beta, p);

	    roh0 = roh1;

	    BaseType::mResidualNorm = sqrt(roh1);

	    BaseType::mIterationsNumber++;

// 	    if (BaseType::mIterationsNumber % 100 == 0)
	      std::cout << "********** iteration = " << BaseType::mIterationsNumber << ", resnorm = " << BaseType::mResidualNorm << std::endl;

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
        
    }; // Class DeflatedCGSolver 

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
				      DeflatedCGSolver<TSparseSpaceType, TDenseSpaceType, 
				      TPreconditionerType, TReordererType>& rThis)
    {
		return IStream;
    }

  /// output stream function
  template<class TSparseSpaceType, class TDenseSpaceType, 
    class TPreconditionerType, 
    class TReordererType>
  inline std::ostream& operator << (std::ostream& OStream, 
				    const DeflatedCGSolver<TSparseSpaceType, TDenseSpaceType, 
				      TPreconditionerType, TReordererType>& rThis)
    {
      rThis.PrintInfo(OStream);
      OStream << std::endl;
      rThis.PrintData(OStream);

      return OStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_DEFLATED_CG_SOLVER_H_INCLUDED  defined 


