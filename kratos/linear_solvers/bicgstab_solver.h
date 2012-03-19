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
//   Date:                $Date: 2008-03-05 09:39:14 $
//   Revision:            $Revision: 1.4 $
//
//


#if !defined(KRATOS_BICGSTAB_SOLVER_H_INCLUDED )
#define  KRATOS_BICGSTAB_SOLVER_H_INCLUDED



// System includes 


// External includes 
#include "boost/smart_ptr.hpp"


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
    class BICGSTABSolver : public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType>
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Counted pointer of BICGSTABSolver
      typedef boost::shared_ptr<BICGSTABSolver> Pointer;

      typedef IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType> BaseType; 
  
      typedef typename TSparseSpaceType::MatrixType SparseMatrixType;
  
      typedef typename TSparseSpaceType::VectorType VectorType;
  
      typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      BICGSTABSolver(){}

      BICGSTABSolver(double NewTolerance) : BaseType(NewTolerance){}

      BICGSTABSolver(double NewTolerance, unsigned int NewMaxIterationsNumber) : BaseType(NewTolerance, NewMaxIterationsNumber){}

      BICGSTABSolver(double NewMaxTolerance, unsigned int NewMaxIterationsNumber, typename TPreconditionerType::Pointer pNewPreconditioner) : 
	BaseType(NewMaxTolerance, NewMaxIterationsNumber, pNewPreconditioner){}

      /// Copy constructor.
    BICGSTABSolver(const BICGSTABSolver& Other) : BaseType(Other){}

      /// Destructor.
      virtual ~BICGSTABSolver(){}
      

      ///@}
      ///@name Operators 
      ///@{
      
      /// Assignment operator.
      BICGSTABSolver& operator=(const BICGSTABSolver& Other)
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

	  //GetTimeTable()->Start(Info());

// KRATOS_WATCH("ln161");
	  BaseType::GetPreconditioner()->Initialize(rA,rX,rB);
// KRATOS_WATCH("ln163");
 	  BaseType::GetPreconditioner()->ApplyInverseRight(rX);
// KRATOS_WATCH("ln165");
	  BaseType::GetPreconditioner()->ApplyLeft(rB);
// KRATOS_WATCH("ln167");
	  bool is_solved = IterativeSolve(rA,rX,rB);
// KRATOS_WATCH("ln169");

 	  BaseType::GetPreconditioner()->Finalize(rX);

	  //GetTimeTable()->Stop(Info());

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
	  //GetTimeTable()->Start(Info());

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

	  //GetTimeTable()->Stop(Info());

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

      /// Return information about this object.
      virtual std::string Info() const
	{
	  std::stringstream buffer;
	  buffer << "Biconjugate gradient stabilized linear solver with " << BaseType::GetPreconditioner()->Info();
	  return  buffer.str();
	}
      
      /// Print information about this object.
      void  PrintInfo(std::ostream& OStream) const
	{
	  OStream << "Biconjugate gradient stabilized linear solver with ";
	  BaseType::GetPreconditioner()->PrintInfo(OStream);
	}

      /// Print object's data.
      void  PrintData(std::ostream& OStream) const 
	{
	  BaseType::PrintData(OStream);
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
// KRATOS_WATCH("ln316");	
	BaseType::mIterationsNumber = 0;
    
	VectorType r(size);
// KRATOS_WATCH(r.size());
// KRATOS_WATCH("ln319");	
// KRATOS_WATCH(rA.size1());
// KRATOS_WATCH(rA.size2());
// KRATOS_WATCH(r.size());
// KRATOS_WATCH(rX.size());
// KRATOS_WATCH(rB.size());

	this->PreconditionedMult(rA,rX,r);
// KRATOS_WATCH("ln321");	
	TSparseSpaceType::ScaleAndAdd(1.00, rB, -1.00, r);
// KRATOS_WATCH("ln322");
	BaseType::mBNorm = TSparseSpaceType::TwoNorm(rB);
// KRATOS_WATCH("ln324");
	VectorType p(r);
	VectorType s(size);
	VectorType q(size);

 	VectorType rs(r); 
 	VectorType qs(size); 
         
	double roh0 = TSparseSpaceType::Dot(r, rs);
	double roh1 = roh0;
	double alpha = 0.00;
	double beta = 0.00;
	double omega = 0.00;
// KRATOS_WATCH("ln337");	
// 	if(roh0 < 1e-30) //we start from the real solution
// 		return  BaseType::IsConverged();

	do
	  {
	    this->PreconditionedMult(rA,p,q);
// KRATOS_WATCH("ln344");
	    alpha = roh0 / TSparseSpaceType::Dot(rs,q);
        
	    TSparseSpaceType::ScaleAndAdd(1.00, r, -alpha, q, s);
// KRATOS_WATCH("ln348");
	    this->PreconditionedMult(rA,s,qs);

	    omega = TSparseSpaceType::Dot(qs,qs);

	    //if(omega == 0.00)
	    if(fabs(omega) <= 1.0e-40)
	      break;
// KRATOS_WATCH("ln356");
	    omega = TSparseSpaceType::Dot(qs,s) / omega;

	    TSparseSpaceType::ScaleAndAdd(alpha, p, 1.00, rX);
	    TSparseSpaceType::ScaleAndAdd(omega, s, 1.00, rX);
	    TSparseSpaceType::ScaleAndAdd(-omega, qs, 1.00, s, r);

	    roh1 = TSparseSpaceType::Dot(r,rs);

	    //if((roh0 == 0.00) || (omega == 0.00))
	    if((fabs(roh0) <= 1.0e-40) || (fabs(omega) <= 1.0e-40))
	      break;
	    
	    beta = (roh1 * alpha) / (roh0 * omega);
// KRATOS_WATCH("ln370");	    
	    TSparseSpaceType::ScaleAndAdd(1.00, p, -omega, q);
	    TSparseSpaceType::ScaleAndAdd(1.00, r, beta, q, p);
	      
	    roh0 = roh1;
        
		BaseType::mResidualNorm =TSparseSpaceType::TwoNorm(r);
		BaseType::mIterationsNumber++;

	  } while(BaseType::IterationNeeded());
	  
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
        
    }; // Class BICGSTABSolver 

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
				      BICGSTABSolver<TSparseSpaceType, TDenseSpaceType, 
				      TPreconditionerType, TReordererType>& rThis)
    {
		return IStream;
    }

  /// output stream function
  template<class TSparseSpaceType, class TDenseSpaceType, 
    class TPreconditionerType, 
    class TReordererType>
  inline std::ostream& operator << (std::ostream& OStream, 
				    const BICGSTABSolver<TSparseSpaceType, TDenseSpaceType, 
				      TPreconditionerType, TReordererType>& rThis)
    {
      rThis.PrintInfo(OStream);
      OStream << std::endl;
      rThis.PrintData(OStream);

      return OStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_BICGSTAB_SOLVER_H_INCLUDED  defined 


