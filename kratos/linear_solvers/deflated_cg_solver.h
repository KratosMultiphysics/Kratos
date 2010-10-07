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
      DeflatedCGSolver() : mSeedPointNo(0) {}

      DeflatedCGSolver(double NewMaxTolerance, int SeedPointNo) : BaseType(NewMaxTolerance), mSeedPointNo(SeedPointNo) {}

      DeflatedCGSolver(double NewMaxTolerance, unsigned int NewMaxIterationsNumber, int SeedPointNo) : BaseType(NewMaxTolerance, NewMaxIterationsNumber), mSeedPointNo(SeedPointNo) {}

      DeflatedCGSolver(double NewMaxTolerance, unsigned int NewMaxIterationsNumber, typename TPreconditionerType::Pointer pNewPreconditioner, int SeedPointNo) : 
      BaseType(NewMaxTolerance, NewMaxIterationsNumber, pNewPreconditioner), mSeedPointNo(SeedPointNo) {}

      /// Copy constructor.
      DeflatedCGSolver(const DeflatedCGSolver& Other) : BaseType(Other), mSeedPointNo(Other.mSeedPointNo) {}


      /// Destructor.
      virtual ~DeflatedCGSolver(){}
      

      ///@}
      ///@name Operators 
      ///@{
      
      /// Assignment operator.
      DeflatedCGSolver& operator=(const DeflatedCGSolver& Other)
      {
        BaseType::operator=(Other);
	mSeedPointNo = Other.mSeedPointNo;
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
      
      //typename LinearSolverType::Pointer  mpLinearSolver;
      int mSeedPointNo;
        
      ///@} 
      ///@name Private Operators
      ///@{ 
        
        
      ///@} 
      ///@name Private Operations
      ///@{ 
        
      bool IterativeSolve(SparseMatrixType& rA, SparseVectorType& rX, SparseVectorType& rB)
      {
	const int size1 = TSparseSpaceType::Size(rX);
	int size2 = size1 / mSeedPointNo;
	if (size1 % mSeedPointNo != 0) size2++;

	std::cout << "********** mSeedPointNo = " << mSeedPointNo << std::endl;
	std::cout << "********** size1 = " << size1 << std::endl;
	std::cout << "********** size2 = " << size2 << std::endl;

	SparseMatrixType W(size1, size2, size1);

	// Building W
	TSparseSpaceType::SetToZero(W);

	std::cout << "********** W set to zero!" << std::endl;

	for (int i = 0; i < size2; i++)
	  for (int j = 0; j < mSeedPointNo; j++)
	    if (i * mSeedPointNo + j < size1)
	      W(i * mSeedPointNo + j, i) = 1;

	std::cout << "********** W built!" << std::endl;

	// Use TSparseSpaceType::size_type?
	std::vector<std::size_t> w(size1);

	// Building w
	for (int i = 0; i < size2; i++)
	  for (int j = 0; j < mSeedPointNo; j++)
	    if (i * mSeedPointNo + j < size1)
	      w[i * mSeedPointNo + j] = i;

	std::cout << "********** w built!" << std::endl;

	// Non-zero structure of Ah
	std::vector<std::set<std::size_t> > AhNZ(size2);

	// Loop over non-zero structure of A and build non-zero structure of Ah
	typename SparseMatrixType::iterator1 a_iterator = rA.begin1();

	for (int i = 0; i < size1; i++)
	{
		#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
		for (typename SparseMatrixType::iterator2 row_iterator = a_iterator.begin() ;
		row_iterator != a_iterator.end() ; ++row_iterator) 
		{
		#else
		for (typename SparseMatrixType::iterator2 row_iterator = begin(a_iterator,
			boost::numeric::ublas::iterator1_tag());
			row_iterator != end(a_iterator,
			boost::numeric::ublas::iterator1_tag()); ++row_iterator )
		{
		#endif
			AhNZ[w[a_iterator.index1()]].insert(w[row_iterator.index2()]);
		}

	   a_iterator++;
	}

	std::cout << "********** NZS built!" << std::endl;

	// Count the number of non-zeros in Ah
	int NZ = 0;
	for (int i = 0; i < size2; i++)
		NZ += AhNZ[i].size();

	std::cout << "********** NZ = " << NZ << std::endl;

	SparseMatrixType Ah(size2, size2, NZ);

	// Insert the non-zero structure into Ah
	for(int i = 0 ; i < size2 ; i++)
	{
		for(std::set<std::size_t>::iterator j = AhNZ[i].begin() ; j != AhNZ[i].end() ; j++)
		{
			Ah.push_back(i,*j, 0.00);
		}
	}

	// Now building Ah
	a_iterator = rA.begin1();

	for (int i = 0; i < size1; i++)
	{
		#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
		for (typename SparseMatrixType::iterator2 row_iterator = a_iterator.begin() ;
		row_iterator != a_iterator.end() ; ++row_iterator) 
		{
		#else
		for (typename SparseMatrixType::iterator2 row_iterator = begin(a_iterator,
			boost::numeric::ublas::iterator1_tag());
			row_iterator != end(a_iterator,
			boost::numeric::ublas::iterator1_tag()); ++row_iterator )
		{
		#endif
			Ah(w[a_iterator.index1()], w[row_iterator.index2()]) += *row_iterator;
		}

	   a_iterator++;
	}
	
	std::cout << "********** W^T * A * W built!" << std::endl;

	// To save some time, we do the factorization once, and do the solve several times.
	// When this is available through the LinearSolver interface, replace this.

	LUSkylineFactorization<TSparseSpaceType, TDenseSpaceType> Factorization;

	//mpLinearSolver = LinearSolverType::Pointer(new SkylineLUFactorizationSolver<TSparseSpaceType, TDenseSpaceType>);

	Factorization.copyFromCSRMatrix(Ah);
        Factorization.factorize();

	std::cout << "********** Factorization done!" << std::endl;

	SparseVectorType r(size1), t(size1), d(size1), p(size1), q(size1);
	SparseVectorType th(size2), dh(size2);

	// r = rA * rX
	this->PreconditionedMult(rA, rX, r);

	// r = rB - r
	TSparseSpaceType::ScaleAndAdd(1.00, rB, -1.00, r);

	std::cout << "********** ||r|| = " << TSparseSpaceType::TwoNorm(r) << std::endl;

	// th = W^T * r
	TSparseSpaceType::TransposeMult(W, r, th);

	// Solve Ah * th = dh
	Factorization.backForwardSolve(size2, th, dh);

	// t = W * dh
	TSparseSpaceType::Mult(W, dh, t);

	// rX = rX + t
	TSparseSpaceType::ScaleAndAdd(1.00, t, 1.00, rX);

	//r = rA * rX
	this->PreconditionedMult(rA, rX, r);

	// r = B - r
	TSparseSpaceType::ScaleAndAdd(1.00, rB, -1.00, r);

	// t = A * r
	//TSparseSpaceType::Mult(rA, r, t);
	this->PreconditionedMult(rA, r, t);

	// th = W^T * t
	TSparseSpaceType::TransposeMult(W, t, th);
	
	// Solve Ah * th = dh
	Factorization.backForwardSolve(size2, th, dh);

	// p = W * dh
	TSparseSpaceType::Mult(W, dh, p);

	// p = r - p
	TSparseSpaceType::ScaleAndAdd(1.00, r, -1.00, p);
	
	// Iteration counter
	BaseType::mIterationsNumber = 0;
	
	BaseType::mBNorm = TSparseSpaceType::TwoNorm(rB);

	double roh0 = TSparseSpaceType::Dot(r, r);
	double roh1 = roh0;
	double beta = 0;

	std::cout << "********** Just before the loop! " << (fabs(roh0) < 1.0e-30) << std::endl;

	if(fabs(roh0) < 1.0e-30) //modification by Riccardo
//	if(roh0 == 0.00)

	  return false;
	    
	do
	  {
	    this->PreconditionedMult(rA, p, q);

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
	    this->PreconditionedMult(rA, r, t);

	    // th = W^T * t
	    TSparseSpaceType::TransposeMult(W, t, th);

	    // Solve Ah * th = dh
	    Factorization.backForwardSolve(size2, th, dh);

	    // t = W * dh
	    TSparseSpaceType::Mult(W, dh, t);

	    // t = r - t
	    TSparseSpaceType::ScaleAndAdd(1.00, r, -1.00, t);

	    // p = beta * p + t
	    TSparseSpaceType::ScaleAndAdd(1.00, t, beta, p);

	    roh0 = roh1;

	    BaseType::mResidualNorm = sqrt(roh1);

	    BaseType::mIterationsNumber++;

	    if (BaseType::mIterationsNumber % 100 == 0)
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


