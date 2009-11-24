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


#if !defined(KRATOS_GPU_BICGSTAB_SOLVER_H_INCLUDED )
#define  KRATOS_GPU_BICGSTAB_SOLVER_H_INCLUDED


// System includes 

//#include <ctime>

// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "linear_solvers/iterative_solver.h"
#include "gpu_sparse.h"
#include "linear_solvers.h"


using namespace Kratos::GPUSparse;

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
  
    class GPUBICGSTABSolver : public IterativeSolver<UblasSpace<double, CompressedMatrix, Vector>, UblasSpace<double, Matrix, Vector>, Preconditioner<UblasSpace<double, CompressedMatrix, Vector>, UblasSpace<double, Matrix, Vector> >, Reorderer<UblasSpace<double, CompressedMatrix, Vector>, UblasSpace<double, Matrix, Vector> > >
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Counted pointer of GPUBICGSTABSolver
      typedef boost::shared_ptr<GPUBICGSTABSolver> Pointer;

      typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;

      typedef UblasSpace<double, Matrix, Vector> DenseSpaceType;
      
      typedef Reorderer<SparseSpaceType, DenseSpaceType > ReordererType;

      typedef IterativeSolver<SparseSpaceType, DenseSpaceType, Preconditioner<SparseSpaceType, DenseSpaceType>, ReordererType > BaseType; 
  
      typedef SparseSpaceType::MatrixType SparseMatrixType;
  
      typedef SparseSpaceType::VectorType VectorType;
  
      typedef DenseSpaceType::MatrixType DenseMatrixType;
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      GPUBICGSTABSolver(){}

      GPUBICGSTABSolver(double NewTolerance) : BaseType(NewTolerance)
	{
		havePreconditioner = false;
		maxIter = 5000;		
		tol = NewTolerance;
		this->preconditioner = 0;
	}

      GPUBICGSTABSolver(double NewTolerance, unsigned int NewMaxIterationsNumber) : BaseType(NewTolerance, NewMaxIterationsNumber)
	{
		havePreconditioner = false;
		maxIter = NewMaxIterationsNumber;		
		tol = NewTolerance;
		this->preconditioner = 0;
	}

      GPUBICGSTABSolver(double NewMaxTolerance, unsigned int NewMaxIterationsNumber, GPUPreconditioner& _preconditioner) : BaseType(NewMaxTolerance, NewMaxIterationsNumber)
	{
		havePreconditioner = true;
		this->preconditioner = &_preconditioner;
		maxIter = NewMaxIterationsNumber;
		tol = NewMaxTolerance;
      	}

//      GPUBICGSTABSolver(double NewMaxTolerance, unsigned int NewMaxIterationsNumber, typename TPreconditionerType::Pointer pNewPreconditioner) : 
//	BaseType(NewMaxTolerance, NewMaxIterationsNumber, pNewPreconditioner){}

      /// Copy constructor.
    GPUBICGSTABSolver(const GPUBICGSTABSolver& Other) : BaseType(Other){}

      /// Destructor.
      virtual ~GPUBICGSTABSolver(){}
      

      ///@}
      ///@name Operators 
      ///@{
      
      /// Assignment operator.
      GPUBICGSTABSolver& operator=(const GPUBICGSTABSolver& Other)
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
	  if(IsNotConsistent(rA, rX, rB))
	    return false;

	  //GetTimeTable()->Start(Info());

//	  BaseType::GetPreconditioner()->Initialize(rA,rX,rB);
// 	  BaseType::GetPreconditioner()->ApplyInverseRight(rX);
//	  BaseType::GetPreconditioner()->ApplyLeft(rB);

	  bool is_solved = IterativeSolve(rA,rX,rB);

// 	  BaseType::GetPreconditioner()->Finalize(rX);

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

//  	  BaseType::GetPreconditioner()->Initialize(rA,rX,rB);

 	  bool is_solved = true;
//	  VectorType x(TDenseSpaceType::Size1(rX));
//	  VectorType b(TDenseSpaceType::Size1(rB));
//	  for(unsigned int i = 0 ; i < TDenseSpaceType::Size2(rX) ; i++)
//	    {
//	      TDenseSpaceType::GetColumn(i,rX, x);
//	      TDenseSpaceType::GetColumn(i,rB, b);
	      
//	      BaseType::GetPreconditioner()->ApplyInverseRight(x);
//	      BaseType::GetPreconditioner()->ApplyLeft(b);

//	      is_solved &= IterativeSolve(rA,x,b);

//	      BaseType::GetPreconditioner()->Finalize(x);
//	    }

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
	  buffer << "GPU Biconjugate gradient stabilized linear solver with " << BaseType::GetPreconditioner()->Info();
	  return  buffer.str();
	}
      
      /// Print information about this object.
      void  PrintInfo(std::ostream& OStream) const
	{
	  OStream << "GPU Biconjugate gradient stabilized linear solver with ";
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
	//This var control preconditioner 
	bool havePreconditioner;
	GPUPreconditioner *preconditioner;
        size_t maxIter;
	double tol;
        
        
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
      

      
#define START_TIMING(t)		t = std::clock()
#define STOP_TIMING(T, t)	T += std::clock() - t
      
      bool IterativeSolve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
      {
	size_t iterNum = 0;
	BICGSTAB_GPU(rA.size1(), rA.size2(), rA.nnz(), &(rA.value_data() [0]), &(rA.index2_data() [0]), &(rA.index1_data() [0]), rB.size(), &(rX[0]), &(rB[0]), 
			tol, maxIter, BaseType::mBNorm, BaseType::mResidualNorm, iterNum, *(this->preconditioner));
	BaseType::mIterationsNumber = (unsigned int) iterNum;
/*	//std::time_t t1 = 0, t2 = 0, t;

	//KRATOS_WATCH("============================ GPU Solver ============================");
	
	//KRATOS_WATCH(rA.nnz());
	//KRATOS_WATCH(rA.size1());

	// Inputs
	GPUCSRMatrix gA(rA.nnz(), rA.size1(), rA.size2(), &(rA.index2_data() [0]), &(rA.index1_data() [0]), &(rA.value_data() [0]));
	GPU_CHECK(gA.GPU_Allocate());
	GPU_CHECK(gA.Copy(CPU_GPU, false));
//KRATOS_WATCH("327");
	
	GPUVector gX(rX.size(), &(rX[0]));
	GPU_CHECK(gX.GPU_Allocate());
	GPU_CHECK(gX.Copy(CPU_GPU));
//KRATOS_WATCH("332");
	
	GPUVector gB(rB.size(), &(rB[0]));
	GPU_CHECK(gB.GPU_Allocate());
	GPU_CHECK(gB.Copy(CPU_GPU));
//KRATOS_WATCH("337");

	const int size = SparseSpaceType::Size(rX);
	
//KRATOS_WATCH("338");
	BaseType::mIterationsNumber = 0;
    
	//VectorType r(size);
	GPUVector r(size);
	GPU_CHECK(r.GPU_Allocate());
//KRATOS_WATCH("344");	
//	PreconditionedMult(rA,rX,r);
	//SparseSpaceType::Mult(rA,rX,r); // r = rA*rX
	GPU_CHECK(GPU_MatrixVectorMultiply(gA, gX, r));
//KRATOS_WATCH("348");	
	//SparseSpaceType::ScaleAndAdd(1.00, rB, -1.00, r); // r = rB - r
	GPU_CHECK(GPU_VectorScaleAndAdd(1.00, gB, -1.00, r));
//KRATOS_WATCH("351");
	//BaseType::mBNorm = SparseSpaceType::TwoNorm(rB); 
	GPU_CHECK(GPU_VectorNorm2(gB, BaseType::mBNorm));
//KRATOS_WATCH("354");
	//VectorType p(r);
	GPUVector p(size);
	GPU_CHECK(p.GPU_Allocate());
	GPU_CHECK(p.CopyGPUValuesFrom(r));
	
	//VectorType s(size);
	GPUVector s(size);
	GPU_CHECK(s.GPU_Allocate());
	
	//VectorType q(size);
	GPUVector q(size);
	GPU_CHECK(q.GPU_Allocate());

 	//VectorType rs(r); 
 	GPUVector rs(size);
 	GPU_CHECK(rs.GPU_Allocate());
 	GPU_CHECK(rs.CopyGPUValuesFrom(r));
 	
 	//VectorType qs(size); 
 	GPUVector qs(size);
 	GPU_CHECK(qs.GPU_Allocate());
//KRATOS_WATCH("376");         
	//double roh0 = SparseSpaceType::Dot(r, rs);
	double roh0;
	GPU_CHECK(GPU_VectorVectorMultiply(r, rs, roh0));
	
	double roh1 = roh0;
	double alpha = 0.00;
	double beta = 0.00;
	double omega = 0.00;
	
// 	if(roh0 < 1e-30) //we start from the real solution
// 		return  BaseType::IsConverged();

	double auxNorm = 0.0;
	do
	  {

  	//START_TIMING(t);
	  	  
	    //PreconditionedMult(rA,p,q);  // q = rA * p
	    GPU_CHECK(GPU_MatrixVectorMultiply(gA, p, q));

		//STOP_TIMING(t1, t);
		
		//START_TIMING(t);

	    //alpha = roh0 / SparseSpaceType::Dot(rs,q);
	    double temp;
	    GPU_CHECK(GPU_VectorVectorMultiply(rs, q, temp));
	    alpha = roh0 / temp;
        
	    //SparseSpaceType::ScaleAndAdd(1.00, r, -alpha, q, s); // s = r - alpha * q
	    GPU_CHECK(GPU_VectorScaleAndAdd(1.00, r, -alpha, q, s));
		//STOP_TIMING(t2, t);
		
		//START_TIMING(t);

	    //PreconditionedMult(rA,s,qs);
	    // ...
	    GPU_CHECK(GPU_MatrixVectorMultiply(gA, s, qs));

		//STOP_TIMING(t1, t);

		//START_TIMING(t);

	    //omega = SparseSpaceType::Dot(qs,qs);
	    GPU_CHECK(GPU_VectorVectorMultiply(qs, qs, omega));

	    //if(omega == 0.00)
	    if(fabs(omega) <= 1.0e-30)
	      break;

	    //omega = SparseSpaceType::Dot(qs,s) / omega;
	    GPU_CHECK(GPU_VectorVectorMultiply(qs, s, temp));
	    omega = temp / omega;

	    //SparseSpaceType::ScaleAndAdd(alpha, p, 1.00, rX);
	    GPU_CHECK(GPU_VectorScaleAndAdd(alpha, p, 1.00, gX));	    
	    
	    //SparseSpaceType::ScaleAndAdd(omega, s, 1.00, rX);
	    GPU_CHECK(GPU_VectorScaleAndAdd(omega, s, 1.00, gX));
	    
	    //SparseSpaceType::ScaleAndAdd(-omega, qs, 1.00, s, r);
	    GPU_CHECK(GPU_VectorScaleAndAdd(-omega, qs, 1.00, s, r));

	    //roh1 = SparseSpaceType::Dot(r,rs);
	    GPU_CHECK(GPU_VectorVectorMultiply(r, rs, roh1));

	    //if((roh0 == 0.00) || (omega == 0.00))
	    if((fabs(roh0) <= 1.0e-30) || (fabs(omega) <= 1.0e-30))
	      break;
	    
	    beta = (roh1 * alpha) / (roh0 * omega);
	    
	    //SparseSpaceType::ScaleAndAdd(1.00, p, -omega, q);
	    GPU_CHECK(GPU_VectorScaleAndAdd(1.00, p, -omega, q));
	    
	    //SparseSpaceType::ScaleAndAdd(1.00, r, beta, q, p);
	    GPU_CHECK(GPU_VectorScaleAndAdd(1.00, r, beta, q, p));
	      
	    roh0 = roh1;
        
		//BaseType::mResidualNorm = SparseSpaceType::TwoNorm(r);
		GPU_CHECK(GPU_VectorNorm2(r, BaseType::mResidualNorm));

		//STOP_TIMING(t2, t);

		BaseType::mIterationsNumber++;


	  } while(BaseType::IterationNeeded());

	  GPU_CHECK(gX.Copy(GPU_CPU));
	  

	GPU_CHECK(gA.GPU_Free());
	GPU_CHECK(gX.GPU_Free());
	GPU_CHECK(gB.GPU_Free());
	GPU_CHECK(r.GPU_Free());
	GPU_CHECK(p.GPU_Free());
	GPU_CHECK(s.GPU_Free());
	GPU_CHECK(q.GPU_Free());
	GPU_CHECK(rs.GPU_Free());
	GPU_CHECK(qs.GPU_Free());

	  //t = t1 + t2;
	  
	  //KRATOS_WATCH(((double) t1) / ((double) t) * 100);
	  //KRATOS_WATCH(((double) t2) / ((double) t) * 100);*/
	  
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
        
    }; // Class GPUBICGSTABSolver 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  
  inline std::istream& operator >> (std::istream& IStream, 
				      GPUBICGSTABSolver& rThis)
    {
		return IStream;
    }

  /// output stream function
  inline std::ostream& operator << (std::ostream& OStream, 
				    const GPUBICGSTABSolver& rThis)
    {
      rThis.PrintInfo(OStream);
      OStream << std::endl;
      rThis.PrintData(OStream);

      return OStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_GPU_BICGSTAB_SOLVER_H_INCLUDED  defined 


