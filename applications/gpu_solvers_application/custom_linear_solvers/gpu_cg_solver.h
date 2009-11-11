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


#if !defined(KRATOS_GPU_CG_SOLVER_H_INCLUDED )
#define  KRATOS_GPU_CG_SOLVER_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <ctime>
using namespace std;
// External includes 


// Project includes
#include "includes/define.h"
#include "linear_solvers/iterative_solver.h"
#include "linear_solvers.h"


//using namespace Kratos::GPUSparse;


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
 /* template<class TSparseSpaceType, class TDenseSpaceType, 
    class TPreconditionerType = Preconditioner<TSparseSpaceType, TDenseSpaceType>, 
    class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
    class CGSolver : public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType>*/
    template<class TSparseSpaceType, class TDenseSpaceType, 
    class TPreconditionerType = Preconditioner<TSparseSpaceType, TDenseSpaceType>, 
    class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
    class GPUCGSolver : public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType>
    {
      public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of CGSolver
      KRATOS_CLASS_POINTER_DEFINITION(GPUCGSolver);

      typedef IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType> BaseType; 
  
      typedef typename TSparseSpaceType::MatrixType SparseMatrixType;
  
      typedef typename TSparseSpaceType::VectorType VectorType;
  
      typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      GPUCGSolver(){}

      GPUCGSolver(double NewMaxTolerance) : BaseType(NewMaxTolerance)
	{
		havePreconditioner = false;
		maxIter = 5000;		
		tol = NewMaxTolerance;
		this->preconditioner = 0;
	}

      GPUCGSolver(double NewMaxTolerance, unsigned int NewMaxIterationsNumber) : BaseType(NewMaxTolerance, NewMaxIterationsNumber){
		havePreconditioner = false;
		maxIter = NewMaxIterationsNumber;		
		tol = NewMaxTolerance;
		this->preconditioner = 0;
      }
	//GPUCG constructor with new preconditioner
      GPUCGSolver(double NewMaxTolerance, unsigned int NewMaxIterationsNumber, GPUPreconditioner& _preconditioner) : BaseType(NewMaxTolerance, NewMaxIterationsNumber)
	{
		havePreconditioner = true;
		this->preconditioner = &_preconditioner;
		maxIter = NewMaxIterationsNumber;
		tol = NewMaxTolerance;
      }

      GPUCGSolver(double NewMaxTolerance, unsigned int NewMaxIterationsNumber, typename TPreconditionerType::Pointer pNewPreconditioner) : 
      BaseType(NewMaxTolerance, NewMaxIterationsNumber, pNewPreconditioner){}

      /// Copy constructor.
      GPUCGSolver(const GPUCGSolver& Other) : BaseType(Other) {}


      /// Destructor.
      virtual ~GPUCGSolver(){  }
      

      ///@}
      ///@name Operators 
      ///@{
      
      /// Assignment operator.
      GPUCGSolver& operator=(const GPUCGSolver& Other)
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
	  {
	    KRATOS_WATCH(rA.size1());
	    KRATOS_WATCH(rA.size2());
	    KRATOS_WATCH(rX.size());
	    KRATOS_WATCH(rB.size());
	    std::cout << "system is not consistent" << std::endl;
	    return false;
	  }
	  
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

      bool IterativeSolve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
      {

	//if(havePreconditioner)
		CG_GPU(rA.size1(), rA.size2(), rA.nnz(), &(rA.value_data() [0]), &(rA.index2_data() [0]), &(rA.index1_data() [0]), rB.size(), &(rX[0]), &(rB[0]), 
			tol, maxIter, BaseType::mBNorm, BaseType::mResidualNorm, BaseType::mIterationsNumber, *(this->preconditioner));
	//else
	//	CG_GPU(rA.size1(), rA.size2(), rA.nnz(), &(rA.value_data() [0]), &(rA.index2_data() [0]), &(rA.index1_data() [0]), rB.size(), &(rX[0]), &(rB[0]), 
	//		tol, maxIter, BaseType::mBNorm, BaseType::mResidualNorm, BaseType::mIterationsNumber, 0);
      
/*	const int size = rX.size();

	//Allocating matrix A
	GPUCSRMatrix gpuA(rA.nnz(), rA.size1(), rA.size2(), &(rA.index2_data() [0]), &(rA.index1_data() [0]), &(rA.value_data() [0]));
	GPU_CHECK(gpuA.GPU_Allocate());
	GPU_CHECK(gpuA.Copy(CPU_GPU, false));

	//Allocating vector b
	GPUVector gpuB(rB.size(), &(rB[0]));
	GPU_CHECK(gpuB.GPU_Allocate());
	GPU_CHECK(gpuB.Copy(CPU_GPU));

	//Allocating vector x
	GPUVector gpuX(rX.size(), &(rX[0]));
	GPU_CHECK(gpuX.GPU_Allocate());
	GPU_CHECK(gpuX.Copy(CPU_GPU));

	if(havePreconditioner){
		double alpha, beta, rho, rho_1 = 1.0;

		//Ini preconditioner
		clock_t s1 = clock();
		this->preconditioner->initialize(gpuA.CPU_RowIndices, gpuA.CPU_Columns, gpuA.CPU_Values,
				gpuA.GPU_RowIndices, gpuA.GPU_Columns, gpuA.GPU_Values,
				gpuA.Size1, gpuA.Size2, gpuA.NNZ, true, true);
		clock_t s2 = clock();
		std::cout << "Time to create hierarchy" << double(s2-s1) / CLOCKS_PER_SEC << "s" << std::endl;


		//Norm(b)
		GPU_CHECK(GPU_VectorNorm2(gpuB, BaseType::mBNorm));

		//r = b - A*x
		GPUVector gpuR(size);
		GPU_CHECK(gpuR.GPU_Allocate());
		GPU_CHECK(GPU_MatrixVectorMultiply(gpuA, gpuX, gpuR));
		GPU_CHECK(GPU_VectorScaleAndAdd(1.00, gpuB, -1.00, gpuR));

		//std::cout << "BNORM " << BaseType::mBNorm << std::endl;
		if(BaseType::mBNorm == 0.0)
			BaseType::mBNorm = 1.0;

		//Norm(r)
		double normr, resid;
		GPU_CHECK(GPU_VectorNorm2(gpuR, normr));
		if((resid = normr/BaseType::mBNorm) <= 1.0e-30){
			std::cout << "Out from outer return" << std::endl;
			return false;
		}
		GPUVector gpuP(size), gpuZ(size), gpuQ(size);
		GPU_CHECK(gpuP.GPU_Allocate());
		GPU_CHECK(gpuZ.GPU_Allocate());
		GPU_CHECK(gpuQ.GPU_Allocate());
		
		BaseType::mIterationsNumber = 0;

		clock_t s3 = 0.0;
		size_t i = 1;
		do{

			GPU_fillWithZeros(size, gpuZ.GPU_Values);
			s1 = clock();
			this->preconditioner->singleStep(gpuR.GPU_Values, gpuZ.GPU_Values);
			s2 = clock();
			s3 += s2-s1;
			rho = GPU_dotProduct(gpuR.Size, gpuR.GPU_Values, 1, gpuZ.GPU_Values, 1);

			if(i == 1)
				//gpuP = gpuZ
				gpuP.CopyGPUValuesFrom(gpuZ);
			else{
				beta = rho / rho_1; 
				GPU_CHECK(GPU_VectorScaleAndAdd(1.00, gpuZ, beta, gpuP));
			}

			GPU_CHECK(GPU_MatrixVectorMultiply(gpuA, gpuP, gpuQ));
			alpha = rho / GPU_dotProduct(gpuP.Size, gpuP.GPU_Values, 1, gpuQ.GPU_Values, 1);
			GPU_CHECK(GPU_VectorScaleAndAdd(alpha, gpuP, 1.00, gpuX));
			GPU_CHECK(GPU_VectorScaleAndAdd(-alpha, gpuQ, 1.00, gpuR));
		
			GPU_CHECK(GPU_VectorNorm2(gpuR, normr));
			BaseType::mResidualNorm = normr;
			if((resid = normr/BaseType::mBNorm) <= 1.0e-30){
				break;
			}
			rho_1 = rho;
			BaseType::mIterationsNumber++;			
			i++;
			KRATOS_WATCH("end iteration");
			KRATOS_WATCH(BaseType::mResidualNorm);
			KRATOS_WATCH(BaseType::mBNorm);
			KRATOS_WATCH(BaseType::mResidualNorm/BaseType::mBNorm);
		}while(BaseType::IterationNeeded());
		std::cout << "Average time for single step" << (double(s3)/(i-1)) / CLOCKS_PER_SEC << "s" << std::endl;

		GPU_CHECK(gpuX.Copy(GPU_CPU));
		this->preconditioner->cleanPreconditioner();
	}else{     		
		double alpha, beta, rho, rho_1 = 1.0;

		//Norm(b)
		GPU_CHECK(GPU_VectorNorm2(gpuB, BaseType::mBNorm));

		//r = b - A*x
		GPUVector gpuR(size);
		GPU_CHECK(gpuR.GPU_Allocate());
		GPU_CHECK(GPU_MatrixVectorMultiply(gpuA, gpuX, gpuR));
		GPU_CHECK(GPU_VectorScaleAndAdd(1.00, gpuB, -1.00, gpuR));

		//std::cout << "BNORM " << BaseType::mBNorm << std::endl
		if(BaseType::mBNorm == 0.0)
			BaseType::mBNorm = 1.0;

		//Norm(r)
		double normr, resid;
		GPU_CHECK(GPU_VectorNorm2(gpuR, normr));
		if((resid = normr/BaseType::mBNorm) <= 1.0e-30){
			std::cout << "Out from outer return" << std::endl;
			return false;
		}
		GPUVector gpuP(size), gpuQ(size);
		GPU_CHECK(gpuP.GPU_Allocate());
		GPU_CHECK(gpuQ.GPU_Allocate());
		
		BaseType::mIterationsNumber = 0;

		size_t i = 1;
		do{

			rho = GPU_dotProduct(gpuR.Size, gpuR.GPU_Values, 1, gpuR.GPU_Values, 1);

			if(i == 1)
				//gpuP = gpuZ
				gpuP.CopyGPUValuesFrom(gpuR);
			else{
				beta = rho / rho_1; 
				GPU_CHECK(GPU_VectorScaleAndAdd(1.00, gpuR, beta, gpuP));
			}

			GPU_CHECK(GPU_MatrixVectorMultiply(gpuA, gpuP, gpuQ));
			alpha = rho / GPU_dotProduct(gpuP.Size, gpuP.GPU_Values, 1, gpuQ.GPU_Values, 1);
			GPU_CHECK(GPU_VectorScaleAndAdd(alpha, gpuP, 1.00, gpuX));
			GPU_CHECK(GPU_VectorScaleAndAdd(-alpha, gpuQ, 1.00, gpuR));
		
			GPU_CHECK(GPU_VectorNorm2(gpuR, normr));
			BaseType::mResidualNorm = normr;
			if((resid = normr/BaseType::mBNorm) <= 1.0e-30){
				break;
			}
			rho_1 = rho;
			BaseType::mIterationsNumber++;			
			i++;
			KRATOS_WATCH("end iteration");
			KRATOS_WATCH(BaseType::mResidualNorm);
			KRATOS_WATCH(BaseType::mBNorm);
			KRATOS_WATCH(BaseType::mResidualNorm/BaseType::mBNorm);
		}while(BaseType::IterationNeeded());

		GPU_CHECK(gpuX.Copy(GPU_CPU));

	}*/
	//std::cout << "Got return with value of convergence " << BaseType::IsConverged() << ", and finalIterations " << BaseType::mIterationsNumber << std::endl;	return BaseType::IsConverged();
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
				      GPUCGSolver<TSparseSpaceType, TDenseSpaceType, 
				      TPreconditionerType, TReordererType>& rThis)
    {
		return IStream;
    }

  /// output stream function
  template<class TSparseSpaceType, class TDenseSpaceType, 
    class TPreconditionerType, 
    class TReordererType>
  inline std::ostream& operator << (std::ostream& OStream, 
				    const GPUCGSolver<TSparseSpaceType, TDenseSpaceType, 
				      TPreconditionerType, TReordererType>& rThis)
    {
      rThis.PrintInfo(OStream);
      OStream << std::endl;
      rThis.PrintData(OStream);

      return OStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_GPU_CG_SOLVER_H_INCLUDED  defined 


