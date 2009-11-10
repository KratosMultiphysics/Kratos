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


#if !defined(KRATOS_AMG_SOLVER_H_INCLUDED )
#define  KRATOS_AMG_SOLVER_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <ctime>
using namespace std;
// External includes 


// Project includes
#include "includes/define.h"
#include "linear_solvers/iterative_solver.h"
#include "Kratos_AMGpreconditioner.h"


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
 /* template<class TSparseSpaceType, class TDenseSpaceType, 
    class TPreconditionerType = Preconditioner<TSparseSpaceType, TDenseSpaceType>, 
    class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
    class CGSolver : public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType>*/
    template<class TSparseSpaceType, class TDenseSpaceType, 
    class TPreconditionerType = Preconditioner<TSparseSpaceType, TDenseSpaceType>, 
    class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
    class AMGSolver : public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType>
    {
      public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of CGSolver
      KRATOS_CLASS_POINTER_DEFINITION(AMGSolver);

      typedef IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType> BaseType; 
  
      typedef typename TSparseSpaceType::MatrixType SparseMatrixType;
  
      typedef typename TSparseSpaceType::VectorType VectorType;
  
      typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      AMGSolver(){}

      AMGSolver(double NewMaxTolerance) : BaseType(NewMaxTolerance)
      {}

      AMGSolver(double NewMaxTolerance, unsigned int NewMaxIterationsNumber, double _W, size_t _numLevelsRoh, bool _assumeZerosForEachStep, size_t _numMaxHierarchyLevels, size_t _minimumSizeAllowed, const Kratos::Vector _preSweeps, const Kratos::Vector _postSweeps) : BaseType(NewMaxTolerance, NewMaxIterationsNumber)
      {
		tol = NewMaxTolerance;
		maxIter = NewMaxIterationsNumber;
		preconditioner = new Kratos_AMGpreconditioner(_W, _numLevelsRoh, _assumeZerosForEachStep, _numMaxHierarchyLevels, _minimumSizeAllowed, _preSweeps, _postSweeps);
      }

      /// Copy constructor.
      AMGSolver(const AMGSolver& Other) : BaseType(Other) {}


      /// Destructor.
      virtual ~AMGSolver(){  delete preconditioner ;}
      

      ///@}
      ///@name Operators 
      ///@{
      
      /// Assignment operator.
      AMGSolver& operator=(const AMGSolver& Other)
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

	  
		/** TODO here preconditioner must solve the problem **/
	//Allocating matrix A
	GPUCSRMatrix gpuA(rA.nnz(), rA.size1(), rA.size2(), &(rA.index2_data() [0]), &(rA.index1_data() [0]), &(rA.value_data() [0]));
	KRATOS_GPU_CHECK(gpuA.GPU_Allocate());
	KRATOS_GPU_CHECK(gpuA.Copy(CPU_GPU, false));

	//Allocating vector b
	GPUVector gpuB(rB.size(), &(rB[0]));
	KRATOS_GPU_CHECK(gpuB.GPU_Allocate());
	KRATOS_GPU_CHECK(gpuB.Copy(CPU_GPU));

	//Allocating vector x
	GPUVector gpuX(rX.size(), &(rX[0]));
	KRATOS_GPU_CHECK(gpuX.GPU_Allocate());
	KRATOS_GPU_CHECK(gpuX.Copy(CPU_GPU));

	this->preconditioner->initialize(gpuA.CPU_RowIndices, gpuA.CPU_Columns, gpuA.CPU_Values,
				gpuA.GPU_RowIndices, gpuA.GPU_Columns, gpuA.GPU_Values,
				gpuA.Size1, gpuA.Size2, gpuA.NNZ, true, true);

	BaseType::mIterationsNumber = this->preconditioner->solve(gpuB.GPU_Values, gpuB.CPU_Values, gpuX.GPU_Values, gpuX.CPU_Values, tol, maxIter);

	KRATOS_GPU_CHECK(gpuX.Copy(GPU_CPU));
	this->preconditioner->cleanPreconditioner();


	  return true;
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

	      //is_solved &= IterativeSolve(rA,x,b);

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

	Kratos_AMGpreconditioner *preconditioner;
	double tol;
	size_t maxIter;
        
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
				      AMGSolver<TSparseSpaceType, TDenseSpaceType, 
				      TPreconditionerType, TReordererType>& rThis)
    {
		return IStream;
    }

  /// output stream function
  template<class TSparseSpaceType, class TDenseSpaceType, 
    class TPreconditionerType, 
    class TReordererType>
  inline std::ostream& operator << (std::ostream& OStream, 
				    const AMGSolver<TSparseSpaceType, TDenseSpaceType, 
				      TPreconditionerType, TReordererType>& rThis)
    {
      rThis.PrintInfo(OStream);
      OStream << std::endl;
      rThis.PrintData(OStream);

      return OStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_AMG_SOLVER_H_INCLUDED  defined 


