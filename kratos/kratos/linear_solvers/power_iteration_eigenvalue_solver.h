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
//   Date:                $Date: 2009-01-15 11:11:35 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_POWER_ITERATION_EIGENVALUE_SOLVERR_H_INCLUDED )
#define  KRATOS_POWER_ITERATION_EIGENVALUE_SOLVERR_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <numeric>
#include <vector>


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
    template<class TSparseSpaceType, class TDenseSpaceType, class TLinearSolverType,
    class TPreconditionerType = Preconditioner<TSparseSpaceType, TDenseSpaceType>, 
    class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
    class PowerIterationEigenvalueSolver : public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType>
    {
      public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of PowerIterationEigenvalueSolver
      KRATOS_CLASS_POINTER_DEFINITION(PowerIterationEigenvalueSolver);

      typedef IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType> BaseType; 
  
      typedef typename TSparseSpaceType::MatrixType SparseMatrixType;
  
      typedef typename TSparseSpaceType::VectorType VectorType;
  
      typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
  
      typedef typename TDenseSpaceType::VectorType DenseVectorType;

      typedef std::size_t SizeType;

      typedef std::size_t IndexType;
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      PowerIterationEigenvalueSolver(){}

       PowerIterationEigenvalueSolver(double NewMaxTolerance, unsigned int NewMaxIterationsNumber, 
			       unsigned int NewRequiredEigenvalueNumber, typename TLinearSolverType::Pointer pLinearSolver) 
      : BaseType(NewMaxTolerance, NewMaxIterationsNumber), mRequiredEigenvalueNumber(NewRequiredEigenvalueNumber), mpLinearSolver(pLinearSolver){}

/*       PowerIterationEigenvalueSolver(double NewMaxTolerance, unsigned int NewMaxIterationsNumber, typename TPreconditionerType::Pointer pNewPreconditioner) :  */
/*       BaseType(NewMaxTolerance, NewMaxIterationsNumber, pNewPreconditioner){} */

      /// Copy constructor.
      PowerIterationEigenvalueSolver(const PowerIterationEigenvalueSolver& Other) : BaseType(Other) {}


      /// Destructor.
      virtual ~PowerIterationEigenvalueSolver(){}
      

      ///@}
      ///@name Operators 
      ///@{
      
      /// Assignment operator.
      PowerIterationEigenvalueSolver& operator=(const PowerIterationEigenvalueSolver& Other)
      {
        BaseType::operator=(Other);
	return *this;
      }
      
      ///@}
      ///@name Operations
      ///@{
      
      static void RandomInitialize(DenseVectorType& R)
      {
	  for(SizeType i = 0 ; i < R.size() ; i++)
	      R[i] = 1.00; //rand();

	  R /= norm_2(R);
      }


      // The power iteration algorithm 
	  void Solve(SparseMatrixType& K, 
		 SparseMatrixType& M,
		 DenseVectorType& Eigenvalues,
		 DenseMatrixType& Eigenvectors)
      {

		using boost::numeric::ublas::trans;
	  
		SizeType size = K.size1();
		SizeType max_iteration = BaseType::GetMaxIterationsNumber();
		double tolerance = BaseType::GetTolerance();

		VectorType x = ZeroVector(size);
		VectorType y = ZeroVector(size);

		RandomInitialize(y);

		if(Eigenvalues.size() < 1)
			Eigenvalues.resize(1,0.00);


		// Starting with first step
		double beta = 0.00;
		double ro = 0.00;
		double old_ro = Eigenvalues[0];
		std::cout << "iteration    beta \t\t ro \t\t convergence norm" << std::endl;
		for(SizeType i = 0 ; i < max_iteration ; i++)
		{
			//K*x = y
			mpLinearSolver->Solve(K,x,y);

			ro = inner_prod(y,x);

			//y = M*x
			noalias(y) = prod(M,x);

			beta = inner_prod(x, y);
			if(beta <= 0.00)
				KRATOS_ERROR(std::invalid_argument, "M is not Positive-definite", "");

			ro = ro / beta;
			beta = sqrt(beta);

			double inverse_of_beta = 1.00 / beta;

			y *= inverse_of_beta;

			if(ro == 0.00)
				KRATOS_ERROR(std::runtime_error, "Perpendicular eigenvector to M", "");

			double convergence_norm = fabs((ro - old_ro) / ro);
			
			std::cout << i << " \t " << beta << " \t " << ro << " \t " << convergence_norm << std::endl;
			//std::cout << "i = " << i << ": beta = " << beta << ", ro = " << ro << ", convergence norm = " << convergence_norm << std::endl;
			
			if(convergence_norm < tolerance)
				break;

			old_ro = ro;



		}

KRATOS_WATCH(ro);
//KRATOS_WATCH(y);

		Eigenvalues[0] = ro;

		if((Eigenvectors.size1() < 1) || (Eigenvectors.size2() < size))
			Eigenvectors.resize(1,size);

		for(SizeType i = 0 ; i < size ; i++)
			Eigenvectors(0,i) = y[i];
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
	  buffer << "Power iteration eigenvalue solver with " << BaseType::GetPreconditioner()->Info();
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


      unsigned int mRequiredEigenvalueNumber;

      typename TLinearSolverType::Pointer mpLinearSolver;
        
      std::vector<DenseVectorType> mQVector;
      std::vector<DenseVectorType> mPVector;
      std::vector<DenseVectorType> mRVector;
        
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
        
    }; // Class PowerIterationEigenvalueSolver 

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
				      PowerIterationEigenvalueSolver<TSparseSpaceType, TDenseSpaceType, 
				      TPreconditionerType, TReordererType>& rThis)
    {
		return IStream;
    }

  /// output stream function
  template<class TSparseSpaceType, class TDenseSpaceType, 
    class TPreconditionerType, 
    class TReordererType>
  inline std::ostream& operator << (std::ostream& OStream, 
				    const PowerIterationEigenvalueSolver<TSparseSpaceType, TDenseSpaceType, 
				      TPreconditionerType, TReordererType>& rThis)
    {
      rThis.PrintInfo(OStream);
      OStream << std::endl;
      rThis.PrintData(OStream);

      return OStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_POWER_ITERATION_EIGENVALUE_SOLVERR_H_INCLUDED defined 































