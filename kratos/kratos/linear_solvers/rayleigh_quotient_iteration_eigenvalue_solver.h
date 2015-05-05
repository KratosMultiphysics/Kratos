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
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2008-03-25 15:55:47 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_RAYLEIGH_QUOTIENT_ITERATION_EIGENVALUE_SOLVER_H_INCLUDED )
#define  KRATOS_RAYLEIGH_QUOTIENT_ITERATION_EIGENVALUE_SOLVER_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <numeric>
#include <vector>


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
    template<class TSparseSpaceType, class TDenseSpaceType, class TLinearSolverType,
    class TPreconditionerType = Preconditioner<TSparseSpaceType, TDenseSpaceType>, 
    class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
    class RayleighQuotientIterationEigenvalueSolver : public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType>
    {
      public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of RayleighQuotientIterationEigenvalueSolver
      KRATOS_CLASS_POINTER_DEFINITION(RayleighQuotientIterationEigenvalueSolver);

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
      RayleighQuotientIterationEigenvalueSolver(){}

       RayleighQuotientIterationEigenvalueSolver(double NewMaxTolerance, unsigned int NewMaxIterationsNumber, 
			       unsigned int NewRequiredEigenvalueNumber, typename TLinearSolverType::Pointer pLinearSolver, double ShiftingConvergence = 0.25) 
      : BaseType(NewMaxTolerance, NewMaxIterationsNumber), mRequiredEigenvalueNumber(NewRequiredEigenvalueNumber), mpLinearSolver(pLinearSolver), mShiftingConvergence(ShiftingConvergence){}

/*       RayleighQuotientIterationEigenvalueSolver(double NewMaxTolerance, unsigned int NewMaxIterationsNumber, typename TPreconditionerType::Pointer pNewPreconditioner) :  */
/*       BaseType(NewMaxTolerance, NewMaxIterationsNumber, pNewPreconditioner){} */

      /// Copy constructor.
      RayleighQuotientIterationEigenvalueSolver(const RayleighQuotientIterationEigenvalueSolver& Other) : BaseType(Other)
		  , mRequiredEigenvalueNumber(Other.mRequiredEigenvalueNumber), mpLinearSolver(Other.mpLinearSolver), mShiftingConvergence(Other.mShiftingConvergence) {}


      /// Destructor.
      virtual ~RayleighQuotientIterationEigenvalueSolver(){}
      

      ///@}
      ///@name Operators 
      ///@{
      
      /// Assignment operator.
      RayleighQuotientIterationEigenvalueSolver& operator=(const RayleighQuotientIterationEigenvalueSolver& Other)
      {
        BaseType::operator=(Other);
	return *this;
      }
      
      ///@}
      ///@name Operations
      ///@{
      
      static void Initialize(DenseVectorType& R, 
		 SparseMatrixType& M)
      {
	  for(SizeType i = 0 ; i < R.size() ; i++)
	      R[i] = M(i,i);

	  if(norm_2(R) == 0.00)
		  KRATOS_THROW_ERROR(std::invalid_argument, "Invalid M matrix. The norm2 of its diagonal is Zero", "");

	  R /= norm_2(R);
      }

	  SizeType SturmSequenceCheck(SparseMatrixType& ShiftedK)
	  {
            // define an object to store skyline matrix and factorization
            LUSkylineFactorization<TSparseSpaceType, TDenseSpaceType> myFactorization;

			// copy myMatrix into skyline format
            myFactorization.copyFromCSRMatrix(ShiftedK);

			// factorize it
            myFactorization.factorize();
			SizeType counter = 0; // number of eigenvalues less than the shift
			for(SizeType i = 0 ; i < ShiftedK.size1() ; i++)
				if(myFactorization.entriesD[i] < 0.00)
					counter++;

			return counter;

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

		Initialize(y,M);

		if(Eigenvalues.size() < 1)
			Eigenvalues.resize(1,0.00);


		const double epsilon = 1.00e-9;
		// Starting with first step
		double beta = 0.00;
		double ro = 0.00;
		double shift_value = 0.00;
		double old_ro = 0.00;//Eigenvalues[0];
		std::cout << "iteration    beta \t ro \t\t convergence norm \t min \t\t max" << std::endl;

		SparseMatrixType shifted_k(K);
		// define an object to store skyline matrix and factorization
		LUSkylineFactorization<TSparseSpaceType, TDenseSpaceType> my_factorization;

		// copy myMatrix into skyline format
		my_factorization.copyFromCSRMatrix(shifted_k);

		// factorize it
		my_factorization.factorize();

		double min_shift_value = 0.00;
		double max_shift_value = 0.00;

		SizeType smaller_eigenvalue_numbers = 0;

		for(SizeType i = 0 ; i < max_iteration ; i++)
		{
			//K*x = y
			//mpLinearSolver->Solve(shifted_k,x,y);
			my_factorization.backForwardSolve(size, y, x);

			ro = inner_prod(y,x);

			//y = M*x
			noalias(y) = prod(M,x);

			beta = inner_prod(x, y);
			if(beta == 0.00)
				KRATOS_THROW_ERROR(std::invalid_argument, "Zero beta norm!", "");

			double delta_ro = (ro / beta);

			ro = delta_ro + shift_value;

			//if(ro < 0.00)
			//	ro = -ro;

			if(ro == 0.00)
				KRATOS_THROW_ERROR(std::runtime_error, "Perpendicular eigenvector to M", "");


			double convergence_norm = fabs((ro - old_ro) / ro);

			if(convergence_norm < mShiftingConvergence) // Start shifting after certain convergence
			{

				// if there are no smaller eigenvalues yet we need to extend the range
				if(smaller_eigenvalue_numbers == 0)
					max_shift_value = ro;

				if((ro > max_shift_value)||(ro < min_shift_value))
					shift_value = (max_shift_value + min_shift_value) / 2.00;
				else
                shift_value = ro;


				noalias(shifted_k) = K - shift_value*M;
				// copy myMatrix into skyline format
				my_factorization.copyFromCSRMatrix(shifted_k);

				// factorize it
				my_factorization.factorize();
					SizeType new_smaller_eigenvalue_numbers = SturmSequenceCheck(shifted_k);

					if(new_smaller_eigenvalue_numbers == 0)
					{
						min_shift_value = shift_value;
            }
					else
					{
						max_shift_value = shift_value;
						smaller_eigenvalue_numbers = new_smaller_eigenvalue_numbers;
					}

				
					unsigned int iteration_number = 0;
					unsigned int max_shift_number = 4;
				while((smaller_eigenvalue_numbers > 1) && (max_shift_value-min_shift_value > epsilon) && (iteration_number++ < max_shift_number))
				{
					shift_value = (max_shift_value + min_shift_value) / 2.00;
					noalias(shifted_k) = K - shift_value*M;
					// copy myMatrix into skyline format
					my_factorization.copyFromCSRMatrix(shifted_k);

					// factorize it
					my_factorization.factorize();

					new_smaller_eigenvalue_numbers = SturmSequenceCheck(shifted_k);

					if(new_smaller_eigenvalue_numbers == 0)
					{
						min_shift_value = shift_value;
						std::cout << "			Finding " << smaller_eigenvalue_numbers << " eigenvalues in [" << min_shift_value << "," << max_shift_value  << "]" << std::endl;
					}
					else
					{
						max_shift_value = shift_value;
						smaller_eigenvalue_numbers = new_smaller_eigenvalue_numbers;
						std::cout << "			Finding " << smaller_eigenvalue_numbers << " eigenvalues in [" << min_shift_value << "," << max_shift_value  << "]" << std::endl;
					}
				}

			}


			if(beta < 0.00)
				beta = -sqrt(-beta);
			else
				//KRATOS_THROW_ERROR(std::invalid_argument, "M is not Positive-definite", "");
			beta = sqrt(beta);

			double inverse_of_beta = 1.00 / beta;

			y *= inverse_of_beta;
			
			std::cout << i << " \t " << beta << " \t " << ro << " \t " << convergence_norm  << " \t\t " <<  min_shift_value << " \t " << max_shift_value << std::endl;
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

		//double y_norm = TSparseSpaceType::TwoNorm(y);

		for(SizeType i = 0 ; i < size ; i++)
			Eigenvectors(0,i) = x[i] / beta;
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

	  double mShiftingConvergence;

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
        
    }; // Class RayleighQuotientIterationEigenvalueSolver 

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
				      RayleighQuotientIterationEigenvalueSolver<TSparseSpaceType, TDenseSpaceType, 
				      TPreconditionerType, TReordererType>& rThis)
    {
		return IStream;
    }

  /// output stream function
  template<class TSparseSpaceType, class TDenseSpaceType, 
    class TPreconditionerType, 
    class TReordererType>
  inline std::ostream& operator << (std::ostream& OStream, 
				    const RayleighQuotientIterationEigenvalueSolver<TSparseSpaceType, TDenseSpaceType, 
				      TPreconditionerType, TReordererType>& rThis)
    {
      rThis.PrintInfo(OStream);
      OStream << std::endl;
      rThis.PrintData(OStream);

      return OStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_RAYLEIGH_QUOTIENT_ITERATION_EIGENVALUE_SOLVER_H_INCLUDED defined 































