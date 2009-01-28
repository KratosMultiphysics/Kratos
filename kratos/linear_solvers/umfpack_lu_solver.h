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
 


#if !defined(KRATOS_UMFPACK_LU_SOLVER_H_INCLUDED )
#define  KRATOS_UMFPACK_LU_SOLVER_H_INCLUDED



// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "linear_solvers/direct_solver.h"
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/ublas_sparse.hpp>
#include <boost/numeric/bindings/umfpack/umfpack.hpp>

namespace umf = boost::numeric::bindings::umfpack;

namespace Kratos
{


	template<class TSparseSpaceType, class TDenseSpaceType, 
	class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
	class UMFpackLUsolver 
		: public DirectSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>
	{
	public:

		/// Counted pointer of UMFpackLUsolver
		typedef boost::shared_ptr<UMFpackLUsolver> Pointer;

		typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType; 

		typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

		typedef typename TSparseSpaceType::VectorType VectorType;

		typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

		/// Default constructor.
		UMFpackLUsolver(){}

		/// Destructor.
		virtual ~UMFpackLUsolver(){}


		/** Normal solve method.
		Solves the linear system Ax=b and puts the result on SystemVector& rX. 
		rX is also th initial guess for iterative methods.
		@param rA. System matrix
		@param rX. Solution vector.
		@param rB. Right hand side vector.
		*/
		bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
		{
			KRATOS_TRY
			if(IsNotConsistent(rA, rX, rB))
				return false;

			umf::symbolic_type< TSparseSpaceType::DataType > Symbolic;
			umf::numeric_type< TSparseSpaceType::DataType > Numeric;

			umf::symbolic(A, Symbolic);
			umf::numeric(A, Symbolic, Numeric);
			umf::symbolic(A,x,b,Numeric);

			return true;

			KRATOS_CATCH("");
		}

		/** Multi solve method for solving a set of linear systems with same coefficient matrix.
		Solves the linear system Ax=b and puts the result on SystemVector& rX. 
		rX is also th initial guess for iterative methods.
		@param rA. System matrix
		@param rX. Solution vector.
		@param rB. Right hand side vector.
		*/
		//bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB)
		//{


		//	return is_solved;
		//}



		/// Print information about this object.
		void  PrintInfo(std::ostream& rOStream) const
		{
			rOStream << "LU factorization solver finished.";
		}

		/// Print object's data.
		void  PrintData(std::ostream& rOStream) const 
		{
		}




	private:



		/// Assignment operator.
		UMFpackLUsolver& operator=(const UMFpackLUsolver& Other);

		/// Copy constructor.
		UMFpackLUsolver(const UMFpackLUsolver& Other);


	}; // Class UMFpackLUsolver 


	/// input stream function
	template<class TSparseSpaceType, class TDenseSpaceType,class TReordererType>
		inline std::istream& operator >> (std::istream& rIStream, UMFpackLUsolver<TSparseSpaceType, 
		TDenseSpaceType, 
		TReordererType>& rThis)
	{
	}

	/// output stream function
	template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
		inline std::ostream& operator << (std::ostream& rOStream, 
		const UMFpackLUsolver<TSparseSpaceType, 
		TDenseSpaceType, 
		TReordererType>& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}


}  // namespace Kratos.

#endif // KRATOS_UMFPACK_LU_SOLVER_H_INCLUDED  defined 


