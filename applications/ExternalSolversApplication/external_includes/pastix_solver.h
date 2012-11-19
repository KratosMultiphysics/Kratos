#if !defined(KRATOS_PastixSolver )
#define  KRATOS_PastixSolver

// External includes
#include "boost/smart_ptr.hpp"

#include "includes/ublas_interface.h"
// #include "external_includes/superlu/superlu.hpp"

//#include "solveARMS.h"

// Project includes
#include "includes/define.h"
#include "linear_solvers/direct_solver.h"

//#include <complex.h>
/* to access functions from the libpastix, respect this order */
//#include "custom_external_libraries/solvePASTIX.c"

namespace ublas = boost::numeric::ublas;

namespace Kratos
{

extern "C"
{
	int solvePASTIX(int verbosity,int mat_size, int nnz, double* AA, size_t* IA, size_t* JA, double* x, double* b, int m_gmres, 
				double tol, int incomplete, int ilu_level_of_fill);
}


template< class TSparseSpaceType, class TDenseSpaceType,
class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class PastixSolver : public DirectSolver< TSparseSpaceType,
	TDenseSpaceType, TReordererType>
{
public:
	/**
	 * Counted pointer of PastixSolver
	 */
	typedef boost::shared_ptr<PastixSolver> Pointer;

	typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

	typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

	typedef typename TSparseSpaceType::VectorType VectorType;

	typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

	/**
	 * Default constructor - uses ILU+GMRES
	 * @param NewMaxTolerance tolerance that will be achieved by the iterative solver
	 * @param NewMaxIterationsNumber this number represents both the number of iterations AND the size of the krylov space
	 * @param level of fill that will be used in the ILU
	 * @param verbosity, a number from 0 (no output) to 2 (maximal output)
	 */
	PastixSolver(double NewMaxTolerance,
	                        int NewMaxIterationsNumber,
	                        int level_of_fill,
							int verbosity)
	{
		std::cout << "setting up pastix for iterative solve " << std::endl;
		mTol = NewMaxTolerance;
		mmax_it = NewMaxIterationsNumber;
		mlevel_of_fill = level_of_fill;
		mincomplete = 1;
		mverbosity=verbosity;
	}

	/**
	 * Direct Solver
	 * @param verbosity, a number from 0 (no output) to 2 (maximal output)
	 */
	PastixSolver(int verbosity)
	{
		std::cout << "setting up pastix for direct solve " <<std::endl;
		mTol = -1;
		mmax_it = -1;
		mlevel_of_fill = -1;
		mincomplete = 0;
		mverbosity=verbosity;


	}

	/**
	 * Destructor
	 */
	virtual ~PastixSolver() {};

	/**
	 * Normal solve method.
	 * Solves the linear system Ax=b and puts the result on SystemVector& rX.
	 * rX is also th initial guess for iterative methods.
	 * @param rA. System matrix
	 * @param rX. Solution vector.
	 * @param rB. Right hand side vector.
	 */
	bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
	{
		KRATOS_WATCH(__LINE__);
				
		int state = solvePASTIX(mverbosity, rA.size1(), rA.value_data().size(), rA.value_data().begin(), &(rA.index1_data()[0]), &(rA.index2_data()[0]), &rX[0], &rB[0]
		,mmax_it,mTol,mincomplete,mlevel_of_fill);

		return state;
	}

	/**
	 * Multi solve method for solving a set of linear systems with same coefficient matrix.
	 * Solves the linear system Ax=b and puts the result on SystemVector& rX.
	 * rX is also th initial guess for iterative methods.
	 * @param rA. System matrix
	 * @param rX. Solution vector.
	 * @param rB. Right hand side vector.
	 */
	bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB)
	{

		bool is_solved = true;


		return is_solved;
	}

	/**
	 * Print information about this object.
	 */
	void  PrintInfo(std::ostream& rOStream) const
	{
		rOStream << "Pastix solver finished.";
	}

	/**
	 * Print object's data.
	 */
	void  PrintData(std::ostream& rOStream) const
	{
	}

private:

	double mTol;
	int mmax_it;
	int mincomplete;
	int mlevel_of_fill;
	int mverbosity;
//     double mDropTol;
//     double mFillTol;
//     double mFillFactor;

	/**
	 * Assignment operator.
	 */
	PastixSolver& operator=(const PastixSolver& Other);

	/**
	 * Copy constructor.
	 */
	PastixSolver(const PastixSolver& Other);

}; // Class PastixSolver


/**
 * input stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType,class TReordererType>
inline std::istream& operator >> (std::istream& rIStream, PastixSolver< TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
	return rIStream;
}

/**
 * output stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const PastixSolver<TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
	rThis.PrintInfo(rOStream);
	rOStream << std::endl;
	rThis.PrintData(rOStream);

	return rOStream;
}

//#undef MPI_COMM_WORLD

}  // namespace Kratos.


#endif // KRATOS_PastixSolver  defined
