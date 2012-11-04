#if !defined(KRATOS_Pastix_Iterative_Solver )
#define  KRATOS_Pastix_Iterative_Solver

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
	int solvePASTIX(int echo_level,int mat_size, int nnz, double* AA, size_t* IA, size_t* JA, double* x, double* b);
}


template< class TSparseSpaceType, class TDenseSpaceType,
class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class Pastix_Iterative_Solver : public DirectSolver< TSparseSpaceType,
	TDenseSpaceType, TReordererType>
{
public:
	/**
	 * Counted pointer of Pastix_Iterative_Solver
	 */
	typedef boost::shared_ptr<Pastix_Iterative_Solver> Pointer;

	typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

	typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

	typedef typename TSparseSpaceType::VectorType VectorType;

	typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

	/**
	 * Default constructor
	 */
	Pastix_Iterative_Solver(double NewMaxTolerance,
	                        int NewMaxIterationsNumber,
	                        int restart)
	{
		mTol = NewMaxTolerance;
		mmax_it = NewMaxIterationsNumber;
		mrestart = restart;
	}

	Pastix_Iterative_Solver()
	{
		mTol = 1e-6;
		mrestart = 150;
		mmax_it = mrestart*3;


	}

	/**
	 * Destructor
	 */
	virtual ~Pastix_Iterative_Solver() {};

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
//		int *index1_vector = new(std::nothrow) int[rA.index1_data().size()];
//		int *index2_vector = new(std::nothrow) int[rA.index2_data().size()];
////
//		for(int unsigned i = 0; i < rA.index1_data().size(); i++)
//			index1_vector[i] = (int)rA.index1_data()[i];
//
//		for(unsigned int i = 0; i < rA.index2_data().size(); i++)
//			index2_vector[i] = (int)rA.index2_data()[i];

		int echo_level = 1; //does not work if we set it to 0 ... should investigate further!

		int state = solvePASTIX(echo_level, rA.size1(), rA.value_data().size(), rA.value_data().begin(), &(rA.index1_data()[0]), &(rA.index2_data()[0]), &rX[0], &rB[0]);
//		int state = solvePASTIX(echo_level, rA.size1(), rA.value_data().size(), rA.value_data().begin(), index1_vector, index2_vector, &rX[0], &rB[0]);

//		delete [] index1_vector;
//		delete [] index2_vector;
		
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
		rOStream << "SuperLU solver finished.";
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
	int mrestart;
//     double mDropTol;
//     double mFillTol;
//     double mFillFactor;

	/**
	 * Assignment operator.
	 */
	Pastix_Iterative_Solver& operator=(const Pastix_Iterative_Solver& Other);

	/**
	 * Copy constructor.
	 */
	Pastix_Iterative_Solver(const Pastix_Iterative_Solver& Other);

}; // Class SkylineLUFactorizationSolver


/**
 * input stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType,class TReordererType>
inline std::istream& operator >> (std::istream& rIStream, Pastix_Iterative_Solver< TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
	return rIStream;
}

/**
 * output stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const Pastix_Iterative_Solver<TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
	rThis.PrintInfo(rOStream);
	rOStream << std::endl;
	rThis.PrintData(rOStream);

	return rOStream;
}

//#undef MPI_COMM_WORLD

}  // namespace Kratos.


#endif // KRATOS_Pastix_Iterative_Solver  defined
