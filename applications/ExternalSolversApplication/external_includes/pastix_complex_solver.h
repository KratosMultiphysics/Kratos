#if !defined(KRATOS_PASTIX_COMPLEX_SOLVER)
#define KRATOS_PASTIX_COMPLEX_SOLVER

#ifdef _OPENMP
#  include <omp.h>
#endif
#include <complex>
#include <vector>
#include <algorithm>
extern "C" {
#include <pastix.h>
}

namespace Kratos
{
template <class TSparseSpaceType, class TDenseSpaceType, class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType>>
class PastixComplexSolver : public DirectSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(PastixComplexSolver);

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    PastixComplexSolver(Parameters& rSettings)
    {
    }

    virtual ~PastixComplexSolver()
	{
	}

	void InitializeSolutionStep(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
		
    }

    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
		pastix_data_t *pastix_data = NULL;
		pastix_int_t nrows = rA.size1();
		std::vector<pastix_int_t> row_ptr(rA.index1_data().size());
		std::vector<pastix_int_t> col_idx(rA.index2_data().size());
		std::vector<pastix_int_t> perm(nrows), invp(nrows);
		pastix_int_t iparm[IPARM_SIZE];
		double dparm[DPARM_SIZE];

		std::copy(std::begin(rA.index1_data()),std::end(rA.index1_data()),std::begin(row_ptr));
		for (auto &i : row_ptr)
			++i;
		std::copy(std::begin(rA.index2_data()),std::end(rA.index2_data()),std::begin(col_idx));
		for (auto &i : col_idx)
			++i;

		// Initialize
		iparm[IPARM_MODIFY_PARAMETER] = API_NO;
		//iparm[IPARM_MATRIX_VERIFICATION] = API_YES;
		iparm[IPARM_START_TASK] = API_TASK_INIT;
        iparm[IPARM_END_TASK  ] = API_TASK_INIT;

		z_pastix(&pastix_data,0,nrows,
			static_cast<pastix_int_t*>(&row_ptr[0]),
			static_cast<pastix_int_t*>(&col_idx[0]),
			static_cast<std::complex<double>*>(&rA.value_data()[0]),
			static_cast<pastix_int_t*>(&perm[0]),
			static_cast<pastix_int_t*>(&invp[0]),
			NULL,
			1,iparm,dparm);

		// Factorize
	    iparm[IPARM_VERBOSE        ] = API_VERBOSE_NOT;
      	iparm[IPARM_RHS_MAKING     ] = API_RHS_B;
        iparm[IPARM_SYM            ] = API_SYM_NO;
        iparm[IPARM_FACTORIZATION  ] = API_FACT_LU;
        iparm[IPARM_TRANSPOSE_SOLVE] = API_YES;
#ifdef _OPENMP
        iparm[IPARM_THREAD_NBR     ] = omp_get_max_threads();
#endif

        iparm[IPARM_START_TASK] = API_TASK_ORDERING;
        iparm[IPARM_END_TASK  ] = API_TASK_NUMFACT;

		z_pastix(&pastix_data,0,nrows,
			static_cast<pastix_int_t*>(&row_ptr[0]),
			static_cast<pastix_int_t*>(&col_idx[0]),
			static_cast<std::complex<double>*>(&rA.value_data()[0]),
			static_cast<pastix_int_t*>(&perm[0]),
			static_cast<pastix_int_t*>(&invp[0]),
			NULL,
			1,iparm,dparm);

		// Solve
		noalias(rX) = rB;
        iparm[IPARM_START_TASK] = API_TASK_SOLVE;
        iparm[IPARM_END_TASK  ] = API_TASK_SOLVE;

		z_pastix(&pastix_data,0,nrows,
			static_cast<pastix_int_t*>(&row_ptr[0]),
			static_cast<pastix_int_t*>(&col_idx[0]),
			static_cast<std::complex<double>*>(&rA.value_data()[0]),
			static_cast<pastix_int_t*>(&perm[0]),
			static_cast<pastix_int_t*>(&invp[0]),
			static_cast<std::complex<double>*>(&rX[0]),
			1,iparm,dparm);

		// Clean
        iparm[IPARM_START_TASK] = API_TASK_CLEAN;
        iparm[IPARM_END_TASK  ] = API_TASK_CLEAN;

		z_pastix(&pastix_data,0,nrows,
			static_cast<pastix_int_t*>(&row_ptr[0]),
			static_cast<pastix_int_t*>(&col_idx[0]),
			static_cast<std::complex<double>*>(&rA.value_data()[0]),
			static_cast<pastix_int_t*>(&perm[0]),
			static_cast<pastix_int_t*>(&invp[0]),
			NULL,
			1,iparm,dparm);

		return true;
    }

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Pastix direct solver finished.";
    }

    void PrintData(std::ostream& rOStream) const override
    {
    }

private:

}; // Class PastixComplexSolver



/**
 * input stream function
 */
template <class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::istream& operator>>(std::istream& rIStream,
                                PastixComplexSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis)
{
    return rIStream;
}

/**
 * output stream function
 */
template <class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator<<(std::ostream& rOStream,
                                const PastixComplexSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.

#endif // KRATOS_PASTIX_COMPLEX_SOLVER defined
