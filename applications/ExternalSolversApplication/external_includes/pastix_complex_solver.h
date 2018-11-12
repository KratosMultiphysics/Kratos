//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//   License:		 BSD License
//   Kratos default license: kratos/license.txt
//
//   Project Name:        $ExternalSolversApplication   $
//   Last modified by:    $Author: michael.andre@tum.de $
//   Date:                $Date:             April 2017 $
//   Revision:            $Revision:                0.0 $
//
//

#if !defined(KRATOS_PASTIX_COMPLEX_SOLVER)
#define KRATOS_PASTIX_COMPLEX_SOLVER

// System includes
#ifdef _OPENMP
  #include <omp.h>
#endif
#include <complex>
#include <vector>
#include <algorithm>

// External includes
extern "C" {
#include <pastix.h>
}

// Project includes
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "includes/kratos_parameters.h"
#include "linear_solvers/linear_solver.h"
#include "linear_solvers/direct_solver.h"

namespace Kratos
{
///@name Kratos Classes
///@{

template <class TSparseSpaceType, class TDenseSpaceType, class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType>>
class PastixComplexSolver : public DirectSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>
{
public:
	///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(PastixComplexSolver);

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

	///@}
    ///@name Life Cycle
    ///@{

    PastixComplexSolver(Parameters rSettings)
    {
        mpPastixData = nullptr;

        Parameters default_settings(R"(
        {
            "solver_type" : "pastix",
            "echo_level" : 0
        })");

        rSettings.ValidateAndAssignDefaults(default_settings);

        mEchoLevel = rSettings["echo_level"].GetInt();
        if (mEchoLevel > 4 || mEchoLevel < 0)
        {
            KRATOS_ERROR << "Invalid echo_level: " << mEchoLevel << std::endl;
        }
    }

    PastixComplexSolver(const PastixComplexSolver& Other) = delete;

    ~PastixComplexSolver() override
	{
	}

	///@}
    ///@name Operators
    ///@{

    PastixComplexSolver& operator=(const PastixComplexSolver& Other) = delete;

	///@}
    ///@name Operations
    ///@{

    /** This function is designed to be called as few times as possible. It creates the data structures
     * that only depend on the connectivity of the matrix (and not on its coefficients)
     * so that the memory can be allocated once and expensive operations can be done only when strictly
     * needed
     * @param rA. System matrix
     * @param rX. Solution vector. it's also the initial guess for iterative linear solvers.
     * @param rB. Right hand side vector.
     */
    void Initialize(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        if (mpPastixData != nullptr)
            this->Clear();

        mNRows = rA.size1();

        if (mRowptr.size() != rA.index1_data().size())
            mRowptr.resize(rA.index1_data().size());

        if (mCol.size() != rA.index2_data().size())
            mCol.resize(rA.index2_data().size());

        if (mPerm.size() != static_cast<unsigned int>(mNRows))
            mPerm.resize(mNRows);
            
        if (mInvp.size() != static_cast<unsigned int>(mNRows))
            mInvp.resize(mNRows);

        // initialize iparm and dparm values
        mIparm[IPARM_MODIFY_PARAMETER] = API_NO;
        mIparm[IPARM_START_TASK      ] = API_TASK_INIT;
        mIparm[IPARM_END_TASK        ] = API_TASK_INIT;

        z_pastix(&mpPastixData, 0, mNRows, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, 1, mIparm, mDparm);
    }

    /** This function is designed to be called every time the coefficients change in the system
     * that is, normally at the beginning of each solve.
     * For example if we are implementing a direct solver, this is the place to do the factorization
     * so that then the backward substitution can be performed effectively more than once
     * @param rA. System matrix
     * @param rX. Solution vector. it's also the initial guess for iterative linear solvers.
     * @param rB. Right hand side vector.
     */
    void InitializeSolutionStep(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        if (mpPastixData != nullptr)
            KRATOS_ERROR << "pastix_data != nullptr upon entering InitializeSolutionStep." << std::endl;

        // 1-based indexing for pastix
        std::copy(std::begin(rA.index1_data()),std::end(rA.index1_data()),std::begin(mRowptr));
        for (auto &i : mRowptr)
            ++i;
        std::copy(std::begin(rA.index2_data()),std::end(rA.index2_data()),std::begin(mCol));
        for (auto &i : mCol)
            ++i;

        // factorize system matrix
        mIparm[IPARM_VERBOSE        ] = static_cast<API_VERBOSE>(mEchoLevel);
        mIparm[IPARM_RHS_MAKING     ] = API_RHS_B; // user-provided rhs
        mIparm[IPARM_SYM            ] = API_SYM_NO; // non-symmetric
        mIparm[IPARM_FACTORIZATION  ] = API_FACT_LU; // LU factorization
        mIparm[IPARM_TRANSPOSE_SOLVE] = API_YES; // solve transpose for csr matrix format
    #ifdef _OPENMP
        mIparm[IPARM_THREAD_NBR     ] = omp_get_max_threads();
    #endif

        mIparm[IPARM_START_TASK] = API_TASK_ORDERING;
        mIparm[IPARM_END_TASK  ] = API_TASK_NUMFACT;

        z_pastix(&mpPastixData, 0, mNRows,
            static_cast<pastix_int_t*>(&mRowptr[0]),
            static_cast<pastix_int_t*>(&mCol[0]),
            static_cast<std::complex<double>*>(&rA.value_data()[0]),
            static_cast<pastix_int_t*>(&mPerm[0]),
            static_cast<pastix_int_t*>(&mInvp[0]),
            nullptr,
            1, mIparm, mDparm);
    }

    /** 
     * @brief This function actually performs the solution work, eventually taking advantage of what was done before in the Initialize and InitializeSolutionStep functions.
     * @param rA. System matrix
     * @param rX. Solution vector. it's also the initial guess for iterative linear solvers.
     * @param rB. Right hand side vector.
     */
    void PerformSolutionStep(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        if (mpPastixData == nullptr)
            KRATOS_ERROR << "pastix_data == nullptr upon entering Solve." << std::endl;

        noalias(rX) = rB;
        mIparm[IPARM_START_TASK] = API_TASK_SOLVE;
        mIparm[IPARM_END_TASK  ] = API_TASK_SOLVE;

        z_pastix(&mpPastixData, 0, mNRows,
            static_cast<pastix_int_t*>(&mRowptr[0]),
            static_cast<pastix_int_t*>(&mCol[0]),
            static_cast<std::complex<double>*>(&rA.value_data()[0]),
            static_cast<pastix_int_t*>(&mPerm[0]),
            static_cast<pastix_int_t*>(&mInvp[0]),
            static_cast<std::complex<double>*>(&rX[0]),
            1, mIparm, mDparm);
    }

    /** 
     * @brief Normal solve method.
     * @details Solves the linear system Ax=b and puts the result on SystemVector& rX. rVectorx is also th initial guess for iterative methods.
     * @param rA. System matrix
     * @param rX. Solution vector. it's also the initial guess for iterative linear solvers.
     * @param rB. Right hand side vector.
     */
    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
	    	InitializeSolutionStep(rA, rX, rB);
        PerformSolutionStep(rA, rX, rB);
        FinalizeSolutionStep(rA, rX, rB);

        return true;
    }

    /** This function is designed to be called at the end of the solve step.
     * for example this is the place to remove any data that we do not want to save for later
     * @param rA. System matrix
     * @param rX. Solution vector. it's also the initial guess for iterative linear solvers.
     * @param rB. Right hand side vector.
     */
    void FinalizeSolutionStep(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        this->Clear();
    }

    /** This function is designed to clean up all internal data in the solver.
     * Clear is designed to leave the solver object as if newly created.
     * After a clear a new Initialize is needed
     */
    void Clear() override
    {
        if (mpPastixData != nullptr)
        {
            mIparm[IPARM_START_TASK] = API_TASK_CLEAN;
            mIparm[IPARM_END_TASK  ] = API_TASK_CLEAN;

            z_pastix(&mpPastixData, 0, mNRows,
                static_cast<pastix_int_t*>(&mRowptr[0]),
                static_cast<pastix_int_t*>(&mCol[0]),
                nullptr,
                static_cast<pastix_int_t*>(&mPerm[0]),
                static_cast<pastix_int_t*>(&mInvp[0]),
                nullptr,
                1, mIparm, mDparm);

            mpPastixData = nullptr;
        }
    }

	///@}
    ///@name Input and output
    ///@{

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Pastix direct solver finished.";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

	///@}

private:
	///@name Member Variables
    ///@{
    
    pastix_data_t* mpPastixData;
    pastix_int_t mNRows;
    std::vector<pastix_int_t> mRowptr;
    std::vector<pastix_int_t> mCol;
    std::vector<pastix_int_t> mPerm;
    std::vector<pastix_int_t> mInvp;
    pastix_int_t mIparm[IPARM_SIZE];
    double mDparm[DPARM_SIZE];
    int mEchoLevel;

    ///@}
}; // Class PastixComplexSolver

///@}

///@name Input and output
///@{

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

///@}

} // namespace Kratos.

#endif // KRATOS_PASTIX_COMPLEX_SOLVER defined
