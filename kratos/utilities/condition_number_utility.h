//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_COND_NUMBER_UTILITY )
#define  KRATOS_COND_NUMBER_UTILITY

/* System includes */

/* External includes */

/* Project includes */
#include "spaces/ublas_space.h"
#include "utilities/math_utils.h"
#ifdef KRATOS_USE_AMATRIX   // This macro definition is for the migration period and to be removed afterward please do not use it
#include "boost/numeric/ublas/matrix.hpp" // for the identity matrix used here.
#else
#endif // KRATOS_USE_AMATRIX

// Linear solvers
#include "linear_solvers/linear_solver.h"
#include "linear_solvers/skyline_lu_factorization_solver.h"
#include "linear_solvers/power_iteration_eigenvalue_solver.h"
#include "linear_solvers/power_iteration_highest_eigenvalue_solver.h"

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

/**
 * @class ConditionNumberUtility
 * @ingroup KratosCore
 * @brief Utility to compute the condition number
 * @details This utility is used in order to compute the condition number of a sparse matrix. Please provide an eigensolver for the maximum and minimum eigenvalue or the power iterations solvers will be considered by default
 * @author Vicente Mataix Ferrandiz
 */
class ConditionNumberUtility
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of ConditionNumberUtility
    KRATOS_CLASS_POINTER_DEFINITION( ConditionNumberUtility );

    /// The sisze type
    typedef std::size_t SizeType;

    /// The index type
    typedef std::size_t IndexType;

    /// The sparse space considered (the one containing the compressed matrix)
    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SparseSpaceType;

    /// The dense space considered
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    /// The compressed matrix
    typedef SparseSpaceType::MatrixType SparseMatrixType;

    /// The vector considered
    typedef SparseSpaceType::VectorType VectorType;

    /// The dense matrix
    typedef LocalSpaceType::MatrixType DenseMatrixType;

    /// The dense vector
    typedef LocalSpaceType::VectorType DenseVectorType;

    /// The definion of the linear solver
    typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;

    /// The reorder considered
    typedef Reorderer<SparseSpaceType,  LocalSpaceType > ReordererType;

    /// Skyline solver definion
    typedef SkylineLUFactorizationSolver<SparseSpaceType,  LocalSpaceType, ReordererType > SkylineLUFactorizationSolverType;

    /// Power iteration solver for the highest eigenvalue
    typedef PowerIterationHighestEigenvalueSolver<SparseSpaceType,  LocalSpaceType, LinearSolverType > PowerIterationHighestEigenvalueSolverType;

    /// Power iteration solver for the lowest eigenvalue
    typedef PowerIterationEigenvalueSolver<SparseSpaceType,  LocalSpaceType, LinearSolverType > PowerIterationEigenvalueSolverType;

    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief Default constructor
     * @details No eigen solvers provided, using the PowerIteration ones as default
     */
    ConditionNumberUtility()
    {
        // definition of empty parameters
        Parameters empty_parameters = Parameters(R"({})" );
        LinearSolverType::Pointer psolver = LinearSolverType::Pointer( new SkylineLUFactorizationSolverType() );
        mpEigenSolverMax = LinearSolverType::Pointer( new PowerIterationHighestEigenvalueSolverType(empty_parameters, psolver) );
        mpEigenSolverMin = LinearSolverType::Pointer( new PowerIterationEigenvalueSolverType(empty_parameters, psolver) );
    }

    /**
     * @brief Default constructor.
     * @param pEigenSolverMax The eigensolver used to determine the highest eigen value
     * @param pEigenSolverMin The eigensolver used to determine the lowest eigen value
     */
    ConditionNumberUtility(
        LinearSolverType::Pointer pEigenSolverMax,
        LinearSolverType::Pointer pEigenSolverMin
        ) :mpEigenSolverMax(pEigenSolverMax),
           mpEigenSolverMin(pEigenSolverMin)
    {}

    /// Destructor
    virtual ~ConditionNumberUtility(){}

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This function computes using the inverse power method the minimal eigenvalue
     * @param rInputMatrix The matrix to compute the eigenvalue
     * @return condition_number The condition number
     */
    double GetConditionNumber(SparseMatrixType& rInputMatrix)
    {
        KRATOS_ERROR_IF(mpEigenSolverMax == nullptr || mpEigenSolverMin == nullptr) << "ERROR:: PLEASE DEFINE THE EigenSolvers" << std::endl;
        return GetConditionNumber(rInputMatrix, mpEigenSolverMax, mpEigenSolverMin);
    }

    /**
     * @brief This function computes using the inverse power method the minimal eigenvalue
     * @param rInputMatrix The matrix to compute the eigenvalue
     * @param pEigenSolverMax The solver to get the maximal eigen value
     * @param pEigenSolverMin The solver to get the minimal eigen value
     * @return condition_number The condition number
     */
    double GetConditionNumber(
        SparseMatrixType& rInputMatrix,
        LinearSolverType::Pointer pEigenSolverMax,
        LinearSolverType::Pointer pEigenSolverMin
        )
    {
        // The eigen system
        DenseVectorType eigen_values;
        DenseMatrixType eigen_vectors;

        const SizeType size_matrix = rInputMatrix.size1();

        SparseMatrixType identity_matrix(size_matrix, size_matrix);
        for (IndexType i = 0; i < size_matrix; ++i)
            identity_matrix.push_back(i, i, 1.0);

        pEigenSolverMax->Solve(rInputMatrix, identity_matrix, eigen_values, eigen_vectors);
        const double max_lambda = eigen_values[0];

        pEigenSolverMin->Solve(rInputMatrix, identity_matrix, eigen_values, eigen_vectors);
        const double min_lambda = eigen_values[0];

        KRATOS_ERROR_IF(min_lambda < std::numeric_limits<double>::epsilon()) << "ERROR:: NOT POSSIBLE TO COMPUTE CONDITION NUMBER. ZERO EIGENVALUE" << std::endl;

        const double condition_number = std::abs(max_lambda)/std::abs(min_lambda);

        return condition_number;
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

    ///@}
    ///@name Friends
    ///@{

private:

    ///@name Private static Member Variables
    ///@{

    ///@}
    ///@name Private member Variables
    ///@{

    LinearSolverType::Pointer mpEigenSolverMax; /// The eigensolver used to determine the highest eigen value
    LinearSolverType::Pointer mpEigenSolverMin; /// The eigensolver used to determine the lowest eigen value

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
    ///@name Private LifeCycle
    ///@{

    ///@}
    ///@name Unaccessible methods
    ///@{

}; /* Class ConditionNumberUtility */

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

}  /* namespace Kratos.*/

#endif /* KRATOS_COND_NUMBER_UTILITY  defined */

