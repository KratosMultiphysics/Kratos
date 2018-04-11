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
#include "includes/ublas_interface.h"
#include "utilities/math_utils.h"
#include "linear_solvers/linear_solver.h"
#include "linear_solvers/power_iteration_eigenvalue_solver.h"

#include "linear_solvers/skyline_lu_factorization_solver.h"
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

///Utility to compute the condition number
/**
 * Utility needed to compute the condition number
 */
class ConditionNumberUtility
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of ConditionNumberUtility
    KRATOS_CLASS_POINTER_DEFINITION( ConditionNumberUtility );

    typedef std::size_t SizeType;
    
    typedef unsigned int IndexType;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    typedef SparseSpaceType::MatrixType SparseMatrixType;

    typedef SparseSpaceType::VectorType VectorType;

    typedef LocalSpaceType::MatrixType DenseMatrixType;

    typedef LocalSpaceType::VectorType DenseVectorType;

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

    /// Default constructor.
    ConditionNumberUtility() 
    {
        //Parameters empty_parameters = Parameters(R"({})");
        Parameters settings_max = Parameters(R"(
        {
            "solver_type"             : "power_iteration_highest_eigenvalue_solver",
            "max_iteration"           : 10000,
            "tolerance"               : 1e-9,
            "required_eigen_number"   : 1,
            "verbosity"               : 0,
            "linear_solver_settings"  : {
                "solver_type"             : "SuperLUSolver",
                "max_iteration"           : 500,
                "tolerance"               : 1e-9,
                "scaling"                 : false,
                "verbosity"               : 0
            }
        }
        )");
        Parameters settings_min = Parameters(R"(
        {
            "solver_type"             : "power_iteration_eigenvalue_solver",
            "max_iteration"           : 10000,
            "tolerance"               : 1e-9,
            "required_eigen_number"   : 1,
            "verbosity"               : 0,
            "linear_solver_settings"  : {
                "solver_type"             : "SuperLUSolver",
                "max_iteration"           : 500,
                "tolerance"               : 1e-9,
                "scaling"                 : false,
                "verbosity"               : 0
            }
        }
        )");
        LinearSolverType::Pointer psolver = LinearSolverType::Pointer ( new SkylineLUFactorizationSolverType() );
        mpEigenSolverMax = LinearSolverType::Pointer( new PowerIterationHighestEigenvalueSolverType(settings_max, psolver));
        mpEigenSolverMin = LinearSolverType::Pointer( new PowerIterationEigenvalueSolverType(settings_min, psolver));
        
    }

    /// Default constructor.
    ConditionNumberUtility(
        LinearSolverType::Pointer pEigenSolverMax,
        LinearSolverType::Pointer pEigenSolverMin
        ) :mpEigenSolverMax(pEigenSolverMax), 
           mpEigenSolverMin(pEigenSolverMin)
    {std::cout << "\n" << "I ENTER second default constructor CONDITION NUMBER"  << std::endl;}

    /// Destructor
    virtual ~ConditionNumberUtility(){}

    ///@}
    ///@name Operations
    ///@{

    /**
     * This function computes using the inverse power method the minimal eigenvalue
     * @param InputMatrix The matrix to compute the eigenvalue
     * @return condition_number The condition number
     */
    
    double GetConditionNumber(SparseMatrixType& InputMatrix)
    {
        KRATOS_ERROR_IF(mpEigenSolverMax == nullptr || mpEigenSolverMin == nullptr) << "ERROR:: PLEASE DEFINE THE EigenSolvers" << std::endl;
        return GetConditionNumber(InputMatrix, mpEigenSolverMax, mpEigenSolverMin);
    }
    
    /**
     * This function computes using the inverse power method the minimal eigenvalue
     * @param InputMatrix: The matrix to compute the eigenvalue
     * @param pEigenSolverMax: The solver to get the maximal eigen value
     * @param pEigenSolverMin: The solver to get the minimal eigen value
     * @return condition_number: The condition number
     */
    
    double GetConditionNumber(
        SparseMatrixType& InputMatrix, 
        LinearSolverType::Pointer pEigenSolverMax,
        LinearSolverType::Pointer pEigenSolverMin
        )
    {
        std::cout << "\n" << "I ENTER GET CONDITION NUMBER"  << std::endl;
        
        DenseVectorType eigen_values;
        DenseMatrixType eigen_vectors;
        
        const SizeType size_matrix = InputMatrix.size1();
        SparseMatrixType identity_matrix = IdentityMatrix(size_matrix, size_matrix);
        
        pEigenSolverMax->Solve(InputMatrix, identity_matrix, eigen_values, eigen_vectors);
        const double max_lambda = eigen_values[0];

        std::cout << "\n" << "max_lambda = " << std::scientific << max_lambda  << std::endl;

        pEigenSolverMin->Solve(InputMatrix, identity_matrix, eigen_values, eigen_vectors);
        const double min_lambda = eigen_values[0];

        std::cout << "\n" << "min_lambda = " << std::scientific << min_lambda  << std::endl;
        
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

    LinearSolverType::Pointer mpEigenSolverMax;
    LinearSolverType::Pointer mpEigenSolverMin;

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

