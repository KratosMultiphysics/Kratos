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

    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Operators
    ///@{

    /// Default constructor.
    ConditionNumberUtility()
        :mpEigenSolverMax(nullptr), 
         mpEigenSolverMin(nullptr)
    {}
    
    /// Default constructor.
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
     * @param InputMatrix The matrix to compute the eigenvalue
     * @param pEigenSolverMax The solver to get the maximal eigen value
     * @param pEigenSolverMin The solver to get the minimal eigen value
     * @return condition_number The condition number
     */
    
    double GetConditionNumber(
        SparseMatrixType& InputMatrix, 
        LinearSolverType::Pointer pEigenSolverMax,
        LinearSolverType::Pointer pEigenSolverMin
        )
    {
        // The eigen system
        DenseVectorType eigen_values;
        DenseMatrixType eigen_vectors;
        
        const SizeType size_matrix = InputMatrix.size1();
        SparseMatrixType identity_matrix = IdentityMatrix(size_matrix, size_matrix);
        
        pEigenSolverMax->Solve(InputMatrix, identity_matrix, eigen_values, eigen_vectors);
        const double max_lambda = eigen_values[0];

        pEigenSolverMin->Solve(InputMatrix, identity_matrix, eigen_values, eigen_vectors);
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

