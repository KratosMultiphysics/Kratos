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

#if !defined(KRATOS_FEAST_CONDITION_NUMBER_UTILITY )
#define  KRATOS_FEAST_CONDITION_NUMBER_UTILITY

// System includes

// External includes

// Project includes
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#ifdef USE_EIGEN_FEAST
    #include "custom_solvers/feast_eigensystem_solver.h"
#endif

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

/// This utility uses the FEAST solver to obtain (estimate) the the condition number of a regular matrix
/**
 * Regular matrix: A*A^H=A^H*A
 */
template<class TSparseSpace = UblasSpace<double, CompressedMatrix, Vector>,
         class TDenseSpace = UblasSpace<double, Matrix, Vector>
         >
class FEASTConditionNumberUtility
{
public:

    ///@name Type Definitions
    ///@{

    /// Definition of the shared pointer of the class
    KRATOS_CLASS_POINTER_DEFINITION(FEASTConditionNumberUtility);

    /// Indexes
    typedef std::size_t                                          SizeType;
    typedef std::size_t                                         IndexType;

    /// Sparse space
    typedef typename TSparseSpace::MatrixType            SparseMatrixType;
    typedef typename TSparseSpace::VectorType            SparseVectorType;

    /// Dense space
    typedef typename TDenseSpace::MatrixType              DenseMatrixType;
    typedef typename TDenseSpace::VectorType              DenseVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Computes the condition number using the maximum and minimum eigenvalue of the system (in moduli)
     * @param InputMatrix: The matrix to obtain the condition number
     * @param pLinearSolver: The complex linear solver considered in the FEAST solver
     * @return condition_number: The condition number obtained
     */
    static inline double GetConditionNumber(const SparseMatrixType& InputMatrix)
    {
#ifdef USE_EIGEN_FEAST
        typedef FEASTEigensystemSolver<true, double, double> FEASTSolverType;

        Parameters this_params(R"(
        {
            "solver_type"                : "feast",
            "symmetric"                  : true,
            "number_of_eigenvalues"      : 3,
            "search_lowest_eigenvalues"  : true,
            "search_highest_eigenvalues" : false,
            "e_min"                      : 0.0,
            "e_max"                      : 1.0,
            "echo_level"                 : 0
        })");

        const std::size_t size_matrix = InputMatrix.size1();

        const double normA = TSparseSpace::TwoNorm(InputMatrix);

        this_params["e_max"].SetDouble(normA);
        this_params["e_min"].SetDouble(-normA);
//         this_params["number_of_eigenvalues"].SetInt(size_matrix * 2/3 - 1);
//         this_params["subspace_size"].SetInt(3/2 * size_matrix + 1);
        SparseMatrixType copy_matrix = InputMatrix;
        SparseMatrixType identity_matrix(size_matrix, size_matrix);
        for (IndexType i = 0; i < size_matrix; ++i)
            identity_matrix.push_back(i, i, 1.0);

        // Create the auxilary eigen system
        DenseMatrixType eigen_vectors(size_matrix, 1);
        DenseVectorType eigen_values(size_matrix);

        // Create the FEAST solver
        FEASTSolverType feast_solver_lowest(this_params);

        // Solve the problem
        feast_solver_lowest.Solve(copy_matrix, identity_matrix, eigen_values, eigen_vectors);

        // Size of the eigen values vector
        int dim_eigen_values = eigen_values.size();

        // We get the moduli of the eigen values
        #pragma omp parallel for
        for (int i = 0; i < dim_eigen_values; i++) {
            eigen_values[i] = std::abs(eigen_values[i]);
        }

        // Now we sort the vector
        std::sort(eigen_values.begin(), eigen_values.end());

        const double lowest_eigen_value = eigen_values[0];

        // Create the FEAST solver
        this_params["search_lowest_eigenvalues"].SetBool(false);
        this_params["search_highest_eigenvalues"].SetBool(true);
        FEASTSolverType feast_solver_highest(this_params);

        // Solve the problem
        copy_matrix = InputMatrix;
        feast_solver_highest.Solve(copy_matrix, identity_matrix, eigen_values, eigen_vectors);

        // Size of the eigen values vector
        dim_eigen_values = eigen_values.size();

        // We get the moduli of the eigen values
        #pragma omp parallel for
        for (int i = 0; i < dim_eigen_values; i++) {
            eigen_values[i] = std::abs(eigen_values[i]);
        }

        // Now we sort the vector
        std::sort(eigen_values.begin(), eigen_values.end());

        const double highest_eigen_value = eigen_values[dim_eigen_values - 1];

        // We compute the eigen value
        const double condition_number = highest_eigen_value/lowest_eigen_value;
#else
        const double condition_number = 0.0;
        KRATOS_ERROR << "YOU MUST COMPILE FEAST IN ORDER TO USE THIS UTILITY" << std::endl;
#endif

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

    FEASTConditionNumberUtility(void);

    FEASTConditionNumberUtility(FEASTConditionNumberUtility& rSource);

}; /* Class FEASTConditionNumberUtility */

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

}  /* namespace Kratos.*/

#endif /* KRATOS_FEAST_CONDITION_NUMBER_UTILITY  defined */

