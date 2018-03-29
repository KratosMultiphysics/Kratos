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


/* System includes */

/* External includes */

/* Project includes */
#include "linear_solvers/linear_solver.h"
#ifdef INCLUDE_FEAST
    #include "external_includes/feast_solver.h"
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
    
    typedef Matrix                                             MatrixType;

    typedef Vector                                             VectorType;

    typedef std::size_t                                          SizeType;
    
    typedef std::size_t                                         IndexType;

    typedef typename TSparseSpace::MatrixType            SparseMatrixType;

    typedef typename TSparseSpace::VectorType            SparseVectorType;

    typedef typename TDenseSpace::MatrixType              DenseMatrixType;

    typedef typename TDenseSpace::VectorType              DenseVectorType;
    
    typedef std::complex<double>                              ComplexType;
    
    typedef compressed_matrix<ComplexType>        ComplexSparseMatrixType;

    typedef matrix<ComplexType>                    ComplexDenseMatrixType;

    typedef vector<ComplexType>                         ComplexVectorType;

    typedef UblasSpace<ComplexType, ComplexSparseMatrixType, ComplexVectorType> ComplexSparseSpaceType;

    typedef UblasSpace<ComplexType, ComplexDenseMatrixType, ComplexVectorType> ComplexDenseSpaceType;

    typedef LinearSolver<ComplexSparseSpaceType, ComplexDenseSpaceType> ComplexLinearSolverType;
    
    ///@}
    ///@name Life Cycle
    ///@{

    /* Constructor */


    /** Destructor */

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    
    /**
     * Computes the condition number using the maximum and minimum eigenvalue of the system (in moduli)
     * @param InputMatrix: The matrix to obtain the condition number
     * @param pLinearSolver: The complex linear solver considered in the FEAST solver
     * @return condition_number: The condition number obtained
     */
    static inline double ConditionNumber(
        const MatrixType& InputMatrix,
        ComplexLinearSolverType::Pointer pLinearSolver = nullptr
        )
    {
#ifdef INCLUDE_FEAST
        typedef FEASTSolver<TSparseSpace, TDenseSpace> FEASTSolverType;
        
        Parameters this_params(R"(
        {
            "solver_type": "FEAST",
            "print_feast_output": false,
            "perform_stochastic_estimate": true,
            "solve_eigenvalue_problem": true,
            "lambda_min": 0.0,
            "lambda_max": 1.0,
            "echo_level": 0,
            "number_of_eigenvalues": 0,
            "search_dimension": 10
        })");
        
        const std::size_t size = InputMatrix.size1();
        
        const double normA = TSparseSpace::TwoNorm(InputMatrix);
        this_params["lambda_max"].SetDouble(normA);
        this_params["lambda_min"].SetDouble(-normA);
        this_params["number_of_eigenvalues"].SetInt(size * 2/3 - 1);
        this_params["search_dimension"].SetInt(3/2 * size + 1);
        SparseMatrixType copy_matrix = InputMatrix;
        SparseMatrixType identity_matrix = IdentityMatrix(size, size);
        
        // Create the auxilary eigen system
        DenseMatrixType eigen_vectors;
        DenseVectorType eigen_values;
        
        // Create the FEAST solver
        FEASTSolverType FEASTSolver(Kratos::make_shared<Parameters>(this_params), pLinearSolver);
        
        // Solve the problem
        FEASTSolver.Solve(copy_matrix, identity_matrix, eigen_values, eigen_vectors);
        
        // Size of the eigen values vector
        const int dim_eigen_values = eigen_values.size();
        
        // We get the moduli of the eigen values
        #pragma omp parallel for 
        for (int i = 0; i < dim_eigen_values; i++)
        {
            eigen_values[i] = std::abs(eigen_values[i]);
        }
        
        // Now we sort the vector
        std::sort(eigen_values.begin(), eigen_values.end());
        
        // We compute the eigen value
        double condition_number = 0.0;
        if (dim_eigen_values > 0) condition_number = eigen_values[dim_eigen_values - 1]/eigen_values[0];
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

