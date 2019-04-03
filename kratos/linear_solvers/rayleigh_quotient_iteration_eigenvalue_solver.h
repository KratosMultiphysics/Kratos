//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//  Collaborator:    Vicente Mataix Ferrandiz
//
//

#if !defined(KRATOS_RAYLEIGH_QUOTIENT_ITERATION_EIGENVALUE_SOLVER_H_INCLUDED )
#define  KRATOS_RAYLEIGH_QUOTIENT_ITERATION_EIGENVALUE_SOLVER_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <numeric>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "linear_solvers/iterative_solver.h"
#include "utilities/sparse_matrix_multiplication_utility.h"
#include "linear_solvers/skyline_lu_factorization_solver.h"

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
 * @class RayleighQuotientIterationEigenvalueSolver
 * @ingroup KratosCore
 * @brief This is a eigen solver based on the Rayleigh quotient iteration algorithm
 * @details Rayleigh quotient iteration is an iterative method, that is, it delivers a sequence of approximate solutions that converges to a true solution in the limit (this is true for all algorithms that compute eigenvalues: since eigenvalues can be irrational numbers, there can be no general method for computing them in a finite number of steps). Very rapid convergence is guaranteed and no more than a few iterations are needed in practice to obtain a reasonable approximation. The Rayleigh quotient iteration algorithm converges cubically for Hermitian or symmetric matrices, given an initial vector that is sufficiently close to an eigenvector of the matrix that is being analyzed.
 * @see https://en.wikipedia.org/wiki/Rayleigh_quotient_iteration
 * @author Pooyan Dadvand
*/
template<class TSparseSpaceType, class TDenseSpaceType, class TLinearSolverType,
         class TPreconditionerType = Preconditioner<TSparseSpaceType, TDenseSpaceType>,
         class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class RayleighQuotientIterationEigenvalueSolver : public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RayleighQuotientIterationEigenvalueSolver
    KRATOS_CLASS_POINTER_DEFINITION(RayleighQuotientIterationEigenvalueSolver);

    /// The base class definition
    typedef IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType> BaseType;

    /// The sparse matrix defintion
    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    /// The vector definition ("sparse")
    typedef typename TSparseSpaceType::VectorType VectorType;

    /// The dense matrix definition
    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    /// The "dense" vector definition
    typedef typename TDenseSpaceType::VectorType DenseVectorType;

    /// The size type definiton
    typedef std::size_t SizeType;

    /// The index type definition
    typedef std::size_t IndexType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    RayleighQuotientIterationEigenvalueSolver()
    {

    }

    /**
     * @brief The "manual" settings constructor
     * @param NewMaxTolerance The tolerance considered
     * @param NewMaxIterationsNumber The maximum number of iterations considered
     * @param NewRequiredEigenvalueNumber The number of eigen values to compute
     * @param ShiftingConvergence The convergence parameter of shifting
     */
    RayleighQuotientIterationEigenvalueSolver(
        double NewMaxTolerance,
        unsigned int NewMaxIterationsNumber,
        unsigned int NewRequiredEigenvalueNumber,
        typename TLinearSolverType::Pointer pLinearSolver,
        double ShiftingConvergence = 0.25
    ): BaseType(NewMaxTolerance, NewMaxIterationsNumber),
       mRequiredEigenvalueNumber(NewRequiredEigenvalueNumber),
       mpLinearSolver(pLinearSolver),
       mShiftingConvergence(ShiftingConvergence)
    {

    }

    /**
     * @brief The parameters constructor
     * @param ThisParameters The input parameters
     * @param pLinearSolver The linear solver considered
     */
    RayleighQuotientIterationEigenvalueSolver(
        Parameters ThisParameters,
        typename TLinearSolverType::Pointer pLinearSolver
        ): mpLinearSolver(pLinearSolver)
    {
        Parameters DefaultParameters = Parameters(R"(
            {
                "solver_type"             : "rayleigh_quotient_iteration_eigenvalue_solver",
                "max_iteration"           : 500,
                "tolerance"               : 1e-9,
                "required_eigen_number"   : 1,
                "shifting_convergence"    : 0.25,
                "verbosity"               : 1,
                "linear_solver_settings"  : {}
            })" );

        ThisParameters.ValidateAndAssignDefaults(DefaultParameters);

        mRequiredEigenvalueNumber = ThisParameters["required_eigen_number"].GetInt();
        mShiftingConvergence = ThisParameters["shifting_convergence"].GetDouble();
        mEchoLevel = ThisParameters["verbosity"].GetInt();
        BaseType::SetTolerance( ThisParameters["tolerance"].GetDouble() );
        BaseType::SetMaxIterationsNumber( ThisParameters["max_iteration"].GetInt() );
    }

    /// Copy constructor.
    RayleighQuotientIterationEigenvalueSolver(
        const RayleighQuotientIterationEigenvalueSolver& Other) :
        BaseType(Other),
        mRequiredEigenvalueNumber(Other.mRequiredEigenvalueNumber), mpLinearSolver(Other.mpLinearSolver),
        mShiftingConvergence(Other.mShiftingConvergence)
    {

    }

    /// Destructor.
    virtual ~RayleighQuotientIterationEigenvalueSolver() {}

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    RayleighQuotientIterationEigenvalueSolver& operator=(const RayleighQuotientIterationEigenvalueSolver& Other)
    {
        BaseType::operator=(Other);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method initializes the system
     * @details It computes the vector R, which contains the components of the diagonal of the M matrix
     * @param rR The vector containing the normalized components of the diagonal
     * @param rM The "mass" matrix
     */
    static void Initialize(
        VectorType& rR,
        const SparseMatrixType& rM
        )
    {
        const SizeType size_m = rM.size1();

        // Resize in case of not same size
        if (rR.size() != size_m)
            rR.resize(size_m);

        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(size_m); ++i) {
            rR[i] = rM(i,i);
        }

        KRATOS_ERROR_IF(norm_2(rR) == 0.0) << "Invalid M matrix. The norm2 of its diagonal is Zero" << std::endl;

        rR /= norm_2(rR);
    }

    /**
     * @brief This method performs a Sturm Sequence Check
     * @param ShiftedK The modified K matrix after apply the M matrix
     */
    SizeType SturmSequenceCheck(SparseMatrixType& ShiftedK)
    {
        // define an object to store skyline matrix and factorization
        LUSkylineFactorization<TSparseSpaceType, TDenseSpaceType> myFactorization;

        // copy myMatrix into skyline format
        myFactorization.copyFromCSRMatrix(ShiftedK);

        // factorize it
        myFactorization.factorize();
        SizeType counter = 0; // number of eigenvalues less than the shift
        for(SizeType i = 0 ; i < ShiftedK.size1() ; i++){
            if(myFactorization.entriesD[i] < 0.0) {
                counter++;
            }
        }

        return counter;
    }

    /**
     * @brief The Rayleigh quotient iteration method
     * @param K: The stiffness matrix
     * @param M: The mass matrix
     * @param Eigenvalues: The vector containing the eigen values
     * @param Eigenvectors: The matrix containing the eigen vectors
     */
    void Solve(
        SparseMatrixType& K,
        SparseMatrixType& M,
        DenseVectorType& Eigenvalues,
        DenseMatrixType& Eigenvectors
    ) override
    {
        SizeType size = K.size1();
        SizeType max_iteration = BaseType::GetMaxIterationsNumber();
        double tolerance = BaseType::GetTolerance();

        VectorType x = boost::numeric::ublas::zero_vector<double>(size);
        VectorType y = boost::numeric::ublas::zero_vector<double>(size);

        Initialize(y,M);

        if(Eigenvalues.size() < 1) {
            Eigenvalues.resize(1,0.0);
        }

        const double epsilon = 1.0e-9;
        // Starting with first step
        double beta = 0.0;
        double ro = 0.0;
        double shift_value = 0.0;
        double old_ro = 0.0;//Eigenvalues[0];

        KRATOS_INFO_IF("RayleighQuotientIterationEigenvalueSolver", mEchoLevel > 1) << "Iteration    beta \t ro \t\t convergence norm \t min \t\t max" << std::endl;

        SparseMatrixType shifted_k(K);
        // define an object to store skyline matrix and factorization
        LUSkylineFactorization<TSparseSpaceType, TDenseSpaceType> my_factorization;

        // Copy myMatrix into skyline format
        my_factorization.copyFromCSRMatrix(shifted_k);

        // Factorize it
        my_factorization.factorize();

        double min_shift_value = 0.0;
        double max_shift_value = 0.0;

        SizeType smaller_eigenvalue_numbers = 0;

        for(SizeType i = 0 ; i < max_iteration ; ++i) {
            //K*x = y
            //mpLinearSolver->Solve(shifted_k,x,y);
            my_factorization.backForwardSolve(size, y, x);

            ro = inner_prod(y,x);

            //y = M*x
            TSparseSpaceType::Mult(M, x, y);

            beta = inner_prod(x, y);
            KRATOS_ERROR_IF(beta == 0.0) << "Zero beta norm!" << std::endl;

            const double delta_ro = (ro / beta);

            ro = delta_ro + shift_value;

            KRATOS_ERROR_IF(ro == 0.0) << "Perpendicular eigenvector to M" << std::endl;

            double convergence_norm = std::abs((ro - old_ro) / ro);

            if(convergence_norm < mShiftingConvergence) { // Start shifting after certain convergence
                // If there are no smaller eigenvalues yet we need to extend the range
                if(smaller_eigenvalue_numbers == 0) {
                    max_shift_value = ro;
                }

                if((ro > max_shift_value)||(ro < min_shift_value)) {
                    shift_value = (max_shift_value + min_shift_value) / 2.0;
                } else {
                    shift_value = ro;
                }

                SparseMatrixMultiplicationUtility::MatrixAdd<SparseMatrixType, SparseMatrixType>(shifted_k, M, - shift_value);

                // Copy myMatrix into skyline format
                my_factorization.copyFromCSRMatrix(shifted_k);

                // Factorize it
                my_factorization.factorize();

                SizeType new_smaller_eigenvalue_numbers = SturmSequenceCheck(shifted_k);

                if(new_smaller_eigenvalue_numbers == 0) {
                    min_shift_value = shift_value;
                } else {
                    max_shift_value = shift_value;
                    smaller_eigenvalue_numbers = new_smaller_eigenvalue_numbers;
                }

                unsigned int iteration_number = 0;
                unsigned int max_shift_number = 4;

                while((smaller_eigenvalue_numbers > 1) && (max_shift_value-min_shift_value > epsilon) && (iteration_number++ < max_shift_number)) {
                    shift_value = (max_shift_value + min_shift_value) / 2.0;
                    shifted_k = K;
                    SparseMatrixMultiplicationUtility::MatrixAdd<SparseMatrixType, SparseMatrixType>(shifted_k, M, - shift_value);

                    // Copy myMatrix into skyline format
                    my_factorization.copyFromCSRMatrix(shifted_k);

                    // Factorize it
                    my_factorization.factorize();

                    new_smaller_eigenvalue_numbers = SturmSequenceCheck(shifted_k);

                    if(new_smaller_eigenvalue_numbers == 0) {
                        min_shift_value = shift_value;
                        KRATOS_INFO_IF("RayleighQuotientIterationEigenvalueSolver", mEchoLevel > 1) << "            Finding " << smaller_eigenvalue_numbers << " eigenvalues in [" << min_shift_value << "," << max_shift_value  << "]" << std::endl;
                    } else {
                        max_shift_value = shift_value;
                        smaller_eigenvalue_numbers = new_smaller_eigenvalue_numbers;
                        KRATOS_INFO_IF("RayleighQuotientIterationEigenvalueSolver", mEchoLevel > 1) <<  "            Finding " << smaller_eigenvalue_numbers << " eigenvalues in [" << min_shift_value << "," << max_shift_value  << "]" << std::endl;
                    }
                }

            }

            if(beta < 0.0) {
                beta = -std::sqrt(-beta);
            } else {
                beta = std::sqrt(beta);
            }

            TSparseSpaceType::InplaceMult(y, 1.0/beta);

            KRATOS_INFO_IF("RayleighQuotientIterationEigenvalueSolver", mEchoLevel > 1) << i << " \t " << beta << " \t " << ro << " \t " << convergence_norm  << " \t\t " <<  min_shift_value << " \t " << max_shift_value << std::endl;

            if(convergence_norm < tolerance) {
                break;
            }

            old_ro = ro;
        }

        KRATOS_INFO_IF("RayleighQuotientIterationEigenvalueSolver", mEchoLevel > 0) << "ro:\n" << ro << std::endl << std::endl << "y:\n" << y << std::endl;

        Eigenvalues[0] = ro;

        if((Eigenvectors.size1() < 1) || (Eigenvectors.size2() < size)) {
            Eigenvectors.resize(1,size);
        }

        for(SizeType i = 0 ; i < size ; i++) {
            Eigenvectors(0,i) = x[i] / beta;
        }
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

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Power iteration eigenvalue solver with " << BaseType::GetPreconditioner()->Info();
        return  buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        BaseType::PrintData(rOStream);
    }

    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    unsigned int mRequiredEigenvalueNumber;

    unsigned int mEchoLevel;

    typename TLinearSolverType::Pointer mpLinearSolver;

    double mShiftingConvergence;

    std::vector<DenseVectorType> mQVector;
    std::vector<DenseVectorType> mPVector;
    std::vector<DenseVectorType> mRVector;

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
    ///@name Un accessible methods
    ///@{


    ///@}

}; // Class RayleighQuotientIterationEigenvalueSolver
///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TSparseSpaceType, class TDenseSpaceType,
         class TPreconditionerType,
         class TReordererType>
inline std::istream& operator >> (std::istream& IStream,
                                  RayleighQuotientIterationEigenvalueSolver<TSparseSpaceType, TDenseSpaceType,
                                  TPreconditionerType, TReordererType>& rThis)
{
    return IStream;
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType,
         class TPreconditionerType,
         class TReordererType>
inline std::ostream& operator << (std::ostream& OStream,
                                  const RayleighQuotientIterationEigenvalueSolver<TSparseSpaceType, TDenseSpaceType,
                                  TPreconditionerType, TReordererType>& rThis)
{
    rThis.PrintInfo(OStream);
    OStream << std::endl;
    rThis.PrintData(OStream);

    return OStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_RAYLEIGH_QUOTIENT_ITERATION_EIGENVALUE_SOLVER_H_INCLUDED defined































