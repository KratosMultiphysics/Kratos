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

/// Short class definition.
/** Detail class definition.
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

    typedef IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    typedef typename TDenseSpaceType::VectorType DenseVectorType;

    typedef std::size_t SizeType;

    typedef std::size_t IndexType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    RayleighQuotientIterationEigenvalueSolver()
    {

    }

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

    static void Initialize(
        DenseVectorType& R,
        SparseMatrixType& M
    )
    {
        for(SizeType i = 0 ; i < R.size() ; i++)
        {
            R[i] = M(i,i);
        }

        if(norm_2(R) == 0.0)
        {
            KRATOS_ERROR << "Invalid M matrix. The norm2 of its diagonal is Zero" << std::endl;
        }

        R /= norm_2(R);
    }

    SizeType SturmSequenceCheck(SparseMatrixType& ShiftedK)
    {
        // define an object to store skyline matrix and factorization
        LUSkylineFactorization<TSparseSpaceType, TDenseSpaceType> myFactorization;

        // copy myMatrix into skyline format
        myFactorization.copyFromCSRMatrix(ShiftedK);

        // factorize it
        myFactorization.factorize();
        SizeType counter = 0; // number of eigenvalues less than the shift
        for(SizeType i = 0 ; i < ShiftedK.size1() ; i++)
        {
            if(myFactorization.entriesD[i] < 0.00)
            {
                counter++;
            }
        }

        return counter;
    }

    /** The Rayleigh quotient iteration method
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
    )
    {
        using boost::numeric::ublas::trans;

        SizeType size = K.size1();
        SizeType max_iteration = BaseType::GetMaxIterationsNumber();
        double tolerance = BaseType::GetTolerance();

        VectorType x = ZeroVector(size);
        VectorType y = ZeroVector(size);

        Initialize(y,M);

        if(Eigenvalues.size() < 1)
        {
            Eigenvalues.resize(1,0.0);
        }

        const double epsilon = 1.0e-9;
        // Starting with first step
        double beta = 0.0;
        double ro = 0.0;
        double shift_value = 0.0;
        double old_ro = 0.0;//Eigenvalues[0];

        if (mEchoLevel > 1)
        {
            std::cout << "Iteration    beta \t ro \t\t convergence norm \t min \t\t max" << std::endl;
        }

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

        for(SizeType i = 0 ; i < max_iteration ; i++)
        {
            //K*x = y
            //mpLinearSolver->Solve(shifted_k,x,y);
            my_factorization.backForwardSolve(size, y, x);

            ro = inner_prod(y,x);

            //y = M*x
            TSparseSpaceType::Mult(M, x, y);

            beta = inner_prod(x, y);
            if(beta == 0.0)
            {
                KRATOS_ERROR << "Zero beta norm!" << std::endl;
            }

            const double delta_ro = (ro / beta);

            ro = delta_ro + shift_value;

            if(ro == 0.0)
            {
                KRATOS_ERROR << "Perpendicular eigenvector to M" << std::endl;
            }

            double convergence_norm = std::abs((ro - old_ro) / ro);

            if(convergence_norm < mShiftingConvergence) // Start shifting after certain convergence
            {
                // If there are no smaller eigenvalues yet we need to extend the range
                if(smaller_eigenvalue_numbers == 0)
                {
                    max_shift_value = ro;
                }

                if((ro > max_shift_value)||(ro < min_shift_value))
                {
                    shift_value = (max_shift_value + min_shift_value) / 2.00;
                }
                else
                {
                    shift_value = ro;
                }

                noalias(shifted_k) = K - shift_value*M;

                // Copy myMatrix into skyline format
                my_factorization.copyFromCSRMatrix(shifted_k);

                // Factorize it
                my_factorization.factorize();

                SizeType new_smaller_eigenvalue_numbers = SturmSequenceCheck(shifted_k);

                if(new_smaller_eigenvalue_numbers == 0)
                {
                    min_shift_value = shift_value;
                }
                else
                {
                    max_shift_value = shift_value;
                    smaller_eigenvalue_numbers = new_smaller_eigenvalue_numbers;
                }

                unsigned int iteration_number = 0;
                unsigned int max_shift_number = 4;

                while((smaller_eigenvalue_numbers > 1) && (max_shift_value-min_shift_value > epsilon) && (iteration_number++ < max_shift_number))
                {
                    shift_value = (max_shift_value + min_shift_value) / 2.00;
                    noalias(shifted_k) = K - shift_value*M;

                    // Copy myMatrix into skyline format
                    my_factorization.copyFromCSRMatrix(shifted_k);

                    // Factorize it
                    my_factorization.factorize();

                    new_smaller_eigenvalue_numbers = SturmSequenceCheck(shifted_k);

                    if(new_smaller_eigenvalue_numbers == 0)
                    {
                        min_shift_value = shift_value;
                        if (mEchoLevel > 1)
                        {
                            std::cout << "            Finding " << smaller_eigenvalue_numbers << " eigenvalues in [" << min_shift_value << "," << max_shift_value  << "]" << std::endl;
                        }
                    }
                    else
                    {
                        max_shift_value = shift_value;
                        smaller_eigenvalue_numbers = new_smaller_eigenvalue_numbers;
                        if (mEchoLevel > 1)
                        {
                            std::cout << "            Finding " << smaller_eigenvalue_numbers << " eigenvalues in [" << min_shift_value << "," << max_shift_value  << "]" << std::endl;
                        }
                    }
                }

            }

            if(beta < 0.0)
            {
                beta = -std::sqrt(-beta);
            }
            else
            {
                beta = std::sqrt(beta);
            }

            TSparseSpaceType::InplaceMult(y, 1.0/beta);

            if (mEchoLevel > 1)
            {
                std::cout << i << " \t " << beta << " \t " << ro << " \t " << convergence_norm  << " \t\t " <<  min_shift_value << " \t " << max_shift_value << std::endl;
            }

            if(convergence_norm < tolerance)
            {
                break;
            }

            old_ro = ro;
        }

        if (mEchoLevel > 0)
        {
            KRATOS_WATCH(ro);
            KRATOS_WATCH(y);
        }

        Eigenvalues[0] = ro;

        if((Eigenvectors.size1() < 1) || (Eigenvectors.size2() < size))
        {
            Eigenvectors.resize(1,size);
        }

        for(SizeType i = 0 ; i < size ; i++)
        {
            Eigenvectors(0,i) = x[i] / beta;
        }
    }

    /**
     * This method returns directly the first eigen value obtained
     * @param K: The stiffness matrix
     * @param M: The mass matrix
     * @return The first eigenvalue
     */
    double GetEigenValue(
        SparseMatrixType& K,
        SparseMatrixType& M
        )
    {
        DenseVectorType eigen_values;
        DenseMatrixType eigen_vectors;
        
        Solve(K, M, eigen_values, eigen_vectors);
        
        return eigen_values[0];
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
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Power iteration eigenvalue solver with " << BaseType::GetPreconditioner()->Info();
        return  buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
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































