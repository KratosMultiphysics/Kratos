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
//

#if !defined(KRATOS_POWER_ITERATION_HIGHEST_EIGENVALUE_SOLVER_H_INCLUDED )
#define  KRATOS_POWER_ITERATION_HIGHEST_EIGENVALUE_SOLVER_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <numeric>
#include <vector>

// External includes

// Project includes
#include "spaces/ublas_space.h"
#include "includes/ublas_interface.h"
#include "processes/process.h"
#include "includes/define.h"
#include "linear_solvers/iterative_solver.h"
#include "utilities/random_initializer_utility.h"

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

/// This class uses the inverted power iteration method to obtain the lowest eigenvalue of a system
/** Basically that
*/
template<class TSparseSpaceType, class TDenseSpaceType, class TLinearSolverType,
         class TPreconditionerType = Preconditioner<TSparseSpaceType, TDenseSpaceType>,
         class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class PowerIterationHighestEigenvalueSolver : public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PowerIterationHighestEigenvalueSolver
    KRATOS_CLASS_POINTER_DEFINITION(PowerIterationHighestEigenvalueSolver);

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
    PowerIterationHighestEigenvalueSolver() {}

    PowerIterationHighestEigenvalueSolver(
        double MaxTolerance,
        unsigned int MaxIterationNumber,
        unsigned int RequiredEigenvalueNumber,
        typename TLinearSolverType::Pointer pLinearSolver
    ): BaseType(MaxTolerance, MaxIterationNumber),   
       mRequiredEigenvalueNumber(RequiredEigenvalueNumber),
       mpLinearSolver(pLinearSolver)
    {

    }

    PowerIterationHighestEigenvalueSolver(
        Parameters ThisParameters,
        typename TLinearSolverType::Pointer pLinearSolver
        ): mpLinearSolver(pLinearSolver)
    {
        Parameters DefaultParameters = Parameters(R"(
        {
            "solver_type"             : "power_iteration_highest_eigenvalue_solver",
            "max_iteration"           : 10000,
            "tolerance"               : 1e-8,
            "required_eigen_number"   : 1,
            "shifting_convergence"    : 0.25,
            "verbosity"               : 1,
            "linear_solver_settings"  : {}
        })" );

        ThisParameters.ValidateAndAssignDefaults(DefaultParameters);

        mRequiredEigenvalueNumber = ThisParameters["required_eigen_number"].GetInt();
        mEchoLevel = ThisParameters["verbosity"].GetInt();
        BaseType::SetTolerance( ThisParameters["tolerance"].GetDouble() );
        BaseType::SetMaxIterationsNumber( ThisParameters["max_iteration"].GetInt() );
    }

    /// Copy constructor.
    PowerIterationHighestEigenvalueSolver(const PowerIterationHighestEigenvalueSolver& Other) : BaseType(Other)
    {

    }


    /// Destructor.
    ~PowerIterationHighestEigenvalueSolver() override
    {

    }


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    PowerIterationHighestEigenvalueSolver& operator=(const PowerIterationHighestEigenvalueSolver& Other)
    {
        BaseType::operator=(Other);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * The power iteration algorithm
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
        using boost::numeric::ublas::trans;

        const SizeType size = K.size1();
        const SizeType max_iteration = BaseType::GetMaxIterationsNumber();
        const double tolerance = BaseType::GetTolerance();

        VectorType x = ZeroVector(size);
        VectorType y = ZeroVector(size);

        RandomInitializeUtility<double>::RandomInitialize(K, y);

        if(Eigenvalues.size() < 1)
        {
            Eigenvalues.resize(1, 0.0);
        }

        // Starting with first step
        double ro = 0.0;
        double old_ro = Eigenvalues[0];
        VectorType y_old = ZeroVector(size);

        if (mEchoLevel > 1)
        {
            std::cout << "Iteration ro \t\t convergence norm" << std::endl;
        }

        for(SizeType i = 0 ; i < max_iteration ; i++)
        {
            // x = K*y
            TSparseSpaceType::Mult(K, y, x);
            
            // y = M*x
            TSparseSpaceType::Mult(M, x, y);
            
            ro = static_cast<double>(*boost::max_element(y));
            
            TSparseSpaceType::InplaceMult(y, 1.0/ro);

            if(ro == 0.0)
            {
                KRATOS_ERROR << "Perpendicular eigenvector to M" << std::endl;
            }

            const double convergence_ro = std::abs((ro - old_ro) / ro);
            const double convergence_norm = TSparseSpaceType::TwoNorm(y - y_old)/TSparseSpaceType::TwoNorm(y);

            if (mEchoLevel > 1)
            {
                std::cout << "Iteration: " << i << "\tro: " << ro << " \tConvergence norm: " << convergence_norm << " \tConvergence ro: " << convergence_ro << std::endl;
            }
            
            if(convergence_norm < tolerance || convergence_ro < tolerance)
            {
                break;
            }

            old_ro = ro;
            TSparseSpaceType::Assign(y_old, 1.0, y);
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
            Eigenvectors(0,i) = y[i];
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

}; // Class PowerIterationHighestEigenvalueSolver

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
                                  PowerIterationHighestEigenvalueSolver<TSparseSpaceType, TDenseSpaceType,
                                  TPreconditionerType, TReordererType>& rThis)
{
    return IStream;
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType,
         class TPreconditionerType,
         class TReordererType>
inline std::ostream& operator << (std::ostream& OStream,
                                  const PowerIterationHighestEigenvalueSolver<TSparseSpaceType, TDenseSpaceType,
                                  TPreconditionerType, TReordererType>& rThis)
{
    rThis.PrintInfo(OStream);
    OStream << std::endl;
    rThis.PrintData(OStream);

    return OStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_POWER_ITERATION_HIGHEST_EIGENVALUE_SOLVER_H_INCLUDED defined
