//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                    
//

#if !defined(KRATOS_POWER_ITERATION_EIGENVALUE_SOLVERR_H_INCLUDED )
#define  KRATOS_POWER_ITERATION_EIGENVALUE_SOLVERR_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <numeric>
#include <vector>


// External includes


// Project includes
#include "includes/define.h"
#include "linear_solvers/iterative_solver.h"


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
class PowerIterationEigenvalueSolver : public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PowerIterationEigenvalueSolver
    KRATOS_CLASS_POINTER_DEFINITION(PowerIterationEigenvalueSolver);

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
    PowerIterationEigenvalueSolver() {}

    PowerIterationEigenvalueSolver(double NewMaxTolerance, unsigned int NewMaxIterationsNumber,
                                   unsigned int NewRequiredEigenvalueNumber, typename TLinearSolverType::Pointer pLinearSolver)
        : BaseType(NewMaxTolerance, NewMaxIterationsNumber), mRequiredEigenvalueNumber(NewRequiredEigenvalueNumber), mpLinearSolver(pLinearSolver) {}

    /*       PowerIterationEigenvalueSolver(double NewMaxTolerance, unsigned int NewMaxIterationsNumber, typename TPreconditionerType::Pointer pNewPreconditioner) :  */
    /*       BaseType(NewMaxTolerance, NewMaxIterationsNumber, pNewPreconditioner){} */

    /// Copy constructor.
    PowerIterationEigenvalueSolver(const PowerIterationEigenvalueSolver& Other) : BaseType(Other) {}


    /// Destructor.
    ~PowerIterationEigenvalueSolver() override {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    PowerIterationEigenvalueSolver& operator=(const PowerIterationEigenvalueSolver& Other)
    {
        BaseType::operator=(Other);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    static void RandomInitialize(DenseVectorType& R)
    {
        for(SizeType i = 0 ; i < R.size() ; i++)
            R[i] = 1.00; //rand();

        R /= norm_2(R);
    }


    // The power iteration algorithm
    void Solve(SparseMatrixType& K,
               SparseMatrixType& M,
               DenseVectorType& Eigenvalues,
               DenseMatrixType& Eigenvectors) override
    {

        using boost::numeric::ublas::trans;

        SizeType size = K.size1();
        SizeType max_iteration = BaseType::GetMaxIterationsNumber();
        double tolerance = BaseType::GetTolerance();

        VectorType x = ZeroVector(size);
        VectorType y = ZeroVector(size);

        RandomInitialize(y);

        if(Eigenvalues.size() < 1)
            Eigenvalues.resize(1,0.00);


        // Starting with first step
        double beta = 0.00;
        double ro = 0.00;
        double old_ro = Eigenvalues[0];
        std::cout << "iteration    beta \t\t ro \t\t convergence norm" << std::endl;
        for(SizeType i = 0 ; i < max_iteration ; i++)
        {
            //K*x = y
            mpLinearSolver->Solve(K,x,y);

            ro = inner_prod(y,x);

            //y = M*x
            noalias(y) = prod(M,x);

            beta = inner_prod(x, y);
            if(beta <= 0.00)
                KRATOS_THROW_ERROR(std::invalid_argument, "M is not Positive-definite", "");

            ro = ro / beta;
            beta = sqrt(beta);

            double inverse_of_beta = 1.00 / beta;

            y *= inverse_of_beta;

            if(ro == 0.00)
                KRATOS_THROW_ERROR(std::runtime_error, "Perpendicular eigenvector to M", "");

            double convergence_norm = fabs((ro - old_ro) / ro);

            std::cout << i << " \t " << beta << " \t " << ro << " \t " << convergence_norm << std::endl;
            //std::cout << "i = " << i << ": beta = " << beta << ", ro = " << ro << ", convergence norm = " << convergence_norm << std::endl;

            if(convergence_norm < tolerance)
                break;

            old_ro = ro;



        }

        KRATOS_WATCH(ro);
//KRATOS_WATCH(y);

        Eigenvalues[0] = ro;

        if((Eigenvectors.size1() < 1) || (Eigenvectors.size2() < size))
            Eigenvectors.resize(1,size);

        for(SizeType i = 0 ; i < size ; i++)
            Eigenvectors(0,i) = y[i];
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

    typename TLinearSolverType::Pointer mpLinearSolver;

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

}; // Class PowerIterationEigenvalueSolver

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
                                  PowerIterationEigenvalueSolver<TSparseSpaceType, TDenseSpaceType,
                                  TPreconditionerType, TReordererType>& rThis)
{
    return IStream;
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType,
         class TPreconditionerType,
         class TReordererType>
inline std::ostream& operator << (std::ostream& OStream,
                                  const PowerIterationEigenvalueSolver<TSparseSpaceType, TDenseSpaceType,
                                  TPreconditionerType, TReordererType>& rThis)
{
    rThis.PrintInfo(OStream);
    OStream << std::endl;
    rThis.PrintData(OStream);

    return OStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_POWER_ITERATION_EIGENVALUE_SOLVERR_H_INCLUDED defined 































