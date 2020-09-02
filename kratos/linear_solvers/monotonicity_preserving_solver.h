//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Daniel Diez
//

#if !defined(KRATOS_MONOTONICITY_PRESRVING_SOLVER_H_INCLUDED )
#define  KRATOS_MONOTONICITY_PRESRVING_SOLVER_H_INCLUDED

// System includes
#include <cmath>
#include <complex>

// External includes

// Project includes
#include "includes/define.h"
#include "factories/linear_solver_factory.h"
#include "linear_solvers/linear_solver.h"

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
 * @class MonotonicityPreservingSolver
 * @ingroup KratosCore
 * @brief
 * @details
 * @author Daniel Diex
 * @tparam TSparseSpaceType The sparse space definition
 * @tparam TDenseSpaceType The dense space definition
 * @tparam TReordererType The reorder considered
 */
template<class TSparseSpaceType, class TDenseSpaceType,
         class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class MonotonicityPreservingSolver
    : public LinearSolver<TSparseSpaceType, TDenseSpaceType,  TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MonotonicityPreservingSolver
    KRATOS_CLASS_POINTER_DEFINITION(MonotonicityPreservingSolver);

    /// Definition of the base type
    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    /// The definition of the spaces (sparse matrix)
    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    /// The definition of the spaces (vector)
    typedef typename TSparseSpaceType::VectorType VectorType;

    /// The definition of the spaces (dense matrix)
    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    /// The definition of the linear solver factory type
    typedef LinearSolverFactory<TSparseSpaceType,TDenseSpaceType> LinearSolverFactoryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MonotonicityPreservingSolver()
    {
    }

    /**
     * @brief Constructor without parameters
     * @param pLinearSolver The linear solver to be scaled
     */
    MonotonicityPreservingSolver(
        typename BaseType::Pointer pLinearSolver
        ) : BaseType (),
            mpLinearSolver(pLinearSolver)
    {
    }

    /**
     * @brief Constructor with parameters
     * @param ThisParameters The configuration parameters of the linear solver
     */
    MonotonicityPreservingSolver(Parameters ThisParameters)
        : BaseType ()
    {
        KRATOS_TRY

        KRATOS_ERROR_IF_NOT(ThisParameters.Has("solver_type")) << "Solver_type must be specified to construct the MonotonicityPreservingSolver" << std::endl;

        mpLinearSolver = LinearSolverFactoryType().Create(ThisParameters);

        KRATOS_CATCH("")
    }

    /// Copy constructor.
    MonotonicityPreservingSolver(const MonotonicityPreservingSolver& Other) : BaseType(Other) {}


    /// Destructor.
    ~MonotonicityPreservingSolver() override {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    MonotonicityPreservingSolver& operator=(const MonotonicityPreservingSolver& Other)
    {
        BaseType::operator=(Other);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{
    /** Some solvers may require a minimum degree of knowledge of the structure of the matrix. To make an example
    * when solving a mixed u-p problem, it is important to identify the row associated to v and p.
    * another example is the automatic prescription of rotation null-space for smoothed-aggregation solvers
    * which require knowledge on the spatial position of the nodes associated to a given dof.
    * This function tells if the solver requires such data
    */
    bool AdditionalPhysicalDataIsNeeded() override
    {
        return mpLinearSolver->AdditionalPhysicalDataIsNeeded();
    }

    /** Some solvers may require a minimum degree of knowledge of the structure of the matrix. To make an example
     * when solving a mixed u-p problem, it is important to identify the row associated to v and p.
     * another example is the automatic prescription of rotation null-space for smoothed-aggregation solvers
     * which require knowledge on the spatial position of the nodes associated to a given dof.
     * This function is the place to eventually provide such data
     */
    void ProvideAdditionalData(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB,
        typename ModelPart::DofsArrayType& rdof_set,
        ModelPart& r_model_part
    ) override
    {
        mpLinearSolver->ProvideAdditionalData(rA,rX,rB,rdof_set,r_model_part);
    }

    void InitializeSolutionStep (SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        mpLinearSolver->InitializeSolutionStep(rA,rX,rB);
    }

    /** This function is designed to be called at the end of the solve step.
     * for example this is the place to remove any data that we do not want to save for later
    @param rA. System matrix
    @param rX. Solution vector. it's also the initial guess for iterative linear solvers.
    @param rB. Right hand side vector.
    */
    void FinalizeSolutionStep (SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        mpLinearSolver->FinalizeSolutionStep(rA,rX,rB);
    }

    /** This function is designed to clean up all internal data in the solver.
     * Clear is designed to leave the solver object as if newly created.
     * After a clear a new Initialize is needed
     */
    void Clear() override
    {
        mpLinearSolver->Clear();
    }

    /**
    @param rA. System matrix
    @param rX. Solution vector. it's also the initial
    guess for iterative linear solvers.
    @param rB. Right hand side vector.
    */
    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        if(this->IsNotConsistent(rA, rX, rB))
            return false;

        double *values_vector = rA.value_data().begin();
        std::size_t *index1_vector = rA.index1_data().begin();
        std::size_t *index2_vector = rA.index2_data().begin();
        const std::size_t matrix_size = rA.size1();

        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(rA.size1()); i++) {
            for (std::size_t k = index1_vector[i]; k < index1_vector[i + 1]; k++) {
                const double value = values_vector[k];
                if (value > 0.0) {
                    const int j = index2_vector[k];
                    if (j > i) {
                        rA(i,i) += value;
                        rA(i,j) -= value;
                        rA(j,i) -= value;
                        rA(j,j) += value;
                        rB[i] += value*rX[j] - value*rX[i];
                        rB[j] += value*rX[i] - value*rX[j];
                    }
                }
            }
        }

        //solve the problem
        bool is_solved = mpLinearSolver->Solve(rA,rX,rB);

        return is_solved;
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
        buffer << "Composite Linear Solver. Uses internally the following linear solver " << mpLinearSolver->Info();
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
    typename LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>::Pointer mpLinearSolver;

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

}; // Class MonotonicityPreservingSolver

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
                                  MonotonicityPreservingSolver<TSparseSpaceType, TDenseSpaceType,
                                  TReordererType>& rThis)
{
    return IStream;
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType,
         class TPreconditionerType,
         class TReordererType>
inline std::ostream& operator << (std::ostream& OStream,
                                  const MonotonicityPreservingSolver<TSparseSpaceType, TDenseSpaceType,
                                  TReordererType>& rThis)
{
    rThis.PrintInfo(OStream);
    OStream << std::endl;
    rThis.PrintData(OStream);

    return OStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_MONOTONICITY_PRESRVING_SOLVER_H_INCLUDED  defined


