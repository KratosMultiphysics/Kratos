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

#pragma once

// System includes

// External includes

// Project includes
#include "factories/linear_solver_factory.h"
#include "linear_solvers/linear_solver.h"
#include "utilities/parallel_utilities.h"

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
 * @brief A linear solver wrapper that ensures monotonicity by modifying the system matrix and right-hand side.
 * @details This class wraps another linear solver (the "inner solver").
 * Before calling the inner solver's Solve() method, this class pre-processes the system matrix 'A' and the right-hand side 'B' within the ProvideAdditionalData() method.
 * The pre-processing ensures that the matrix 'A' has non-positive off-diagonal entries.
 * It iterates over all entries A(i,j) of the matrix. If a positive off-diagonal
 * entry is found (A(i,j) > 0 with i != j), this positive value is:
 * 1. Subtracted from the off-diagonal entries A(i,j) and A(j,i).
 * 2. Added to the diagonal entries A(i,i) and A(j,j).
 * To keep the system solution unchanged, the right-hand side vector 'B' is also
 * modified based on the current solution values (dofs_values) associated with
 * the degrees of freedom:
 * B[i] += A(i,j)_original * (dofs_values[j] - dofs_values[i])
 * B[j] += A(i,j)_original * (dofs_values[i] - dofs_values[j])
 * This transformation enforces a condition required for monotonicity (Discrete Maximum Principle) often needed in transport or convection-diffusion problems, without altering the final solution of the linear system. The modified system is then passed to the inner solver.
 * The standard, linear conforming FEM, while robust, may not guarantee the DMP on general or distorted meshes for anisotropic diffusion. The DMP is a desirable qualitative property stating that if the source term is non-negative, the solution should attain its minimum on the boundary and be everywhere non-negative.
 * A key requirement for a scheme to satisfy the DMP is for the LHS matrix A to be an M-matrix. An M-matrix has non-positive off-diagonal entries and a non-negative inverse. For standard FEM, the off-diagonal entries of the LHS matrix are generally not guaranteed to be non-positive, leading to a loss of monotonicity.
 * @author Daniel Diez
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
    using BaseType = LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>;

    /// The definition of the spaces (sparse matrix)
    using SparseMatrixType = typename TSparseSpaceType::MatrixType;

    /// The definition of the spaces (vector)
    using VectorType = typename TSparseSpaceType::VectorType;

    /// The definition of the spaces (dense matrix)
    using DenseMatrixType = typename TDenseSpaceType::MatrixType;

    /// The definition of the linear solver factory type
    using LinearSolverFactoryType = LinearSolverFactory<TSparseSpaceType,TDenseSpaceType>;

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

        Parameters default_parameters( R"(
        {
            "solver_type"                    : "monotonicity_preserving",
            "inner_solver_settings"          : {
                "preconditioner_type"            : "amg",
                "solver_type"                    : "AMGCL",
                "smoother_type"                  : "ilu0",
                "krylov_type"                    : "lgmres",
                "coarsening_type"                : "aggregation",
                "max_iteration"                  : 100,
                "provide_coordinates"            : false,
                "gmres_krylov_space_dimension"   : 100,
                "verbosity"                      : 1,
                "tolerance"                      : 1e-6,
                "scaling"                        : false,
                "block_size"                     : 1,
                "use_block_matrices_if_possible" : true,
                "coarse_enough"                  : 1000,
                "max_levels"                     : -1,
                "pre_sweeps"                     : 1,
                "post_sweeps"                    : 1,
                "use_gpgpu"                      : false
            }

        }  )" );

        // Now validate agains defaults -- this also ensures no type mismatch
        ThisParameters.ValidateAndAssignDefaults(default_parameters);

        mpLinearSolver = LinearSolverFactoryType().Create(ThisParameters["inner_solver_settings"]);

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

    /** 
     * @brief Some solvers may require a minimum degree of knowledge of the structure of the matrix.
     * @details To make an example when solving a mixed u-p problem, it is important to identify the row associated to v and p. another example is the automatic prescription of rotation null-space for smoothed-aggregation solvers which require knowledge on the spatial position of the nodes associated to a given dof.
     * @return This function tells if the solver requires such data, always true
     */
    bool AdditionalPhysicalDataIsNeeded() override
    {
        return true;
    }

    /** 
     * @brief Some solvers may require a minimum degree of knowledge of the structure of the matrix. 
     * @details To make an example when solving a mixed u-p problem, it is important to identify the row associated to v and p. another example is the automatic prescription of rotation null-space for smoothed-aggregation solvers which require knowledge on the spatial position of the nodes associated to a given dof. This function is the place to eventually provide such data
     * @param rA. System matrix
     * @param rX. Solution vector. it's also the initial guess for iterative linear solvers.
     * @param rB. Right hand side vector.
     * @param rDofSet The degrees of freedom set.
     * @param rModelPart The modelpart to be considered.
     */
    void ProvideAdditionalData(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB,
        typename ModelPart::DofsArrayType& rDofSet,
        ModelPart& rModelPart
    ) override
    {
        VectorType dofs_values(rDofSet.size());
        TSparseSpaceType::SetToZero(dofs_values);

        block_for_each(rDofSet, [&](Dof<double>& rDof){
            const std::size_t id = rDof.EquationId();
            dofs_values[id] = rDof.GetSolutionStepValue();
        });
        typename SparseMatrixType::value_type *values_vector = rA.value_data().begin();
        std::size_t *index1_vector = rA.index1_data().begin();
        std::size_t *index2_vector = rA.index2_data().begin();

        IndexPartition<std::size_t>(rA.size1()).for_each(
            [&](std::size_t i)
            {
                for (std::size_t k = index1_vector[i]; k < index1_vector[i + 1]; k++) {
                    const typename SparseMatrixType::value_type value = values_vector[k];
                    if (value > 0.0) {
                        const auto j = index2_vector[k];
                        if (j > i) {
                            rA(i,j) -= value;
                            rA(j,i) -= value;
                            // Values conflicting with other threads
                            auto& r_aii = rA(i,i).ref();
                            AtomicAdd(r_aii, value);
                            auto& r_ajj = rA(j,j).ref();
                            AtomicAdd(r_ajj, value);
                            auto& r_bi = rB[i];
                            AtomicAdd(r_bi, value*dofs_values[j] - value*dofs_values[i]);
                            auto& r_bj = rB[j];
                            AtomicAdd(r_bj, value*dofs_values[i] - value*dofs_values[j]);
                        }
                    }
                }
            }
        );

        // Invoque ProvideAdditionalData from base solver if required
        if (mpLinearSolver->AdditionalPhysicalDataIsNeeded()) {
            mpLinearSolver->ProvideAdditionalData(rA,rX,rB,rDofSet,rModelPart);
        }
    }

    /**
     * @brief Function called at the beginning of the solution step
     * @details Invoques the InitializeSolutionStep from the base linear solver
     * @param rA. System matrix
     * @param rX. Solution vector. it's also the initial guess for iterative linear solvers.
     * @param rB. Right hand side vector.
     */
    void InitializeSolutionStep(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB
        ) override
    {
        mpLinearSolver->InitializeSolutionStep(rA,rX,rB);
    }

    /** 
     * @brief This function is designed to be called at the end of the solve step.
     * @details For example this is the place to remove any data that we do not want to save for later
     * @param rA. System matrix
     * @param rX. Solution vector. it's also the initial guess for iterative linear solvers.
     * @param rB. Right hand side vector.
     */
    void FinalizeSolutionStep(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB
        ) override
    {
        mpLinearSolver->FinalizeSolutionStep(rA,rX,rB);
    }

    /** 
     * @brief This function is designed to clean up all internal data in the solver.
     * @details Clear is designed to leave the solver object as if newly created.
     * After a clear a new Initialize is needed
     */
    void Clear() override
    {
        mpLinearSolver->Clear();
    }

    /**
     * @brief Solves the system of equations
     * @details Invoques the Solve method from the base linear solver 
     * @param rA. System matrix
     * @param rX. Solution vector. it's also the initial guess for iterative linear solvers.
     * @param rB. Right hand side vector.
     */
    bool Solve(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB
        ) override
    {
        return mpLinearSolver->Solve(rA,rX,rB);
    }

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
private:
    ///@name Member Variables
    ///@{

    typename LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>::Pointer mpLinearSolver;

    ///@}
}; // Class MonotonicityPreservingSolver

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