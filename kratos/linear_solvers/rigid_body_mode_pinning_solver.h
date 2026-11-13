//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

#pragma once

// System includes
#include <sstream>

// Project includes
#include "factories/linear_solver_factory.h"
#include "linear_solvers/linear_solver.h"
#include "linear_solvers/rigid_body_mode_pinning_kernel.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @class RigidBodyModePinningSolver
 * @ingroup KratosCore
 * @brief Linear solver wrapper that removes the (effective) rigid body modes of the
 *        unconstrained components of the system by pinning a minimal set of DOFs.
 * @details The algorithm is implemented by RigidBodyModePinningKernel (see its
 *          documentation); this class only wires it into the current linear solver
 *          interface: the pinning is applied in InitializeSolutionStep (before the
 *          wrapped solver factorizes) and the rigid body content is filtered out of
 *          the solution in PerformSolutionStep.
 * @details The Parameters constructor builds the wrapped solver from the nested
 *          "inner_solver_settings" parameter through the LinearSolverFactory
 *          (default {"solver_type": "skyline_lu_factorization"}); any solver registered in
 *          the factory may be used, and any inner-solver parameter may be set inside it.
 *          With "echo_level" > 0 the number of connected components and the number of pinned
 *          DOFs is printed; with "echo_level" > 1 the pinned DOFs (dof id and node id) are
 *          listed.
 * @note The rows of A are modified in place (as ScalingSolver does): the matrix must be
 *       re-assembled before the next call to InitializeSolutionStep.
 * @tparam TSparseSpaceType The sparse space definition
 * @tparam TDenseSpaceType The dense space definition
 * @tparam TReordererType The reorder considered
 */
template<class TSparseSpaceType, class TDenseSpaceType,
         class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class RigidBodyModePinningSolver
    : public LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RigidBodyModePinningSolver
    KRATOS_CLASS_POINTER_DEFINITION(RigidBodyModePinningSolver);

    /// Definition of the base type
    using BaseType = LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>;

    /// The definition of the spaces (sparse matrix)
    using SparseMatrixType = typename TSparseSpaceType::MatrixType;

    /// The definition of the spaces (vector)
    using VectorType = typename TSparseSpaceType::VectorType;

    /// The definition of the linear solver factory type
    using LinearSolverFactoryType = LinearSolverFactory<TSparseSpaceType, TDenseSpaceType>;

    /// The index type definition to be consistent
    using IndexType = typename BaseType::IndexType;

    /// The pinning algorithm (shared with the future interface wrapper)
    using KernelType = RigidBodyModePinningKernel<SparseMatrixType, VectorType>;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor without parameters
     * @param pLinearSolver The linear solver to be wrapped
     * @param Tolerance Relative tolerance used to detect the unconstrained rigid body modes
     */
    RigidBodyModePinningSolver(
        typename BaseType::Pointer pLinearSolver,
        const double Tolerance = 1.0e-6
        ) : BaseType(),
            mpLinearSolver(pLinearSolver),
            mKernel(Tolerance)
    {
        mEchoLevel = 0;
    }

    /**
     * @brief Constructor from parameters
     * @details The wrapped solver is selected through the nested "inner_solver_settings"
     *          parameter (default {"solver_type": "skyline_lu_factorization"}) and built with
     *          the LinearSolverFactory; any inner-solver parameter may be set inside it.
     * @param ThisParameters The configuration parameters
     */
    RigidBodyModePinningSolver(Parameters ThisParameters = Parameters(R"({})"))
        : BaseType()
    {
        KRATOS_TRY

        // Validate and assign default parameters
        Parameters default_parameters(R"(
        {
            "solver_type"           : "rigid_body_mode_pinning",
            "inner_solver_settings" : {
                "solver_type" : "skyline_lu_factorization"
            },
            "tolerance"             : 1.0e-6,
            "echo_level"            : 0
        }
        )");

        ThisParameters.ValidateAndAssignDefaults(default_parameters);

        mTolerance = ThisParameters["tolerance"].GetDouble();
        mKernel.SetTolerance(mTolerance);
        mEchoLevel = ThisParameters["echo_level"].GetInt();

        // Build the wrapped solver from the nested settings, as for the other factory-based
        // wrapper solvers (e.g. MonotonicityPreservingSolver): "inner_solver_settings" already
        // carries the inner "solver_type" and any inner-solver parameter, so it is delegated
        // to the factory as is
        mpLinearSolver = LinearSolverFactoryType().Create(ThisParameters["inner_solver_settings"]);

        KRATOS_CATCH("")
    }

    /// Copy constructor.
    RigidBodyModePinningSolver(const RigidBodyModePinningSolver& Other) = delete;

    /// Destructor.
    ~RigidBodyModePinningSolver() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    RigidBodyModePinningSolver& operator=(const RigidBodyModePinningSolver& Other) = delete;

    ///@}
    ///@name Operations
    ///@{

    /// Delegates the initialization to the wrapped solver
    void Initialize(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        mpLinearSolver->Initialize(rA, rX, rB);
    }

    /**
     * @brief Pins the rigid body modes of the unconstrained components in place and
     *        delegates the initialization of the solution step to the wrapped solver.
     * @details The pinning is applied before delegating, so a wrapped direct solver
     *          factorizes the already-pinned matrix in its own InitializeSolutionStep.
     * @param rA. System matrix (modified in place: pinned rows/columns are zeroed off-diagonal)
     * @param rX. Solution vector (pinned entries are set to zero)
     * @param rB. Right hand side vector (pinned entries are set to zero)
     */
    void InitializeSolutionStep(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        if (this->IsConsistent(rA, rX, rB)) {
            mKernel.PinRigidBodyModes(rA, rX, rB, mEchoLevel);
        }

        mpLinearSolver->InitializeSolutionStep(rA, rX, rB);
    }

    /**
     * @brief Delegates the solution of the system to the wrapped solver and removes the
     *        rigid body content of the solution.
     * @details After a successful solve, the rigid body content of the unconstrained
     *          components is projected out of the solution, so the returned solution is
     *          the unique minimal-norm one, independent of the DOFs that were pinned.
     */
    bool PerformSolutionStep(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        const bool success = mpLinearSolver->PerformSolutionStep(rA, rX, rB);
        if (success && mKernel.HasRigidBodyModes()) {
            mKernel.FilterRigidBodyContent(rX);
        }
        return success;
    }

    /// Delegates the finalization of the solution step to the wrapped solver
    void FinalizeSolutionStep(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        mpLinearSolver->FinalizeSolutionStep(rA, rX, rB);
    }

    /// Delegates the cleanup to the wrapped solver
    void Clear() override
    {
        mpLinearSolver->Clear();
    }

    /**
     * @brief Builds the rigid body mode matrix Z (n x M) from the node coordinates and DOF types.
     * @details M = 3 in 2D and M = 6 in 3D (see RigidBodyModePinningKernel for the column layout).
     */
    void ProvideAdditionalData(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB,
        typename ModelPart::DofsArrayType& rdof_set,
        ModelPart& r_model_part
    ) override
    {
        mpLinearSolver->ProvideAdditionalData(rA, rX, rB, rdof_set, r_model_part);

        mKernel.BuildRigidBodyModes(rA, r_model_part);
    }

    /// The rigid body modes need the spatial position of the nodes
    bool AdditionalPhysicalDataIsNeeded() override
    {
        return true;
    }

    ///@}
    ///@name Access
    ///@{

    void SetTolerance(double NewTolerance) override
    {
        mTolerance = NewTolerance;
        mKernel.SetTolerance(NewTolerance);
    }

    double GetTolerance() override
    {
        return mTolerance;
    }

    ///@}
    ///@name Inquiry
    ///@{

    IndexType GetIterationsNumber() override
    {
        return mpLinearSolver->GetIterationsNumber();
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Rigid body mode pinning solver. Uses internally the following linear solver " << mpLinearSolver->Info();
        return buffer.str();
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

    typename BaseType::Pointer mpLinearSolver; /// Pointer to the internal linear solver
    double mTolerance = 1.0e-6;                /// Relative tolerance to detect unconstrained rigid body modes
    int mEchoLevel = 0;                        /// Echo level (0: silent, 1: pinning summary, 2: pinned dof list)
    KernelType mKernel;                        /// The pinning algorithm (shared with the future interface wrapper)

    ///@}

}; // Class RigidBodyModePinningSolver

///@}

}  // namespace Kratos.
