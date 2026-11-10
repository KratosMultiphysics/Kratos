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
#include "includes/variables.h"
#include "future/linear_solvers/eigen_sparse_lu_solver.h"
#include "future/linear_solvers/linear_solver.h"
#include "future/linear_solvers/skyline_lu_factorization_solver.h"
#include "linear_solvers/rigid_body_mode_pinning_kernel.h"

namespace Kratos::Future
{
///@name Kratos Classes
///@{

/**
 * @class RigidBodyModePinningSolver
 * @ingroup KratosCore
 * @brief Linear solver wrapper that removes the (effective) rigid body modes of the
 *        unconstrained components of the system by pinning a minimal set of DOFs.
 * @details The algorithm is implemented by RigidBodyModePinningKernel (see its
 *          documentation); this class only wires it into the future linear solver
 *          interface: the pinning is applied in InitializeSolutionStep (before the
 *          wrapped solver factorizes) and the rigid body content is filtered out of
 *          the solution in PerformSolutionStep. The model part is taken from the
 *          additional data of the linear system (PrepareAdditionalData).
 * @details The Parameters constructor builds the wrapped solver from the nested
 *          "inner_solver_settings" parameter (default {"solver_type":
 *          "skyline_lu_factorization"}), using the same settings spelling as the
 *          current-interface solver (which delegates the nested settings to the
 *          LinearSolverFactory). Since the future interface has no solver factory yet,
 *          the concrete solver is selected here from the nested "solver_type"
 *          (currently "skyline_lu_factorization" and "eigen_sparse_lu"), and the nested
 *          settings are passed to it as is, so any inner-solver parameter may be set
 *          inside it.
 *          With "echo_level" > 0 the number of connected components and the number of pinned
 *          DOFs is printed; with "echo_level" > 1 the pinned DOFs (dof id and node id) are
 *          listed.
 * @note The rows of A are modified in place (as ScalingSolver does): the matrix must be
 *       re-assembled before the next call to InitializeSolutionStep.
 * @tparam TLinearAlgebra The version of the linear algebra to be used
 */
template<class TLinearAlgebra>
class RigidBodyModePinningSolver
    : public LinearSolver<TLinearAlgebra>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RigidBodyModePinningSolver
    KRATOS_CLASS_POINTER_DEFINITION(RigidBodyModePinningSolver);

    /// Definition of the base type
    using BaseType = LinearSolver<TLinearAlgebra>;

    /// The definition of the linear system type
    using LinearSystemType = typename BaseType::LinearSystemType;

    /// The definition of the spaces (sparse matrix)
    using MatrixType = typename TLinearAlgebra::MatrixType;

    /// The definition of the spaces (vector)
    using VectorType = typename TLinearAlgebra::VectorType;

    /// The index type definition to be consistent
    using IndexType = typename BaseType::IndexType;

    /// Dense vector tag type definition
    using DenseVectorTag = typename LinearSystemTags::DenseVectorTag;

    /// Sparse matrix tag type definition
    using SparseMatrixTag = typename LinearSystemTags::SparseMatrixTag;

    /// The pinning algorithm (shared with the current interface wrapper)
    using KernelType = Kratos::RigidBodyModePinningKernel<MatrixType, VectorType>;

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
        this->mEchoLevel = 0;
    }

    /**
     * @brief Constructor from parameters
     * @details The wrapped solver is selected through the nested "inner_solver_settings"
     *          parameter (default {"solver_type": "skyline_lu_factorization"}); the nested
     *          settings are passed to it as is, so any inner-solver parameter may be set
     *          inside it.
     * @param Settings The configuration parameters
     */
    RigidBodyModePinningSolver(Parameters Settings = Parameters(R"({})"))
        : BaseType()
    {
        // Validate and assign default parameters
        Settings.ValidateAndAssignDefaults(GetDefaultParameters());

        // Assign the validated settings
        this->AssignSettings(Settings);

        mKernel.SetTolerance(Settings["tolerance"].GetDouble());

        // Build the wrapped solver from the nested settings, with the same spelling as the
        // current-interface solver: "inner_solver_settings" already carries the inner
        // "solver_type" and any inner-solver parameter, so it is passed to the concrete
        // solver as is (the future interface has no solver factory yet). The linear system
        // tags and solve flags of the validated outer settings are propagated to the inner
        // solver, unless they were explicitly set inside "inner_solver_settings"
        Parameters inner_settings = Settings["inner_solver_settings"];
        for (const std::string& r_key : {"dx_tag", "rhs_tag", "lhs_tag", "multiple_solve", "echo_level"}) {
            if (!inner_settings.Has(r_key)) {
                inner_settings.AddValue(r_key, Settings.GetValue(r_key));
            }
        }
        const std::string inner_solver_type = inner_settings["solver_type"].GetString();
        if (inner_solver_type == "skyline_lu_factorization") {
            mpLinearSolver = Kratos::make_shared<SkylineLUFactorizationSolver<TLinearAlgebra>>(inner_settings);
        } else if (inner_solver_type == "eigen_sparse_lu") {
            mpLinearSolver = Kratos::make_shared<EigenSparseLUSolver<TLinearAlgebra>>(inner_settings);
        } else {
            KRATOS_ERROR << "RigidBodyModePinningSolver: unknown \"solver_type\" in \"inner_solver_settings\": \""
                         << inner_solver_type
                         << "\" (expected \"skyline_lu_factorization\" or \"eigen_sparse_lu\")." << std::endl;
        }
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
    void Initialize(LinearSystemType& rLinearSystem) override
    {
        mpLinearSolver->Initialize(rLinearSystem);
    }

    /**
     * @brief Pins the rigid body modes of the unconstrained components in place and
     *        delegates the initialization of the solution step to the wrapped solver.
     * @details The pinning is applied before delegating, so a wrapped direct solver
     *          factorizes the already-pinned matrix in its own InitializeSolutionStep.
     * @param rLinearSystem The linear system to be solved (the matrix and the vectors are
     *                      modified in place: pinned rows/columns are zeroed off-diagonal,
     *                      pinned entries of the rhs and of the initial guess are set to zero)
     */
    void InitializeSolutionStep(LinearSystemType& rLinearSystem) override
    {
        if (rLinearSystem.IsConsistent(SparseMatrixTag::LHS, DenseVectorTag::RHS, DenseVectorTag::Dx)) {
            auto& r_lin_op = *rLinearSystem.pGetLinearOperator(SparseMatrixTag::LHS);
            KRATOS_ERROR_IF(r_lin_op.IsMatrixFree()) << "RigidBodyModePinningSolver cannot be used with matrix-free linear operators." << std::endl;
            auto& rA = r_lin_op.GetMatrix();
            auto& rX = *rLinearSystem.pGetVector(DenseVectorTag::Dx);
            auto& rB = *rLinearSystem.pGetVector(DenseVectorTag::RHS);

            mKernel.PinRigidBodyModes(rA, rX, rB, this->mEchoLevel);
        }

        mpLinearSolver->InitializeSolutionStep(rLinearSystem);
    }

    /**
     * @brief Delegates the solution of the system to the wrapped solver and removes the
     *        rigid body content of the solution.
     * @details After a successful solve, the rigid body content of the unconstrained
     *          components is projected out of the solution, so the returned solution is
     *          the unique minimal-norm one, independent of the DOFs that were pinned.
     */
    bool PerformSolutionStep(LinearSystemType& rLinearSystem) override
    {
        const bool success = mpLinearSolver->PerformSolutionStep(rLinearSystem);
        if (success && mKernel.HasRigidBodyModes()) {
            mKernel.FilterRigidBodyContent(*rLinearSystem.pGetVector(DenseVectorTag::Dx));
        }
        return success;
    }

    /// Delegates the finalization of the solution step to the wrapped solver
    void FinalizeSolutionStep(LinearSystemType& rLinearSystem) override
    {
        mpLinearSolver->FinalizeSolutionStep(rLinearSystem);
    }

    /// Delegates the cleanup to the wrapped solver
    void Clear() override
    {
        mpLinearSolver->Clear();
    }

    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters(R"({
            "solver_type" : "rigid_body_mode_pinning",
            "inner_solver_settings" : {
                "solver_type" : "skyline_lu_factorization"
            },
            "tolerance" : 1.0e-6
        })");
        default_parameters.AddMissingParameters(BaseType::GetDefaultParameters());

        return default_parameters;
    }

    /**
     * @brief Builds the rigid body mode matrix Z (n x M) from the node coordinates and DOF types.
     * @details M = 3 in 2D and M = 6 in 3D (see RigidBodyModePinningKernel for the column layout).
     *          The model part is retrieved from the additional data of the linear system.
     */
    void PrepareAdditionalData(LinearSystemType& rLinearSystem) override
    {
        mpLinearSolver->PrepareAdditionalData(rLinearSystem);

        auto& r_model_part = *rLinearSystem.pGetModelPart();
        auto& r_A = *rLinearSystem.pGetMatrix(SparseMatrixTag::LHS);

        mKernel.BuildRigidBodyModes(r_A, r_model_part);
    }

    ///@}
    ///@name Access
    ///@{

    /// The rigid body modes need the spatial position of the nodes
    bool RequiresAdditionalData() const override
    {
        return true;
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
    KernelType mKernel;                        /// The pinning algorithm (shared with the current interface wrapper)

    ///@}

}; // Class RigidBodyModePinningSolver

///@}

}  // namespace Kratos::Future.
