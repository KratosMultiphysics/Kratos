//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

// Project includes
#include "solving_strategies/builder_and_solvers/p_multigrid_builder_and_solver.hpp" // PMultigridBuilderAndSolver
#include "includes/variables.h" // NL_ITERATION_NUMBER
#include "spaces/ublas_space.h" // TUblasSparseSpace, TUblasDenseSpace
#include "linear_solvers/linear_solver.h" // LinearSolver
#include "includes/model_part.h" // ModelPart
#include "utilities/dof_utilities/block_build_dof_array_utility.h" // BlockBuildDofArrayUtility
#include "utilities/sparse_matrix_multiplication_utility.h" // SparseMatrixMultiplicationUtility
#include "utilities/atomic_utilities.h" // AtomicAdd
#include "utilities/proxies.h" // MakeProxy
#include "utilities/profiler.h" // KRATOS_PROFILE_SCOPE, KRATOS_PROFILE_SCOPE_MILLI
#include "solving_strategies/builder_and_solvers/p_multigrid_utilities.hpp" // MakePRestrictionOperator, MakeSparseTopology

// System includes
#include <algorithm> // std::lower_bound, std::sort, std::transform
#include <optional> // std::optional
#include <limits> // std::numeric_limits
#include <unordered_set> // std::unordered_set
#include <unordered_map> // std::unordered_map
#include <sstream> // std::stringstream
#include <numeric> // std::iota


namespace Kratos {


DiagonalScaling ParseDiagonalScaling(Parameters Settings)
{
    KRATOS_TRY
    const auto diagonal_scaling_strategy = Settings["diagonal_scaling"].Get<std::string>();
    if (diagonal_scaling_strategy == "none") {
        return DiagonalScaling::None;
    } else if (diagonal_scaling_strategy == "abs_max") {
        return DiagonalScaling::AbsMax;
    } else if (diagonal_scaling_strategy == "norm") {
        return DiagonalScaling::Norm;
    } else {
        KRATOS_ERROR << "unsupported setting for \"diagonal_scaling\": "
                        << diagonal_scaling_strategy << ". Options are:\n"
                        << "- \"none\"\n"
                        << "- \"abs_max\"\n"
                        << "- \"norm\"\n";
    }
    KRATOS_CATCH("")
}


enum class ConstraintImposition
{
    None                = 0,
    MasterSlave         = 1,
    Lagrange            = 2,
    AugmentedLagrange   = 3,
    Penalty             = 4
}; // enum class ConstraintImposition


template <class TSparse, class TDense>
class ConstraintAssembler : public DataValueContainer
{
public:
    ConstraintAssembler() noexcept
        : ConstraintAssembler(ConstraintImposition::MasterSlave)
    {}

    ConstraintAssembler(ConstraintImposition Method)
        : DataValueContainer()
    {
        std::string method_name;

        switch (Method) {
            case ConstraintImposition::MasterSlave:
                method_name = "master_slave";
                break;
            case ConstraintImposition::AugmentedLagrange:
                method_name = "augmented_lagrange";
                break;
            default:
                KRATOS_ERROR << "Unsupported constraint imposition method: " << (int)Method;
        } // switch Method

        this->SetValue(ConstraintAssembler::GetImpositionVariable(), method_name);
    }

    /// @details Define an overriding virtual destructor to ensure compile time errors
    //           if the base class' destructor turns non-virtual.
    virtual ~ConstraintAssembler() override = default;

    /// @brief Allocate memory for the constraint gap vector and relation matrix, and compute its topology.
    /// @details This function is responsible for large memory allocations, as well as computing
    ///          the sparsity pattern of the relation matrix. It must also modify the provided
    ///          left hand side matrix, solution vector, right hand side vector, and DoF list such that
    ///          these containers will not require reallocation during later stages of the solution process.
    /// @param rConstraints Constraint set of the related @ref ModelPart.
    /// @param rProcessInfo Current @ref ProcessInfo of the related @ref ModelPart.
    /// @param rLhs Unconstrained left hand side matrix' topology.
    /// @param rDofSet Unconstrained set of @ref Dof "DoFs".
    /// @note This function should be invoked @b after the unconstrained system is allocated, but @b before
    ///       it is assembled.
    virtual void Allocate(const ModelPart::MasterSlaveConstraintContainerType& rConstraints,
                          const ProcessInfo& rProcessInfo,
                          typename TSparse::MatrixType& rLhs,
                          typename TSparse::VectorType& rSolution,
                          typename TSparse::VectorType& rRhs,
                          ModelPart::DofsArrayType& rDofSet)
    {
    }

    /// @brief Compute and assemble constraint contributions into the preallocated relation matrix and constraint gap vector.
    /// @details This function is responsible for computing the entries of the relation matrix
    ///          as well as the constraint gap vector.
    /// @param rConstraints Constraint set of the related @ref ModelPart.
    /// @param rProcessInfo @ref ProcessInfo of the related @ref ModelPart.
    /// @param rDofSet Unconstrained set of @ref Dof "DoFs".
    /// @note This function must be preceded by a call to @ref ConstraintAssembler::Allocate, and should not make large scale
    ///       reallocations.
    virtual void Assemble(const ModelPart::MasterSlaveConstraintContainerType& rConstraints,
                          const ProcessInfo& rProcessInfo,
                          const ModelPart::DofsArrayType& rDofSet)
    {
    }

    /// @brief Prepare the linear system for the solution loop.
    /// @details This function is supposed to perform tasks on the linear system
    ///          that are required only once, before calls to the linear solver begin.
    ///          Constraint imposition methods that do not require a solution loop
    ///          (for example, master-slave elimination one-shots the constraints)
    ///          should manipulate the system here. If the set of @ref Dof "DoFs" has
    ///          to be changed, it should also be carried out here.
    /// @param rConstraints Constraint set of the related @ref ModelPart.
    /// @param rProcessInfo @ref ProcessInfo of the related @ref ModelPart.
    /// @param rLhs Unconstrained left hand side matrix with topology to accomodate constraint imposition.
    /// @param rRhs Unconstrained right hand side vector with space to accomodate constraint imposition.
    /// @param rDofSet Unconstrained set of @ref Dof "DoFs" with space to accomodate constraint imposition.
    virtual void Initialize(const ModelPart::MasterSlaveConstraintContainerType& rConstraints,
                            const ProcessInfo& rProcessInfo,
                            typename TSparse::MatrixType& rLhs,
                            typename TSparse::VectorType& rRhs,
                            ModelPart::DofsArrayType& rDofSet)
    {
    }

    /// @brief Manipulate the linear system before invoking the linear solver in the solution loop's current iteration.
    /// @param rConstraints Constraints of the related @ref ModelPart.
    /// @param rProcessInfo @ref ProcessInfo of the related @ref ModelPart.
    /// @param rLhs Left hand side matrix.
    /// @param rSolution Unconverged solution vector.
    /// @param rRhs Right hand side vector.
    /// @param rDofSet @ref Dof "DoFs" to solve for.
    /// @param iIteration 1-based index of the current iteration in the solution loop.
    virtual void InitializeSolutionStep(const ModelPart::MasterSlaveConstraintContainerType& rConstraints,
                                        const ProcessInfo& rProcessInfo,
                                        typename TSparse::MatrixType& rLhs,
                                        typename TSparse::VectorType& rSolution,
                                        typename TSparse::VectorType& rRhs,
                                        const ModelPart::DofsArrayType& rDofSet,
                                        const std::size_t iIteration)
    {
    }

    /// @brief Return type of @ref ConstraintAssembler::FinalizeSolutionStep.
    /// @details This class indicates
    ///          - whether the solution loop converged, and
    ///          - whether the constraint imposition is finished.
    ///          The members of this class can be in 3 valid configurations:
    ///          - constraints have not converged (@p converged is @p false) and more iterations are requested (@p finished is @p false).
    ///          - constraints have not converged (@p converged is @p false) and no more iterations are requested (@p finished is @p true).
    ///          - constraints have converged (@p converged is @p true) and no more iterations are requested (@p finished is @p true).
    struct Status {
        bool finished;      ///< @brief Indicates whether constraint imposition is finished (no more iterations are requested).
        bool converged;     ///< @brief Indicates whether constraint imposition converged.
    }; // struct Status

    /// @brief Perform constraint-related tasks after invoking the linear solver in the current iteration of the solution loop.
    /// @details This function is supposed to evaluate the convergence of constraint imposition,
    ///          decide whether to request more iterations in the solution loop.
    /// @param rConstraints Constraints of the related @ref ModelPart.
    /// @param rProcessInfo @ref ProcessInfo of the related @ref ModelPart.
    /// @param rLhs Constrained left hand side matrix.
    /// @param rSolution Converged solution vector.
    /// @param rRhs Constrained right hand side vector.
    /// @param rDofSet Constrained set of @ref Dof "DoFs".
    /// @param iIteration 1-based index of the current iteration in the solution loop.
    /// @warning The solution loop will continue indefinitely unless this function eventually
    ///          returns a @ref ConstraintImposition::Status whose @ref ConstraintImposition::Status::finished "finished"
    ///          is @p true.
    virtual Status FinalizeSolutionStep(const ModelPart::MasterSlaveConstraintContainerType& rConstraints,
                                        const ProcessInfo& rProcessInfo,
                                        typename TSparse::MatrixType& rLhs,
                                        typename TSparse::VectorType& rSolution,
                                        typename TSparse::VectorType& rRhs,
                                        const ModelPart::DofsArrayType& rDofSet,
                                        const std::size_t iIteration) = 0;

    /// @brief Perform tasks related to constraint imposition after constraints converged.
    /// @param rConstraints Constraints of the related @ref ModelPart.
    /// @param rProcessInfo @ref ProcessInfo of the related @ref ModelPart.
    /// @param rLhs Constrained left hand side matrix.
    /// @param rSolution Converged solution vector.
    /// @param rRhs Constrained right hand side vector.
    /// @param rDofSet Constrained set of @ref Dof "DoFs".
    virtual void Finalize(const ModelPart::MasterSlaveConstraintContainerType& rConstraints,
                          const ProcessInfo& rProcessInfo,
                          typename TSparse::MatrixType& rLhs,
                          typename TSparse::VectorType& rSolution,
                          typename TSparse::VectorType& rRhs,
                          ModelPart::DofsArrayType& rDofSet)
    {
    }

    /// @brief Release memory related to the linear system and constraints.
    /// @details Derived classes must call the @ref ConstraintAssembler::Clear "Clear"
    ///          function of their parents at some point.
    virtual void Clear()
    {
        mRelationMatrix = typename TSparse::MatrixType();
        mConstraintGapVector = typename TSparse::VectorType();
    }

    ConstraintImposition GetImposition() const
    {
        const std::string method_name = this->GetValue(ConstraintAssembler::GetImpositionVariable());
        if (method_name == "master_slave") {
            return ConstraintImposition::MasterSlave;
        } else if (method_name == "augmented_lagrange") {
            return ConstraintImposition::AugmentedLagrange;
        } else {
            KRATOS_ERROR << "Unsupported constraint imposition: \"" << method_name << "\"";
        }
    }

    const typename TSparse::MatrixType& GetRelationMatrix() const noexcept
    {
        return mRelationMatrix;
    }

    const typename TSparse::VectorType& GetConstraintGapVector() const noexcept
    {
        return mConstraintGapVector;
    }

protected:
    typename TSparse::MatrixType& GetRelationMatrix() noexcept
    {
        return mRelationMatrix;
    }

    typename TSparse::VectorType& GetConstraintGapVector() noexcept
    {
        return mConstraintGapVector;
    }

private:
    static const Variable<std::string>& GetImpositionVariable() noexcept
    {
        return IDENTIFIER;
    }

    typename TSparse::MatrixType mRelationMatrix;

    typename TSparse::VectorType mConstraintGapVector;
}; // class ConstraintImposition


template <class TSparse, class TDense>
class MasterSlaveConstraintAssembler : public ConstraintAssembler<TSparse,TDense>
{
public:
    using Base = ConstraintAssembler<TSparse,TDense>;

    MasterSlaveConstraintAssembler() noexcept
        : MasterSlaveConstraintAssembler(Parameters())
    {}

    MasterSlaveConstraintAssembler(Parameters Settings)
        : Base(ConstraintImposition::MasterSlave),
          mSlaveIds(),
          mMasterIds(),
          mInactiveSlaveIds(),
          mDiagonalScaling(DiagonalScaling::None)
    {
        KRATOS_TRY
        Settings.ValidateAndAssignDefaults(MasterSlaveConstraintAssembler::GetDefaultParameters());
        mDiagonalScaling = ParseDiagonalScaling(Settings);
        KRATOS_CATCH("")
    }

    /// @copydoc Base::Allocate
    void Allocate(const ModelPart::MasterSlaveConstraintContainerType& rConstraints,
                  const ProcessInfo& rProcessInfo,
                  typename TSparse::MatrixType& rLhs,
                  typename TSparse::VectorType& rSolution,
                  typename TSparse::VectorType& rRhs,
                  ModelPart::DofsArrayType& rDofSet) override
    {
        KRATOS_TRY
        std::vector<std::unordered_set<IndexType>> indices(rDofSet.size());
        std::vector<LockObject> mutexes(indices.size());
        const auto it_const_begin = rConstraints.begin();

        #pragma omp parallel
        {
            Element::EquationIdVectorType slave_ids;
            Element::EquationIdVectorType master_ids;
            std::unordered_map<IndexType, std::unordered_set<IndexType>> temp_indices;

            #pragma omp for nowait
            for (int i_const = 0; i_const < static_cast<int>(rConstraints.size()); ++i_const) {
                auto it_const = it_const_begin + i_const;
                it_const->EquationIdVector(slave_ids, master_ids, rProcessInfo);

                // Slave DoFs
                for (auto &id_i : slave_ids) {
                    temp_indices[id_i].insert(master_ids.begin(), master_ids.end());
                }
            }

            // Merging temporary indices into global rows.
            for (auto& pair_temp_indices : temp_indices) {
                std::scoped_lock<LockObject> lock(mutexes[pair_temp_indices.first]);
                indices[pair_temp_indices.first].insert(pair_temp_indices.second.begin(), pair_temp_indices.second.end());
            }
        }

        mSlaveIds.clear();
        mMasterIds.clear();
        for (int i = 0; i < static_cast<int>(indices.size()); ++i) {
            if (indices[i].size() == 0) // Master dof!
                mMasterIds.push_back(i);
            else // Slave dof
                mSlaveIds.push_back(i);
            indices[i].insert(i); // Ensure that the diagonal is there in T
        }

        MakeSparseTopology<false,typename TSparse::DataType>(indices,
                                                             indices.size(),
                                                             this->GetRelationMatrix());
        this->GetConstraintGapVector().resize(indices.size(), false);

        KRATOS_CATCH("")

        KRATOS_TRY
        rLhs = SparseMatrixMultiplicationUtility::MergeMatrices<typename TSparse::DataType>(rLhs, this->GetRelationMatrix());
        KRATOS_CATCH("")
    }

    /// @copydoc Base::Assemble
    void Assemble(const ModelPart::MasterSlaveConstraintContainerType& rConstraints,
                  const ProcessInfo& rProcessInfo,
                  const ModelPart::DofsArrayType& rDofSet) override
    {
        KRATOS_TRY

        // Init.
        mInactiveSlaveIds.clear();
        TSparse::SetToZero(this->GetRelationMatrix());
        TSparse::SetToZero(this->GetConstraintGapVector());
        std::vector<LockObject> mutexes(rDofSet.size());

        // Declare local containers.
        typename TDense::MatrixType local_relation_matrix;
        typename TDense::VectorType local_constraint_gap_vector;
        Element::EquationIdVectorType slave_equation_ids, master_equation_ids;
        const int number_of_constraints = static_cast<int>(rConstraints.size());

        #pragma omp parallel firstprivate(local_relation_matrix, local_constraint_gap_vector, slave_equation_ids, master_equation_ids)
        {
            std::unordered_set<IndexType> tls_inactive_slave_dofs;

            #pragma omp for
            for (int i_const = 0; i_const < number_of_constraints; ++i_const) {
                auto it_const = rConstraints.begin() + i_const;

                // If the constraint is active
                if (it_const->IsActive()) {
                    it_const->EquationIdVector(slave_equation_ids,
                                               master_equation_ids,
                                               rProcessInfo);
                    it_const->CalculateLocalSystem(local_relation_matrix,
                                                   local_constraint_gap_vector,
                                                   rProcessInfo);

                    for (IndexType i = 0; i < slave_equation_ids.size(); ++i) {
                        const IndexType i_global = slave_equation_ids[i];

                        // Assemble matrix row.
                        {
                            std::scoped_lock<LockObject> lock(mutexes[i_global]);
                            MapRowContribution<TSparse,TDense>(this->GetRelationMatrix(),
                                                               local_relation_matrix,
                                                               i_global,
                                                               i,
                                                               master_equation_ids);
                        }

                        // Assemble constant vector
                        KRATOS_ERROR_IF_NOT(i_global < this->GetConstraintGapVector().size())
                            << "constraint gap vector size " << this->GetConstraintGapVector().size()
                            << " <= insertion index " << i_global << "\n";
                        const double constant_value = local_constraint_gap_vector[i];
                        AtomicAdd(this->GetConstraintGapVector()[i_global], constant_value);
                    }
                } else { // Taking into account inactive constraints
                    tls_inactive_slave_dofs.insert(slave_equation_ids.begin(), slave_equation_ids.end());
                }
            }

            // We merge all the sets in one thread
            #pragma omp critical
            {
                mInactiveSlaveIds.insert(tls_inactive_slave_dofs.begin(), tls_inactive_slave_dofs.end());
            }
        }

        // Setting the master dofs into the T and C system
        for (auto eq_id : mMasterIds) {
            std::scoped_lock<LockObject> lock(mutexes[eq_id]);
            this->GetConstraintGapVector()[eq_id] = 0.0;
            this->GetRelationMatrix()(eq_id, eq_id) = 1.0;
        }

        // Setting inactive slave dofs in the T and C system
        for (auto eq_id : mInactiveSlaveIds) {
            std::scoped_lock<LockObject> lock(mutexes[eq_id]);
            this->GetConstraintGapVector()[eq_id] = 0.0;
            this->GetRelationMatrix()(eq_id, eq_id) = 1.0;
        }

        KRATOS_CATCH("")
    }

    /// @copydoc Base::Initialize
    void Initialize(const ModelPart::MasterSlaveConstraintContainerType& rConstraints,
                    const ProcessInfo& rProcessInfo,
                    typename TSparse::MatrixType& rLhs,
                    typename TSparse::VectorType& rRhs,
                    ModelPart::DofsArrayType& rDofSet) override
    {
        KRATOS_TRY
        // Compute the transposed matrix of the global relation matrix
        {
            // Storage for an intermediate matrix transpose(relation_matrix) * lhs_matrix
            typename TSparse::MatrixType left_multiplied_lhs(this->GetRelationMatrix().size2(), rLhs.size2());
            {
                // Transpose the relation matrix
                typename TSparse::MatrixType transposed_relation_matrix(this->GetRelationMatrix().size2(), this->GetRelationMatrix().size1());
                SparseMatrixMultiplicationUtility::TransposeMatrix(transposed_relation_matrix, this->GetRelationMatrix(), 1.0);

                typename TSparse::VectorType b_modified(rRhs.size());
                TSparse::Mult(transposed_relation_matrix, rRhs, b_modified);
                TSparse::Copy(b_modified, rRhs);

                SparseMatrixMultiplicationUtility::MatrixMultiplication(transposed_relation_matrix, rLhs, left_multiplied_lhs);
            } // deallocate transposed_relation_matrix

            SparseMatrixMultiplicationUtility::MatrixMultiplication(left_multiplied_lhs, this->GetRelationMatrix(), rLhs);
        } // deallocate left_multiplied_lhs

        // Compute the scale factor for slave DoFs.
        typename TSparse::DataType diagonal_scale_factor = 0;
        GetDiagonalScaleFactor<TSparse>(diagonal_scale_factor,
                                        rLhs,
                                        this->mDiagonalScaling,
                                        rProcessInfo);

        // Apply diagonal values on slaves.
        block_for_each(this->mSlaveIds, [this, &rLhs, &rRhs, diagonal_scale_factor](const auto iSlave){
            if (this->mInactiveSlaveIds.find(iSlave) == this->mInactiveSlaveIds.end()) {
                rLhs(iSlave, iSlave) = diagonal_scale_factor;
                rRhs[iSlave] = 0.0;
            }
        });
        KRATOS_CATCH("")
    }

    /// @copydoc Base::FinalizeSolutionStep
    typename Base::Status FinalizeSolutionStep(const ModelPart::MasterSlaveConstraintContainerType& rConstraints,
                                               const ProcessInfo& rProcessInfo,
                                               typename TSparse::MatrixType& rLhs,
                                               typename TSparse::VectorType& rSolution,
                                               typename TSparse::VectorType& rRhs,
                                               const ModelPart::DofsArrayType& rDofSet,
                                               const std::size_t iIteration) override
    {
        return typename Base::Status {/*finished=*/true, /*converged=*/true};
    }

    /// @copydoc Base::Finalize
    void Finalize(const ModelPart::MasterSlaveConstraintContainerType& rConstraints,
                  const ProcessInfo& rProcessInfo,
                  typename TSparse::MatrixType& rLhs,
                  typename TSparse::VectorType& rSolution,
                  typename TSparse::VectorType& rRhs,
                  ModelPart::DofsArrayType& rDofSet) override
    {
        KRATOS_TRY
        typename TSparse::VectorType original_solution(rSolution.size());
        TSparse::SetToZero(original_solution);
        TSparse::Mult(this->GetRelationMatrix(),
                      rSolution,
                      original_solution);
        rSolution.swap(original_solution);
        KRATOS_CATCH("")
    }

    /// @copydoc Base::Clear
    void Clear() override
    {
        Base::Clear();
        mSlaveIds = decltype(mSlaveIds)();
        mMasterIds = decltype(mMasterIds)();
        mInactiveSlaveIds = decltype(mInactiveSlaveIds)();
    }

    static Parameters GetDefaultParameters()
    {
        return Parameters(R"({
    "method" : "master_slave",
    "diagonal_scaling" : "none"
})");
    }

private:
    std::vector<std::size_t> mSlaveIds;

    std::vector<std::size_t> mMasterIds;

    std::unordered_set<std::size_t> mInactiveSlaveIds;

    DiagonalScaling mDiagonalScaling;
}; // class MasterSlaveConstraintAssembler


template <class TSparse, class TDense>
class AugmentedLagrangeConstraintAssembler : public ConstraintAssembler<TSparse,TDense>
{
public:
    using Base = ConstraintAssembler<TSparse,TDense>;

    AugmentedLagrangeConstraintAssembler() noexcept
        : AugmentedLagrangeConstraintAssembler(Parameters())
    {}

    AugmentedLagrangeConstraintAssembler(Parameters Settings)
        : Base(ConstraintImposition::AugmentedLagrange),
          mSlaveToConstraintMap(),
          mTransposeRelationMatrix(),
          mVerbosity(1)
    {
        Settings.ValidateAndAssignDefaults(AugmentedLagrangeConstraintAssembler::GetDefaultParameters());

        Vector algorithmic_parameters(3);
        algorithmic_parameters[0] = Settings["penalty_factor"].Get<double>();
        algorithmic_parameters[1] = Settings["initial_lagrange_multiplier"].Get<double>();
        algorithmic_parameters[2] = Settings["tolerance"].Get<double>();
        this->SetValue(this->GetAlgorithmicParametersVariable(), algorithmic_parameters);
        this->SetValue(NL_ITERATION_NUMBER, Settings["max_iterations"].Get<int>());
        this->mVerbosity = Settings["verbosity"].Get<int>();
    }

    /// @copydoc Base::Allocate
    void Allocate(const ModelPart::MasterSlaveConstraintContainerType& rConstraints,
                  const ProcessInfo& rProcessInfo,
                  typename TSparse::MatrixType& rLhs,
                  typename TSparse::VectorType& rSolution,
                  typename TSparse::VectorType& rRhs,
                  ModelPart::DofsArrayType& rDofSet) override
    {
        KRATOS_TRY

        // Construct a map that associates slave indices with constraint equation indices.
        mSlaveToConstraintMap.clear();

        {
            MasterSlaveConstraint::IndexType i_constraint = 0;
            MasterSlaveConstraint::EquationIdVectorType slave_ids, master_ids;
            for (const auto& r_constraint : rConstraints) {
                r_constraint.EquationIdVector(slave_ids, master_ids, rProcessInfo);
                for (const auto i_slave : slave_ids) {
                    const auto emplace_result = mSlaveToConstraintMap.emplace(i_slave, i_constraint);
                    if (emplace_result.second) ++i_constraint;
                } // for i_slave in slave_ids
            } // for r_constraint in rModelPart.MasterSlaveConstraints
        }

        {
            struct TLS {
                Element::EquationIdVectorType master_ids, slave_ids;
            }; // struct TLS

            std::vector<std::unordered_set<IndexType>> indices(mSlaveToConstraintMap.size());
            std::vector<LockObject> mutexes(indices.size());

            block_for_each(rConstraints,
                           TLS(),
                           [&mutexes, &indices, &rProcessInfo, this](const auto& r_constraint, TLS& r_tls){
                r_constraint.EquationIdVector(r_tls.slave_ids,
                                              r_tls.master_ids,
                                              rProcessInfo);

                for (const auto i_slave : r_tls.slave_ids) {
                    const auto i_constraint = mSlaveToConstraintMap[i_slave];
                    std::scoped_lock<LockObject> lock(mutexes[i_constraint]);
                    indices[i_constraint].insert(r_tls.master_ids.begin(), r_tls.master_ids.end());
                    indices[i_constraint].insert(i_slave);
                } // for i_slave in slave_ids
            }); // for r_constraint in rModelPart.MasterSlaveConstraints()

            MakeSparseTopology<false,typename TSparse::DataType>(indices,
                                                                 rDofSet.size(),
                                                                 this->GetRelationMatrix());

            this->GetConstraintGapVector().resize(mSlaveToConstraintMap.size(), false);
        }

        {
                typename TSparse::MatrixType product;
            {
                typename TSparse::MatrixType transpose;
                SparseMatrixMultiplicationUtility::TransposeMatrix(transpose, this->GetRelationMatrix());
                SparseMatrixMultiplicationUtility::MatrixMultiplication(transpose, this->GetRelationMatrix(), product);
            }
            rLhs = SparseMatrixMultiplicationUtility::MergeMatrices<typename TSparse::DataType>(rLhs, product);
        }

        KRATOS_CATCH("")
    }

    /// @copydoc Base::Assemble
    void Assemble(const ModelPart::MasterSlaveConstraintContainerType& rConstraints,
                  const ProcessInfo& rProcessInfo,
                  const ModelPart::DofsArrayType& rDofSet) override
    {
        KRATOS_TRY

        // Function-wide variables.
        std::vector<LockObject> mutexes(mSlaveToConstraintMap.size());

        // Init.
        TSparse::SetToZero(this->GetRelationMatrix());
        TSparse::SetToZero(this->GetConstraintGapVector());

        // Thread-local storage for constraint assembly.
        struct TLS {
            MasterSlaveConstraint::EquationIdVectorType slave_ids, master_ids;
            typename TDense::MatrixType local_relation_matrix;
            typename TDense::VectorType local_constraint_gap_vector;
        }; // struct TLS

        // Constraint assembly.
        block_for_each(rConstraints,
                       TLS(),
                       [&mutexes, &rProcessInfo, this](const MasterSlaveConstraint& r_constraint, TLS& r_tls){
            if (r_constraint.IsActive()) {
                r_constraint.EquationIdVector(r_tls.slave_ids,
                                              r_tls.master_ids,
                                              rProcessInfo);
                r_constraint.CalculateLocalSystem(r_tls.local_relation_matrix,
                                                  r_tls.local_constraint_gap_vector,
                                                  rProcessInfo);

                // Slaves also need to be part of the relation matrix in this formulation,
                // so the local relation matrix one must be appended.
                if (r_tls.local_relation_matrix.size2()) {
                    r_tls.local_relation_matrix.resize(r_tls.local_relation_matrix.size1(),
                                                       r_tls.local_relation_matrix.size2() + 1,
                                                       true);
                    for (std::size_t i_row=0; i_row<r_tls.local_relation_matrix.size1(); ++i_row) {
                        r_tls.local_relation_matrix(i_row, r_tls.local_relation_matrix.size2() - 1) = -1;
                    } // for i_row in range(local_relation_matrix.size2)
                } // if local_relation_matrix.size2

                // Assemble local rows into the global relation matrix.
                for (std::size_t i_row=0ul; i_row<r_tls.slave_ids.size(); ++i_row) {
                    const auto it_constraint_index = mSlaveToConstraintMap.find(r_tls.slave_ids[i_row]);
                    KRATOS_ERROR_IF(it_constraint_index == mSlaveToConstraintMap.end())
                        << "slave " << i_row << " (Dof " << r_tls.slave_ids[i_row] << ") "
                        << "in constraint " << r_constraint.Id() << " "
                        << "is not associated with a constraint index.";

                    const std::size_t i_constraint = it_constraint_index->second;
                    KRATOS_ERROR_IF(this->GetRelationMatrix().size1() <= i_constraint);

                    auto dof_ids = r_tls.master_ids;
                    dof_ids.push_back(r_tls.slave_ids[i_row]);
                    KRATOS_ERROR_IF_NOT(dof_ids.size() == r_tls.local_relation_matrix.size2());

                    // Indirect sort the local relation matrix' row based on DoF IDs.
                    // This step is required because the global relation matrix is in CSR format,
                    // and Impl::MapRowContribution expects the column indices to be sorted.
                    {
                        std::vector<std::size_t> index_array(dof_ids.size());
                        std::iota(index_array.begin(), index_array.end(), 0ul);
                        std::sort(index_array.begin(),
                                  index_array.end(),
                                  [&dof_ids](const std::size_t i_left, const std::size_t i_right) -> bool {
                                      return dof_ids[i_left] < dof_ids[i_right];
                                  });

                        // Reorder the current row in the local relation matrix.
                        std::vector<MasterSlaveConstraint::MatrixType::value_type> row(r_tls.local_relation_matrix.size2());
                        std::transform(index_array.begin(),
                                       index_array.end(),
                                       row.begin(),
                                       [&r_tls, i_row](const std::size_t i_column){
                                        return r_tls.local_relation_matrix(i_row, i_column);
                                       });
                        std::copy(row.begin(),
                                  row.end(),
                                  r_tls.local_relation_matrix.data().begin() + i_row * r_tls.local_relation_matrix.size2());

                        // Reorder DoF indices.
                        auto swap_dof_ids = dof_ids;
                        std::transform(index_array.begin(),
                                       index_array.end(),
                                       swap_dof_ids.begin(),
                                       [&dof_ids](const auto i_column){
                                        return dof_ids[i_column];
                                       });
                        dof_ids.swap(swap_dof_ids);
                    }

                    {
                        std::scoped_lock<LockObject> lock(mutexes[i_constraint]);
                        MapRowContribution<TSparse,TDense>(this->GetRelationMatrix(),
                                                           r_tls.local_relation_matrix,
                                                           i_constraint,
                                                           i_row,
                                                           dof_ids);
                    }
                } // for i_row in range(r_tls.slave_ids.size)
            } // if r_constraint.IsActive
        }); // for r_constraint in rModelPart.MasterSlaveCosntraints
        KRATOS_CATCH("")

        KRATOS_TRY
        SparseMatrixMultiplicationUtility::TransposeMatrix(mTransposeRelationMatrix, this->GetRelationMatrix());
        KRATOS_CATCH("")
    }

    /// @copydoc Base::Initialize
    void Initialize(const ModelPart::MasterSlaveConstraintContainerType& rConstraints,
                    const ProcessInfo& rProcessInfo,
                    typename TSparse::MatrixType& rLhs,
                    typename TSparse::VectorType& rRhs,
                    ModelPart::DofsArrayType& rDofSet) override
    {
        KRATOS_TRY

        using SparseUtils = SparseMatrixMultiplicationUtility;
        const typename TSparse::DataType penalty_factor = this->GetPenaltyFactor();
        const typename TSparse::DataType initial_lagrange_multiplier = this->GetInitialLagrangeMultiplier();

        {
            KRATOS_ERROR_IF(mTransposeRelationMatrix.size1() != this->GetRelationMatrix().size2()
                            || mTransposeRelationMatrix.size2() != this->GetRelationMatrix().size1()
                            || mTransposeRelationMatrix.nnz() != this->GetRelationMatrix().nnz())
                << "the transpose of the relation matrix is uninitialized or out of date";

            typename TSparse::MatrixType relation_product;
            SparseUtils::MatrixMultiplication(mTransposeRelationMatrix,
                                              this->GetRelationMatrix(),
                                              relation_product);

            // Add terms to the LHS matrix.
            SparseUtils::InPlaceMatrixAdd(rLhs,
                                          relation_product,
                                          penalty_factor);

            // Add terms to the RHS vector.
            typename TSparse::VectorType rhs_term(rRhs.size()),
                                         lagrange_multipliers(this->GetRelationMatrix().size1(), initial_lagrange_multiplier),
                                         constraint_gaps = this->GetConstraintGapVector();
            KRATOS_ERROR_IF_NOT(constraint_gaps.size() == lagrange_multipliers.size());

            TSparse::UnaliasedAdd(lagrange_multipliers,
                                  -penalty_factor,
                                  constraint_gaps);
            TSparse::Mult(mTransposeRelationMatrix,
                          lagrange_multipliers,
                          rhs_term);
            TSparse::UnaliasedAdd(rRhs,
                                  -1.0,
                                  rhs_term);
        }

        KRATOS_CATCH("")
    }

    /// @copydoc Base::FinalizeSolutionStep
    typename Base::Status FinalizeSolutionStep(const ModelPart::MasterSlaveConstraintContainerType& rConstraints,
                                               const ProcessInfo& rProcessInfo,
                                               typename TSparse::MatrixType& rLhs,
                                               typename TSparse::VectorType& rSolution,
                                               typename TSparse::VectorType& rRhs,
                                               const ModelPart::DofsArrayType& rDofSet,
                                               const std::size_t iIteration) override
    {
        KRATOS_TRY
        const int max_iterations = this->GetValue(NL_ITERATION_NUMBER);
        // Compute the constraint residuals.
        typename TSparse::VectorType constraint_residual(this->GetConstraintGapVector().size());
        TSparse::SetToZero(constraint_residual);

        KRATOS_TRY
        TSparse::Mult(this->GetRelationMatrix(), rSolution, constraint_residual);
        TSparse::UnaliasedAdd(constraint_residual, -1.0, this->GetConstraintGapVector());
        KRATOS_CATCH("")

        // Decide whether to update the RHS or terminate the loop.
        const auto constraint_norm = TSparse::TwoNorm(constraint_residual);
        if (3 <= this->mVerbosity) {
            std::stringstream residual_stream;
            residual_stream << std::setprecision(8) << std::scientific<< constraint_norm;
            std::cout << "PMultigridBuilderAndSolver: iteration " << iIteration
                      << " constraint residual " << residual_stream.str() << "\n";
        }

        if (this->GetTolerance() < constraint_norm && static_cast<int>(iIteration) < max_iterations) {
            KRATOS_TRY
            typename TSparse::VectorType rhs_update(rRhs.size());
            TSparse::SetToZero(rhs_update);
            TSparse::InplaceMult(constraint_residual, -this->GetPenaltyFactor());
            TSparse::Mult(mTransposeRelationMatrix, constraint_residual, rhs_update);
            TSparse::UnaliasedAdd(rRhs, 1.0, rhs_update);
            return typename Base::Status {/*finished=*/false, /*converged=*/false};
            KRATOS_CATCH("")
        } else {
            return typename Base::Status {/*finished=*/true, /*converged=*/constraint_norm <= this->GetTolerance()};
        }
        KRATOS_CATCH("")
    }

    /// @copydoc Base::Clear
    void Clear() override
    {
        Base::Clear();
        mSlaveToConstraintMap = decltype(mSlaveToConstraintMap)();
        mTransposeRelationMatrix = decltype(mTransposeRelationMatrix)();
    }

    static Parameters GetDefaultParameters()
    {
        return Parameters(R"({
    "method" : "augmented_lagrange",
    "penalty_factor" : 1e3,
    "initial_lagrange_multiplier" : 0.0,
    "tolerance" : 1e-8,
    "max_iterations" : 1e1,
    "verbosity" : 1
})");
    }

    typename TSparse::DataType GetPenaltyFactor() const
    {
        return this->GetValue(AugmentedLagrangeConstraintAssembler::GetAlgorithmicParametersVariable())[0];
    }

    typename TSparse::DataType GetInitialLagrangeMultiplier() const
    {
        return this->GetValue(AugmentedLagrangeConstraintAssembler::GetAlgorithmicParametersVariable())[1];
    }

    typename TSparse::DataType GetTolerance() const
    {
        return this->GetValue(AugmentedLagrangeConstraintAssembler::GetAlgorithmicParametersVariable())[2];
    }

    static const Variable<Vector>& GetAlgorithmicParametersVariable() noexcept
    {
        return SHAPE_FUNCTIONS_VECTOR;
    }

private:
    /// @brief A map associating slave IDs with constraint indices.
    std::unordered_map<std::size_t,std::size_t> mSlaveToConstraintMap;

    typename TSparse::MatrixType mTransposeRelationMatrix;

    int mVerbosity;
}; // class AugmentedLagrangeConstraintAssembler


// --------------------------------------------------------- //
// PMG Hierarchy
// --------------------------------------------------------- //


template <class TSparse,
          class TDense,
          class TSolver>
class PGrid
{
public:
    PGrid(const unsigned Level)
        : mpRestrictionOperator(new typename TSparse::MatrixType),
          mpLhs(new typename TSparse::MatrixType),
          mpRhs(new typename TSparse::VectorType)
    {}


    PGrid()
        : PGrid(0u)
    {}


    template <class TParentSparse>
    void MakeLhsTopology(const ModelPart& rModelPart,
                         const typename TParentSparse::MatrixType& rParentLhs)
    {
        KRATOS_TRY
        // The restriction operator immediately constructs the linear equivalent,
        // because one-level coarsening strategies are not supported yet. The problem
        // is that when going from a generic polynomial level q to some other lower polynomial
        // level p!=1, new nodes must be introduced that do not exist in the original fine mesh.
        // This wouldn't be a huge problem by itself, but deciding on DoF indices on the coarse
        // level would be quite painful and require to keep track of the coarse grid's topology
        // in some manner.
        MakePRestrictionOperator<
            std::numeric_limits<unsigned>::max(),
            typename TSparse::DataType>(rModelPart,
                                        *mpRestrictionOperator,
                                        rParentLhs.size1());
        KRATOS_CATCH("")
    }


    template <bool AssembleLHS,
              bool AssembleRHS,
              class TParentSparse>
    void Assemble(const ModelPart& rModelPart,
                  const typename TParentSparse::MatrixType* pParentLhs,
                  const typename TParentSparse::VectorType* pParentRhs)
    {
        KRATOS_TRY

        // Assemble the coarse LHS matrix if requested.
        if constexpr (AssembleLHS) {
            KRATOS_ERROR_IF_NOT(pParentLhs);

            const auto check_matrix_multiplication = [](const typename TSparse::MatrixType& r_left,
                                                        const typename TSparse::MatrixType& r_right,
                                                        const typename TSparse::MatrixType& r_out) {
                KRATOS_ERROR_IF((r_left.size2() != r_right.size1()) ||
                                (r_out.size1() != r_left.size1()) ||
                                (r_out.size2() != r_right.size2()))
                    << "invalid sizes for matrix product: "
                    << "(" << r_left.size1() << "x" << r_left.size2() << ")"
                    << " x "
                    << "(" << r_right.size1() << "x" << r_right.size2() << ")"
                    << " => "
                    << "(" << r_out.size1() << "x" << r_out.size2() << ")";
            };

            typename TSparse::MatrixType prolongation_operator(mpRestrictionOperator->size2(), mpRestrictionOperator->size1()),
                                         left_multiplied_lhs(mpRestrictionOperator->size1(), pParentLhs->size2());

            check_matrix_multiplication(*mpRestrictionOperator, *pParentLhs, left_multiplied_lhs);
            SparseMatrixMultiplicationUtility::MatrixMultiplication(*mpRestrictionOperator, *pParentLhs, left_multiplied_lhs);
            SparseMatrixMultiplicationUtility::TransposeMatrix(prolongation_operator, *mpRestrictionOperator, 1.0);

            mpLhs->resize(mpRestrictionOperator->size1(), mpRestrictionOperator->size1(), false);
            check_matrix_multiplication(left_multiplied_lhs, prolongation_operator, *mpLhs);
            SparseMatrixMultiplicationUtility::MatrixMultiplication(left_multiplied_lhs, prolongation_operator, *mpLhs);
        } // if AssembleLHS

        // Assemble the coarse RHS vector if requested.
        if constexpr (AssembleRHS) {
            // Sanity checks
            KRATOS_ERROR_IF_NOT(pParentRhs);
            KRATOS_ERROR_IF_NOT(mpRestrictionOperator->size2() == pParentRhs->size())
                << "expecting an RHS vector of size " << mpRestrictionOperator->size2()
                << " but got " << pParentRhs->size();

            mpRhs->resize(mpRestrictionOperator->size1(), false);
            TSparse::Mult(*mpRestrictionOperator, *pParentRhs, *mpRhs);
        } // if AssembleRHS

        KRATOS_CATCH("")
    }

private:
    std::shared_ptr<typename TSparse::MatrixType> mpRestrictionOperator;

    std::shared_ptr<typename TSparse::MatrixType> mpLhs;

    std::shared_ptr<typename TSparse::VectorType> mpRhs;

    std::optional<std::unique_ptr<PGrid>> mMaybeChild;
}; // class PGrid


// --------------------------------------------------------- //
// PIMPL
// --------------------------------------------------------- //


template <class TSparse, class TDense, class TSolver>
struct PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::Impl
{
    using Interface = PMultigridBuilderAndSolver<TSparse,TDense,TSolver>;

    // --------------------------------------------------------- //
    // Member Variables
    // --------------------------------------------------------- //

    Interface* mpInterface;

    std::shared_ptr<ConstraintAssembler<TSparse,TDense>> mpConstraintAssembler;

    typename TSparse::DataType mDiagonalScaleFactor;

    PGrid<TSparse,TDense,TSolver> mHierarchy;

    DiagonalScaling mDiagonalScaling;

    int mMaxDepth;

    int mVerbosity;

    // --------------------------------------------------------- //
    // Special Member Functions
    // --------------------------------------------------------- //

    Impl(Interface* pInterface)
        : mpInterface(pInterface),
          mpConstraintAssembler(new MasterSlaveConstraintAssembler<TSparse,TDense>),
          mDiagonalScaleFactor(1),
          mHierarchy(),
          mDiagonalScaling(DiagonalScaling::None),
          mMaxDepth(-1),
          mVerbosity(0)
    {}

    // --------------------------------------------------------- //
    // Solution
    // --------------------------------------------------------- //

    /// @brief Initialize the linear solver and solve the provided system.
    void Solve(typename Interface::TSystemMatrixType& rLhs,
               typename Interface::TSystemVectorType& rSolution,
               typename Interface::TSystemVectorType& rRhs,
               ModelPart& rModelPart,
               Interface& rInterface)
    {
        KRATOS_TRY
        if (rInterface.GetLinearSolver().AdditionalPhysicalDataIsNeeded()) {
            rInterface.GetLinearSolver().ProvideAdditionalData(rLhs,
                                                               rSolution,
                                                               rRhs,
                                                               rInterface.GetDofSet(),
                                                               rModelPart);
        }
        KRATOS_CATCH("")

        std::size_t i_iteration = 0ul;
        typename ConstraintAssembler<TSparse,TDense>::Status status;
        do {
            mpConstraintAssembler->InitializeSolutionStep(rModelPart.MasterSlaveConstraints(),
                                                          rModelPart.GetProcessInfo(),
                                                          rLhs,
                                                          rSolution,
                                                          rRhs,
                                                          rInterface.GetDofSet(),
                                                          i_iteration);
            rInterface.GetLinearSolver().Solve(rLhs, rSolution, rRhs);
            status = mpConstraintAssembler->FinalizeSolutionStep(rModelPart.MasterSlaveConstraints(),
                                                                 rModelPart.GetProcessInfo(),
                                                                 rLhs,
                                                                 rSolution,
                                                                 rRhs,
                                                                 rInterface.GetDofSet(),
                                                                 ++i_iteration);
        } while (not status.finished);

        if (1 <= mVerbosity and not status.converged) {
            std::cerr << "PMultigridBuilderAndSolver: failed to converge constraints\n";
        }

        mpConstraintAssembler->Finalize(rModelPart.MasterSlaveConstraints(),
                                        rModelPart.GetProcessInfo(),
                                        rLhs,
                                        rSolution,
                                        rRhs,
                                        mpInterface->GetDofSet());
    }


    // --------------------------------------------------------- //
    // Topology
    // --------------------------------------------------------- //


    template <class TProxy>
    static void CollectDoFs(const TProxy& rEntities,
                            const ProcessInfo& rProcessInfo,
                            typename Interface::TSchemeType& rScheme,
                            LockObject* pLockBegin,
                            [[maybe_unused]] LockObject* pLockEnd,
                            std::unordered_set<std::size_t>* pRowSetBegin,
                            [[maybe_unused]] std::unordered_set<std::size_t>* pRowSetEnd)
    {
        KRATOS_TRY
        KRATOS_PROFILE_SCOPE(KRATOS_CODE_LOCATION);
        using TLS = Element::EquationIdVectorType;
        block_for_each(rEntities,
                       TLS(),
                       [&rScheme, &rProcessInfo, pLockBegin, pRowSetBegin](const auto& rEntity, TLS& rTls){
            if (rEntity.GetEntity().IsActive()) {
                rScheme.EquationId(rEntity.GetEntity(), rTls, rProcessInfo);
                for (const auto equation_id : rTls) {
                    [[maybe_unused]] std::scoped_lock<LockObject> lock(pLockBegin[equation_id]);
                    auto& r_row_indices = pRowSetBegin[equation_id];
                    r_row_indices.insert(rTls.begin(), rTls.end());
                } // for equation_id in rTls
            } // if rEntity.IsActive
        });
        KRATOS_CATCH("")
    }


    void MakeLhsTopology(const typename Interface::TSchemeType::Pointer& rpScheme,
                         typename Interface::TSystemMatrixType& rLhs,
                         ModelPart& rModelPart)
    {
        KRATOS_PROFILE_SCOPE_MILLI(KRATOS_CODE_LOCATION);
        KRATOS_TRY

        std::vector<std::unordered_set<std::size_t> > indices(mpInterface->GetEquationSystemSize());

        {
            std::vector<LockObject> mutexes(mpInterface->GetEquationSystemSize());

            // Collect DoFs from elements.
            Impl::CollectDoFs(MakeProxy<Globals::DataLocation::Element>(rModelPart),
                              rModelPart.GetProcessInfo(),
                              *rpScheme,
                              mutexes.data(),
                              mutexes.data() + mutexes.size(),
                              indices.data(),
                              indices.data() + indices.size());

            // Collect DoFs from conditions.
            Impl::CollectDoFs(MakeProxy<Globals::DataLocation::Condition>(rModelPart),
                              rModelPart.GetProcessInfo(),
                              *rpScheme,
                              mutexes.data(),
                              mutexes.data() + mutexes.size(),
                              indices.data(),
                              indices.data() + indices.size());
        }

        // Compute and allocate LHS topology.
        MakeSparseTopology<false,typename TSparse::DataType>(indices,
                                                             indices.size(),
                                                             rLhs);

        // Construct the coarse hierarhy's topology.
        mHierarchy.template MakeLhsTopology<TSparse>(rModelPart, rLhs);
        KRATOS_CATCH("")
    }


    // --------------------------------------------------------- //
    // Assembly
    // --------------------------------------------------------- //


    /// @brief Compute and assemble local contributions from elements and conditions into the global system.
    /// @details This function body mainly handles looping over elements/conditions and its parallelization.
    ///          The actual local system calculation, as well as assembly, is deferred to @ref MapEntityContribution.
    template <bool AssembleLHS,
              bool AssembleRHS>
    void Assemble(ModelPart& rModelPart,
                  typename Interface::TSchemeType& rScheme,
                  std::optional<typename Interface::TSystemMatrixType*> pMaybeLhs,
                  std::optional<typename Interface::TSystemVectorType*> pMaybeRhs)
    {
        KRATOS_TRY

        // Sanity checks.
        if constexpr (AssembleLHS) {
            KRATOS_ERROR_IF_NOT(pMaybeLhs.has_value() && pMaybeLhs.value() != nullptr)
                << "Requested an assembly of the left hand side, but no matrix is provided to assemble into.";
        }

        if constexpr (AssembleRHS) {
            KRATOS_ERROR_IF_NOT(pMaybeRhs.has_value() && pMaybeRhs.value() != nullptr)
                << "Requested an assembly of the right hand side, but no vector is provided to assemble into.";
        }

        // Global variables.
        const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
        typename Interface::TSystemMatrixType* pLhs = pMaybeLhs.has_value() ? pMaybeLhs.value() : nullptr;
        typename Interface::TSystemVectorType* pRhs = pMaybeRhs.has_value() ? pMaybeRhs.value() : nullptr;
        auto p_locks = std::make_unique<std::vector<LockObject>>(AssembleLHS ? pLhs->size1() : 0ul);

        [[maybe_unused]] const int element_count = rModelPart.Elements().size();
        [[maybe_unused]] const int condition_count = rModelPart.Conditions().size();
        [[maybe_unused]] const auto it_element_begin = rModelPart.Elements().begin();
        [[maybe_unused]] const auto it_condition_begin = rModelPart.Conditions().begin();

        // Collect contributions from constitutive entities.
        #pragma omp parallel
        {
            // Thread-local variables.
            Element::EquationIdVectorType equation_indices;
            typename Interface::LocalSystemMatrixType lhs_contribution;
            typename Interface::LocalSystemVectorType rhs_contribution;

            // Collect contributions from elements.
            #pragma omp for schedule(guided, 512) nowait
            for (int i_entity=0; i_entity<element_count; ++i_entity) {
                MapEntityContribution<TSparse,TDense,AssembleLHS,AssembleRHS>(
                    *(it_element_begin + i_entity),
                    rScheme,
                    r_process_info,
                    equation_indices,
                    &lhs_contribution,
                    &rhs_contribution,
                    pLhs,
                    pRhs,
                    p_locks->data());
            } // pragma omp for

            // Collect contributions from conditions.
            #pragma omp for schedule(guided, 512)
            for (int i_entity=0; i_entity<condition_count; ++i_entity) {
                MapEntityContribution<TSparse,TDense,AssembleLHS,AssembleRHS>(
                    *(it_condition_begin + i_entity),
                    rScheme,
                    r_process_info,
                    equation_indices,
                    &lhs_contribution,
                    &rhs_contribution,
                    pLhs,
                    pRhs,
                    p_locks->data());
            } // pragma omp for
        } // pragma omp parallel

        // Assemble coarse hierarchy.
        mHierarchy.template Assemble<AssembleLHS, AssembleRHS,TSparse>(rModelPart,
                                                                       pLhs,
                                                                       pRhs);

        KRATOS_CATCH("")
    }

}; // class PMultigridBuilderAndSolver::Impl


// --------------------------------------------------------- //
// Constructors
// --------------------------------------------------------- //


template <class TSparse, class TDense, class TSolver>
PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::~PMultigridBuilderAndSolver() = default;


template <class TSparse, class TDense, class TSolver>
PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::PMultigridBuilderAndSolver()
    : Interface(),
      mpImpl(new Impl(this))
{
}


template <class TSparse, class TDense, class TSolver>
PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::PMultigridBuilderAndSolver(const typename TSolver::Pointer& pSolver,
                                                                               Parameters Settings)
    : Interface(pSolver),
      mpImpl(new Impl(this))
{
    Settings.ValidateAndAssignDefaults(this->GetDefaultParameters());
    this->AssignSettings(Settings);
}


template <class TSparse, class TDense, class TSolver>
typename PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::Interface::Pointer
PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::Create(typename TSolver::Pointer pSolver,
                                                           Parameters Settings) const
{
    KRATOS_TRY
    return typename Interface::Pointer(new PMultigridBuilderAndSolver(pSolver, Settings));
    KRATOS_CATCH("")
}


// --------------------------------------------------------- //
// Allocation and Initialization
// --------------------------------------------------------- //


template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::SetUpDofSet(typename Interface::TSchemeType::Pointer pScheme,
                                                                     ModelPart& rModelPart)
{
    KRATOS_TRY;

    BlockBuildDofArrayUtility::SetUpDofArray(rModelPart,
                                             this->GetDofSet(),
                                             this->GetEchoLevel(),
                                             this->GetCalculateReactionsFlag());

    KRATOS_CATCH("");
}


template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::SetUpSystem(ModelPart& rModelPart)
{
    this->mEquationSystemSize = this->GetDofSet().size();

    KRATOS_TRY
    // Set equation indices of DoFs.
    IndexPartition<std::size_t>(this->GetDofSet().size()).for_each([&, this](std::size_t Index){
        (this->GetDofSet().begin() + Index)->SetEquationId(Index);
    });
    KRATOS_CATCH("")
}


template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::ResizeAndInitializeVectors(typename Interface::TSchemeType::Pointer pScheme,
                                                                                    typename Interface::TSystemMatrixPointerType& rpLhs,
                                                                                    typename Interface::TSystemVectorPointerType& rpSolution,
                                                                                    typename Interface::TSystemVectorPointerType& rpRhs,
                                                                                    ModelPart& rModelPart)
{
    KRATOS_TRY

    // Construct empty containers if necessary.
    if (!rpLhs)
        rpLhs.reset(new typename Interface::TSystemMatrixType);

    if (!rpSolution)
        rpSolution.reset(new typename Interface::TSystemVectorType);
    if (rpSolution->size() != this->mEquationSystemSize)
        rpSolution->resize(this->mEquationSystemSize, false);
    TSparse::SetToZero(*rpSolution);

    if (!rpRhs)
        rpRhs.reset(new typename Interface::TSystemVectorType);
    if (rpRhs->size() != this->mEquationSystemSize)
        rpRhs->resize(this->mEquationSystemSize, false);
    TSparse::SetToZero(*rpRhs);

    // Construct LHS topology if necessary or requested.
    if (rpLhs->size1() == 0 || this->GetReshapeMatrixFlag() == true) {
        rpLhs->resize(this->mEquationSystemSize, this->mEquationSystemSize, false);
        mpImpl->MakeLhsTopology(pScheme, *rpLhs, rModelPart);

        // Make constraint topology.
        mpImpl->mpConstraintAssembler->Allocate(rModelPart.MasterSlaveConstraints(),
                                                rModelPart.GetProcessInfo(),
                                                *rpLhs,
                                                *rpSolution,
                                                *rpRhs,
                                                this->GetDofSet());
    } else {
        if (rpLhs->size1() != this->mEquationSystemSize || rpLhs->size2() != this->mEquationSystemSize) {
            KRATOS_ERROR <<"The equation system size has changed during the simulation. This is not permitted."<<std::endl;
        }
    }

    KRATOS_CATCH("")
}


// --------------------------------------------------------- //
// Hooks
// --------------------------------------------------------- //


template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse, TDense, TSolver>::InitializeSolutionStep(ModelPart& rModelPart,
                                                                                  typename Interface::TSystemMatrixType& rLhs,
                                                                                  typename Interface::TSystemVectorType& rSolution,
                                                                                  typename Interface::TSystemVectorType& rRhs)
{
    KRATOS_TRY
    //
    KRATOS_CATCH("")
}


// --------------------------------------------------------- //
// Assembly
// --------------------------------------------------------- //


template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::Build(typename Interface::TSchemeType::Pointer pScheme,
                                                               ModelPart& rModelPart,
                                                               typename Interface::TSystemMatrixType& rLhs,
                                                               typename Interface::TSystemVectorType& rRhs)
{
    KRATOS_PROFILE_SCOPE_MILLI(KRATOS_CODE_LOCATION);
    KRATOS_ERROR_IF(!pScheme) << "missing scheme" << std::endl;
    KRATOS_TRY
    mpImpl->template Assemble</*AssembleLHS=*/true,/*AssembleRHS=*/true>(rModelPart,
                                                                         *pScheme,
                                                                         &rLhs,
                                                                         &rRhs);
    KRATOS_CATCH("")
}


template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::BuildLHS(typename Interface::TSchemeType::Pointer pScheme,
                                                                  ModelPart& rModelPart,
                                                                  typename Interface::TSystemMatrixType& rLhs)
{
    KRATOS_PROFILE_SCOPE_MILLI(KRATOS_CODE_LOCATION);
    KRATOS_TRY
    KRATOS_ERROR_IF(!pScheme) << "missing scheme" << std::endl;
    mpImpl->template Assemble</*AssembleLHS=*/true,/*AssembleRHS=*/false>(rModelPart,
                                                                          *pScheme,
                                                                          &rLhs,
                                                                          nullptr);
    KRATOS_CATCH("")
}


template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::BuildRHS(typename Interface::TSchemeType::Pointer pScheme,
                                                                  ModelPart& rModelPart,
                                                                  typename Interface::TSystemVectorType& rRhs)
{
    KRATOS_PROFILE_SCOPE_MILLI(KRATOS_CODE_LOCATION);
    KRATOS_TRY
    mpImpl->template Assemble</*AssembleLHS=*/false,/*AssembleRHS=*/true>(rModelPart,
                                                                          *pScheme,
                                                                          nullptr,
                                                                          &rRhs);
    block_for_each(this->GetDofSet(), [&rRhs](Dof<double>& rDof){
        if (rDof.IsFixed())
            rRhs[rDof.EquationId()] = 0.0;
    });

    KRATOS_CATCH("")
}


// --------------------------------------------------------- //
// Constraint Imposition
// --------------------------------------------------------- //


template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::ApplyDirichletConditions(typename Interface::TSchemeType::Pointer pScheme,
                                                                                  ModelPart& rModelPart,
                                                                                  typename Interface::TSystemMatrixType& rLhs,
                                                                                  typename Interface::TSystemVectorType& rSolution,
                                                                                  typename Interface::TSystemVectorType& rRhs)
{
    KRATOS_PROFILE_SCOPE(KRATOS_CODE_LOCATION);
    const std::size_t system_size = rLhs.size1();
    Vector scaling_factors (system_size);

    const auto it_dof_iterator_begin = this->GetDofSet().begin();

    // NOTE: dofs are assumed to be numbered consecutively in the BlockBuilderAndSolver
    IndexPartition<std::size_t>(this->GetDofSet().size()).for_each([&](std::size_t Index){
        auto it_dof_iterator = it_dof_iterator_begin + Index;
        if (it_dof_iterator->IsFixed()) {
            scaling_factors[Index] = 0.0;
        } else {
            scaling_factors[Index] = 1.0;
        }
    });

    //GetDiagonalScaleFactor<TSparse>(mpImpl->mDiagonalScaleFactor,
    //                                rLhs,
    //                                mpImpl->mDiagonalScaling,
    //                                rModelPart.GetProcessInfo());

    double* Avalues = rLhs.value_data().begin();
    std::size_t* Arow_indices = rLhs.index1_data().begin();
    std::size_t* Acol_indices = rLhs.index2_data().begin();

    IndexPartition<std::size_t>(system_size).for_each([&](std::size_t Index){
        const std::size_t col_begin = Arow_indices[Index];
        const std::size_t col_end = Arow_indices[Index+1];
        const double k_factor = scaling_factors[Index];
        if (k_factor == 0.0) {
            // Zero out the whole row, except the diagonal
            for (std::size_t j = col_begin; j < col_end; ++j)
                if (Acol_indices[j] != Index )
                    Avalues[j] = 0.0;
            // Zero out the RHS
            rRhs[Index] = 0.0;
        } else {
            // Zero out the column which is associated with the zero'ed row
            for (std::size_t j = col_begin; j < col_end; ++j)
                if(scaling_factors[ Acol_indices[j] ] == 0 )
                    Avalues[j] = 0.0;
        }
    });
}


template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::ApplyConstraints(typename Interface::TSchemeType::Pointer pScheme,
                                                                          ModelPart& rModelPart,
                                                                          typename Interface::TSystemMatrixType& rLhs,
                                                                          typename Interface::TSystemVectorType& rRhs)
{
    KRATOS_PROFILE_SCOPE_MILLI(KRATOS_CODE_LOCATION);
    KRATOS_TRY
    mpImpl->mpConstraintAssembler->Assemble(
        rModelPart.MasterSlaveConstraints(),
        rModelPart.GetProcessInfo(),
        this->GetDofSet());
    mpImpl->mpConstraintAssembler->Initialize(
        rModelPart.MasterSlaveConstraints(),
        rModelPart.GetProcessInfo(),
        rLhs,
        rRhs,
        this->GetDofSet());
    KRATOS_CATCH("")
}


// --------------------------------------------------------- //
// Compound Assembly and Solution
// --------------------------------------------------------- //


template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::BuildAndSolve(typename Interface::TSchemeType::Pointer pScheme,
                                                                       ModelPart& rModelPart,
                                                                       typename Interface::TSystemMatrixType& rLhs,
                                                                       typename Interface::TSystemVectorType& rSolution,
                                                                       typename Interface::TSystemVectorType& rRhs)
{
    KRATOS_TRY

    // Assemble unconstrained system.
    Build(pScheme, rModelPart, rLhs, rRhs);

    // Apply multifreedom constraints.
    ApplyConstraints(pScheme, rModelPart, rLhs, rRhs);

    // Apply Dirichlet conditions.
    ApplyDirichletConditions(pScheme, rModelPart, rLhs, rSolution, rRhs);

    // Solve constrained assembled system.
    mpImpl->Solve(rLhs, rSolution, rRhs, rModelPart, *this);
    KRATOS_CATCH("")
}


// --------------------------------------------------------- //
// Postprocessing
// --------------------------------------------------------- //



template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::CalculateReactions(typename Interface::TSchemeType::Pointer pScheme,
                                                                            ModelPart& rModelPart,
                                                                            typename Interface::TSystemMatrixType& rLhs,
                                                                            typename Interface::TSystemVectorType& rSolution,
                                                                            typename Interface::TSystemVectorType& rRhs)
{
    KRATOS_PROFILE_SCOPE_MILLI(KRATOS_CODE_LOCATION);
    TSparse::SetToZero(rRhs);
    mpImpl->template Assemble</*AssembleLHS=*/false,/*AssembleRHS=*/true>(rModelPart,
                                                                          *pScheme,
                                                                          nullptr,
                                                                          &rRhs);
    block_for_each(this->GetDofSet(), [&rRhs](Dof<double>& rDof){
        rDof.GetSolutionStepReactionValue() = -rRhs[rDof.EquationId()];
    });
}


// --------------------------------------------------------- //
// Misc
// --------------------------------------------------------- //


template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::AssignSettings(const Parameters Settings)
{
    KRATOS_TRY
    Parameters(Settings).ValidateAndAssignDefaults(this->GetDefaultParameters());
    KRATOS_CATCH("")

    KRATOS_TRY
    Interface::AssignSettings(Settings);
    KRATOS_CATCH("")

    // Set the scaling strategy for the diagonal entries of constrained DoFs.
    mpImpl->mDiagonalScaling = ParseDiagonalScaling(Settings);

    // Set multifreedom constraint imposition strategy.
    KRATOS_TRY
    Parameters constraint_imposition_settings = Settings["constraint_imposition"];
    const std::string constraint_imposition_name = constraint_imposition_settings["method"].Get<std::string>();
    if (constraint_imposition_name == "master_slave_elimination") {
        mpImpl->mpConstraintAssembler = std::make_shared<MasterSlaveConstraintAssembler<TSparse,TDense>>(constraint_imposition_settings);
    } else if (constraint_imposition_name == "augmented_lagrange") {
        mpImpl->mpConstraintAssembler = std::make_shared<AugmentedLagrangeConstraintAssembler<TSparse,TDense>>(constraint_imposition_settings);
    } else {
        std::stringstream message;
        message << "Unsupported constraint imposition \"" << constraint_imposition_name << "\". Options are:\n";
        message << "\t\"master_slave_elimination\"\n";
        message << "\t\"augmented_lagrange\"";
        KRATOS_ERROR << message.str();
    }
    KRATOS_CATCH("")

    // Other settings.
    KRATOS_TRY
    mpImpl->mMaxDepth = Settings["max_depth"].Get<int>();
    mpImpl->mVerbosity = Settings["verbosity"].Get<int>();
    KRATOS_CATCH("")
}


template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::Clear()
{
    Interface::Clear();
    mpImpl->mpConstraintAssembler->Clear();
}


template <class TSparse, class TDense, class TSolver>
Parameters PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::GetDefaultParameters() const
{
    Parameters parameters = Parameters(R"({
"name"              : "p_multigrid",
"diagonal_scaling"  : "abs_max",
"max_depth"         : -1,
"smoother_settings" : {
    "solver_type" : "gauss_seidel"
},
"solver_settings"   : {
    "solver_type" : "amgcl"
},
"constraint_imposition" : {
    "method" : "master_slave_elimination"
},
"verbosity"         : 0
})");
    parameters.RecursivelyAddMissingParameters(Interface::GetDefaultParameters());
    return parameters;
}


template <class TSparse, class TDense, class TSolver>
std::string PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::Info() const
{
    return "PMultigridBuilderAndSolver";
}


template <class TSparse, class TDense, class TSolver>
std::size_t PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::GetEquationSystemSize() const noexcept
{
    return Interface::mEquationSystemSize;
}

template <class TSparse, class TDense, class TSolver>
TSolver& PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::GetLinearSolver() noexcept
{
    return *Interface::mpLinearSystemSolver;
}


// --------------------------------------------------------- //
// Template Instantiations
// --------------------------------------------------------- //


template class PMultigridBuilderAndSolver<TUblasSparseSpace<double>,
                                          TUblasDenseSpace<double>,
                                          LinearSolver<TUblasSparseSpace<double>,
                                                       TUblasDenseSpace<double>>>;


} // namespace Kratos
