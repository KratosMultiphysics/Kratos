//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Raul Bravo
//

#pragma once

// System includes
#include <limits>
#include <vector>

// External includes
#include <Eigen/Core>
#include <Eigen/Dense>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"

// Future extensions
#include "future/containers/define_linear_algebra_serial.h"
#include "future/containers/linear_system_tags.h"
#include "future/solving_strategies/schemes/implicit_scheme.h"
#include "future/solving_strategies/strategies/implicit_strategy_data.h"

namespace Kratos::Future
{

///@addtogroup RomApplication
///@{

///@name Kratos Classes
///@{

/**
 * @class RomProjector
 * @ingroup RomApplication
 * @brief Builds, projects and solves the reduced system of a projection-based ROM using the Future schemes
 * @details The class is decoder-agnostic. The full order system is built through a Future implicit scheme
 * and projected onto the basis provided by the caller in each call, which is the Jacobian of the decoder.
 * The rows of the basis follow the effective DOF set (row i belongs to the DOF with effective equation id i).
 * @tparam TLinearAlgebra The struct containing the linear algebra types
 * @author Raul Bravo
 */
template<class TLinearAlgebra = SerialLinearAlgebraTraits>
class RomProjector
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RomProjector
    KRATOS_CLASS_POINTER_DEFINITION(RomProjector);

    /// Index type definition
    using IndexType = typename TLinearAlgebra::IndexType;

    /// Matrix type definition
    using MatrixType = typename TLinearAlgebra::MatrixType;

    /// Vector type definition
    using VectorType = typename TLinearAlgebra::VectorType;

    /// Scheme type definition
    using SchemeType = ImplicitScheme<TLinearAlgebra>;

    /// Strategy data container type definition
    using StrategyDataType = ImplicitStrategyData<TLinearAlgebra>;

    /// DoF type definition
    using DofType = Dof<double>;

    /// Dense types of the basis and of the reduced system
    using EigenDynamicMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using EigenDynamicVector = Eigen::Matrix<double, Eigen::Dynamic, 1>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// @brief Constructor
    /// @param pScheme Scheme used to build the full order system
    /// @param pStrategyData Container with the linear system arrays and DOF sets, already initialized by the scheme
    RomProjector(
        typename SchemeType::Pointer pScheme,
        typename StrategyDataType::Pointer pStrategyData)
        : mpScheme(pScheme)
        , mpStrategyData(pStrategyData)
    {
    }

    virtual ~RomProjector() = default;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Builds the full order system with the constraints applied
     * This method sets the system arrays to zero, builds them from the current database and applies the
     * master-slave and Dirichlet constraints, leaving the result in the effective linear system
     */
    void BuildEffectiveSystem()
    {
        KRATOS_TRY

        auto p_linear_system = mpStrategyData->pGetLinearSystem();
        auto& r_lhs = *(p_linear_system->pGetMatrix(LinearSystemTags::SparseMatrixTag::LHS));
        auto& r_rhs = *(p_linear_system->pGetVector(LinearSystemTags::DenseVectorTag::RHS));
        auto& r_dx = *(p_linear_system->pGetVector(LinearSystemTags::DenseVectorTag::Dx));
        r_lhs.SetValue(0.0);
        r_rhs.SetValue(0.0);
        r_dx.SetValue(0.0);

        mpScheme->Build(r_lhs, r_rhs);
        mpScheme->BuildLinearSystemConstraints(*mpStrategyData);
        mpScheme->ApplyLinearSystemConstraints(*mpStrategyData);

        KRATOS_CATCH("")
    }

    /**
     * @brief Galerkin projection of the effective system
     * This method computes Phi^T A Phi and Phi^T b, taking as zero the rows of Phi of the fixed DOFs
     * @param rPhi Basis (Jacobian of the decoder), with one row per effective DOF and one column per ROM DOF
     */
    void Project(const Eigen::Ref<const EigenDynamicMatrix>& rPhi)
    {
        KRATOS_TRY

        const auto p_eff_linear_system = mpStrategyData->pGetEffectiveLinearSystem();
        const auto& r_eff_lhs = *(p_eff_linear_system->pGetMatrix(LinearSystemTags::SparseMatrixTag::LHS));
        const auto& r_eff_rhs = *(p_eff_linear_system->pGetVector(LinearSystemTags::DenseVectorTag::RHS));
        const auto& r_eff_dof_set = *(mpStrategyData->pGetEffectiveDofSet());

        const IndexType system_size = r_eff_lhs.size1();
        const IndexType n_rom_dofs = rPhi.cols();
        KRATOS_ERROR_IF(static_cast<IndexType>(rPhi.rows()) != system_size) << "The basis has " << rPhi.rows()
            << " rows but the effective system size is " << system_size << "." << std::endl;

        // Set the free DOFs vector (0 means fixed / 1 means free)
        std::vector<uint8_t> free_dofs_vector(system_size, 1);
        const auto dof_begin = r_eff_dof_set.begin();
        IndexPartition<IndexType>(r_eff_dof_set.size()).for_each([&](IndexType Index){
            const auto it_dof = dof_begin + Index;
            if (it_dof->IsFixed()) {
                free_dofs_vector[it_dof->EffectiveEquationId()] = 0;
            }
        });

        // Compute A*Phi and the RHS without the fixed DOFs
        // Note that as the rows of the fixed DOFs are zero in A*Phi there is no need to modify Phi in the left product
        if (static_cast<IndexType>(mLhsTimesPhi.rows()) != system_size || static_cast<IndexType>(mLhsTimesPhi.cols()) != n_rom_dofs) {
            mLhsTimesPhi.resize(system_size, n_rom_dofs);
        }
        EigenDynamicVector free_rhs(system_size);
        const auto& r_row_indices = r_eff_lhs.index1_data();
        const auto& r_col_indices = r_eff_lhs.index2_data();
        const auto& r_values = r_eff_lhs.value_data();
        IndexPartition<IndexType>(system_size).for_each([&](IndexType i){
            auto lhs_times_phi_row = mLhsTimesPhi.row(i);
            lhs_times_phi_row.setZero();
            if (free_dofs_vector[i]) {
                for (IndexType k = r_row_indices[i]; k < r_row_indices[i+1]; ++k) {
                    const IndexType j = r_col_indices[k];
                    if (free_dofs_vector[j]) {
                        lhs_times_phi_row.noalias() += r_values[k] * rPhi.row(j);
                    }
                }
                free_rhs[i] = r_eff_rhs[i];
            } else {
                free_rhs[i] = 0.0;
            }
        });

        mReducedLhs.noalias() = rPhi.transpose() * mLhsTimesPhi;
        mReducedRhs.noalias() = rPhi.transpose() * free_rhs;

        KRATOS_CATCH("")
    }

    /**
     * @brief Solves the reduced system obtained in the last projection
     * @return EigenDynamicVector The increment of the reduced coordinates
     */
    EigenDynamicVector SolveReduced() const
    {
        KRATOS_TRY

        return mReducedLhs.colPivHouseholderQr().solve(mReducedRhs);

        KRATOS_CATCH("")
    }

    /**
     * @brief Sets the provided solution in the free DOFs
     * The solution is imposed through its increment with respect to the current one, so that the update
     * of the constraints and of the mesh is the one of the scheme
     * @param rSolution Solution vector, with one entry per effective DOF
     */
    void SetSolution(const Eigen::Ref<const EigenDynamicVector>& rSolution)
    {
        KRATOS_TRY

        auto& r_eff_dx = *(mpStrategyData->pGetEffectiveLinearSystem()->pGetVector(LinearSystemTags::DenseVectorTag::Dx));
        auto& r_eff_dof_set = *(mpStrategyData->pGetEffectiveDofSet());
        KRATOS_ERROR_IF(static_cast<IndexType>(rSolution.size()) != r_eff_dof_set.size()) << "The solution has " << rSolution.size()
            << " entries but there are " << r_eff_dof_set.size() << " effective DOFs." << std::endl;

        block_for_each(r_eff_dof_set, [&](const DofType& rDof){
            const IndexType eff_eq_id = rDof.EffectiveEquationId();
            r_eff_dx[eff_eq_id] = rSolution[eff_eq_id] - rDof.GetSolutionStepValue();
        });

        mpScheme->Update(*mpStrategyData);

        KRATOS_CATCH("")
    }

    /**
     * @brief Projects the residual of each element and condition onto the provided basis
     * Row k is Phi_k^T r_k, with r_k the local RHS of the k-th entity of the model part of the scheme (the elements first, then
     * the conditions) and Phi_k the rows of the basis of its DOFs, taken as zero for the fixed DOFs as in Project. The rows of the
     * inactive entities are zero. These are the residual snapshots of the HROM training (see HRomTrainingUtility).
     * Note that the master-slave constraints are not supported.
     * @param rPhi Basis (Jacobian of the decoder), with one row per effective DOF and one column per ROM DOF
     * @return EigenDynamicMatrix The projected residuals, with one row per entity and one column per ROM DOF
     */
    EigenDynamicMatrix ProjectEntityResiduals(const Eigen::Ref<const EigenDynamicMatrix>& rPhi) const
    {
        KRATOS_TRY

        auto& r_model_part = mpScheme->GetModelPart();
        const auto& r_dof_set = *(mpStrategyData->pGetDofSet());
        const auto& r_eff_dof_set = *(mpStrategyData->pGetEffectiveDofSet());
        KRATOS_ERROR_IF(static_cast<IndexType>(rPhi.rows()) != r_eff_dof_set.size()) << "The basis has " << rPhi.rows()
            << " rows but there are " << r_eff_dof_set.size() << " effective DOFs." << std::endl;

        // Row of the basis of each DOF, by its equation id (the one of the local systems), or -1 for the fixed DOFs
        constexpr auto non_effective_id = std::numeric_limits<typename DofType::EquationIdType>::max();
        std::vector<int> basis_rows(r_dof_set.size());
        const auto dof_begin = r_dof_set.begin();
        IndexPartition<IndexType>(r_dof_set.size()).for_each([&](IndexType Index){
            const auto it_dof = dof_begin + Index;
            KRATOS_ERROR_IF(it_dof->EffectiveEquationId() == non_effective_id) << "The DOF " << it_dof->GetVariable().Name() << " of node "
                << it_dof->Id() << " is not an effective one. The master-slave constraints are not supported." << std::endl;
            basis_rows[it_dof->EquationId()] = it_dof->IsFixed() ? -1 : static_cast<int>(it_dof->EffectiveEquationId());
        });

        const IndexType n_elements = r_model_part.NumberOfElements();
        const IndexType n_conditions = r_model_part.NumberOfConditions();
        EigenDynamicMatrix projected_residuals = EigenDynamicMatrix::Zero(n_elements + n_conditions, rPhi.cols());
        const auto& r_process_info = r_model_part.GetProcessInfo();

        struct TLS
        {
            Vector LocalVector;
            Element::EquationIdVectorType LocalEqIds;
        };
        const auto project_entity_residual = [&](auto& rEntity, const IndexType Row, TLS& rTLS){
            if (rEntity.IsActive()) {
                rEntity.CalculateRightHandSide(rTLS.LocalVector, r_process_info);
                rEntity.EquationIdVector(rTLS.LocalEqIds, r_process_info);
                auto projected_residual = projected_residuals.row(Row);
                for (IndexType i = 0; i < rTLS.LocalEqIds.size(); ++i) {
                    const int basis_row = basis_rows[rTLS.LocalEqIds[i]];
                    if (basis_row >= 0) {
                        projected_residual.noalias() += rTLS.LocalVector[i] * rPhi.row(basis_row);
                    }
                }
            }
        };
        const auto elements_begin = r_model_part.ElementsBegin();
        IndexPartition<IndexType>(n_elements).for_each(TLS(), [&](IndexType Index, TLS& rTLS){
            project_entity_residual(*(elements_begin + Index), Index, rTLS);
        });
        const auto conditions_begin = r_model_part.ConditionsBegin();
        IndexPartition<IndexType>(n_conditions).for_each(TLS(), [&](IndexType Index, TLS& rTLS){
            project_entity_residual(*(conditions_begin + Index), n_elements + Index, rTLS);
        });

        return projected_residuals;

        KRATOS_CATCH("")
    }

    ///@}
    ///@name Access
    ///@{

    const EigenDynamicMatrix& GetReducedLhs() const
    {
        return mReducedLhs;
    }

    const EigenDynamicVector& GetReducedRhs() const
    {
        return mReducedRhs;
    }

    ///@}
private:
    ///@name Member Variables
    ///@{

    typename SchemeType::Pointer mpScheme;

    typename StrategyDataType::Pointer mpStrategyData;

    EigenDynamicMatrix mLhsTimesPhi; // Product of the effective LHS and the basis

    EigenDynamicMatrix mReducedLhs;

    EigenDynamicVector mReducedRhs;

    ///@}
}; // Class RomProjector

///@}
///@} addtogroup block

} // namespace Kratos::Future
