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
//  PROTOTYPE (study work/future-rom/strategy-plugin): not part of any PR
//

#pragma once

// System includes
#include <limits>
#include <utility>
#include <vector>

// External includes
#include <Eigen/Core>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

// Application includes
#include "future/solving_strategies/strategies/implicit_strategy_data.h"

namespace Kratos::Future
{

/**
 * @class EntityResidualsScheme
 * @ingroup RomApplication
 * @brief Any implicit scheme, plus the projection of the residual of each element and condition (HROM training)
 * @details The residual of each entity is the one of the scheme it derives from (its CalculateLocalSystemContribution), so it
 * includes whatever the scheme adds to the entity (e.g. the time integration terms of a dynamic scheme, or the ones of the legacy
 * scheme of a LegacySchemeAdapter). Everything else is the base scheme.
 * @tparam TScheme The implicit scheme it derives from (e.g. StaticScheme, LegacySchemeAdapter)
 */
template<class TScheme>
class EntityResidualsScheme : public TScheme
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(EntityResidualsScheme);

    using BaseType = TScheme;

    /// Thread local storage of the build helpers of the base scheme (as in HRomScheme)
    using TLSType = typename TScheme::BaseType::TLSType;

    using IndexType = std::size_t;

    using EigenDynamicMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    template<class... TArguments>
    explicit EntityResidualsScheme(
        ModelPart& rModelPart,
        Parameters ThisParameters,
        TArguments&&... rArguments)
        : BaseType(rModelPart, ThisParameters, std::forward<TArguments>(rArguments)...)
    {
    }

    virtual ~EntityResidualsScheme() = default;

    /**
     * @brief Projects the residual of each element and condition onto the provided basis
     * Row k is Phi_k^T r_k, with r_k the residual of the k-th entity of the model part (the elements first, then the conditions) as
     * contributed through the scheme, and Phi_k the rows of the basis of its DOFs, taken as zero for the fixed DOFs. The rows of the
     * inactive entities are zero. Note that the master-slave constraints are not supported.
     * @param rPhi Basis (Jacobian of the decoder), with one row per effective DOF and one column per ROM DOF
     * @param rImplicitStrategyData Strategy data initialized by this scheme
     */
    template<class TStrategyData>
    EigenDynamicMatrix ProjectEntityResiduals(
        const Eigen::Ref<const EigenDynamicMatrix>& rPhi,
        TStrategyData& rImplicitStrategyData)
    {
        KRATOS_TRY

        auto& r_model_part = this->GetModelPart();
        const auto& r_dof_set = *(rImplicitStrategyData.pGetDofSet());
        const auto& r_eff_dof_set = *(rImplicitStrategyData.pGetEffectiveDofSet());
        KRATOS_ERROR_IF(static_cast<IndexType>(rPhi.rows()) != r_eff_dof_set.size()) << "The basis has " << rPhi.rows()
            << " rows but there are " << r_eff_dof_set.size() << " effective DOFs." << std::endl;

        // Row of the basis of each DOF, by its equation id (the one of the local systems), or -1 for the fixed DOFs
        std::vector<int> basis_rows(r_dof_set.size());
        const auto dof_begin = r_dof_set.begin();
        IndexPartition<IndexType>(r_dof_set.size()).for_each([&](IndexType Index){
            const auto it_dof = dof_begin + Index;
            basis_rows[it_dof->EquationId()] = it_dof->IsFixed() ? -1 : static_cast<int>(it_dof->EffectiveEquationId());
        });

        const int n_elements = static_cast<int>(r_model_part.NumberOfElements());
        const int n_conditions = static_cast<int>(r_model_part.NumberOfConditions());
        EigenDynamicMatrix projected_residuals = EigenDynamicMatrix::Zero(n_elements + n_conditions, rPhi.cols());
        const auto& r_process_info = r_model_part.GetProcessInfo();
        const auto elements_begin = r_model_part.ElementsBegin();
        const auto conditions_begin = r_model_part.ConditionsBegin();

        #pragma omp parallel
        {
            TLSType aux_tls;
            const auto project_entity_residual = [&](auto& rEntity, const IndexType Row){
                // The contribution of the entity through the scheme (the local system: some entities do not compute the RHS alone)
                if (this->CalculateLocalSystemContribution(rEntity, aux_tls, r_process_info)) {
                    auto projected_residual = projected_residuals.row(Row);
                    for (IndexType i = 0; i < aux_tls.LocalEqIds.size(); ++i) {
                        const int basis_row = basis_rows[aux_tls.LocalEqIds[i]];
                        if (basis_row >= 0) {
                            projected_residual.noalias() += aux_tls.LocalVector[i] * rPhi.row(basis_row);
                        }
                    }
                }
            };
            #pragma omp for schedule(guided, 64)
            for (int k = 0; k < n_elements; ++k) {
                project_entity_residual(*(elements_begin + k), k);
            }
            #pragma omp for schedule(guided, 64)
            for (int k = 0; k < n_conditions; ++k) {
                project_entity_residual(*(conditions_begin + k), n_elements + k);
            }
        }

        return projected_residuals;

        KRATOS_CATCH("")
    }
}; // Class EntityResidualsScheme

} // namespace Kratos::Future
