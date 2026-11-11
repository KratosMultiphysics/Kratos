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

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "spaces/ublas_space.h"
#include "solving_strategies/schemes/scheme.h"

// Application includes
#include "future/containers/define_linear_algebra_serial.h"
#include "future/solving_strategies/schemes/static_scheme.h"
#include "future/solving_strategies/strategies/implicit_strategy_data.h"

namespace Kratos::Future
{

/**
 * @class LegacySchemeAdapter
 * @ingroup RomApplication
 * @brief Future implicit scheme whose time integration is the one of a legacy scheme
 * @details The DOFs, the system and its assembly are the ones of the Future scheme (no builder and solver). The contribution of each
 * element and condition, the prediction, the update and the operations of each step and iteration are forwarded to the legacy scheme,
 * e.g. the one a legacy solver creates for its problem (Bossak, BDF, the rotation scheme of the fluids...).
 * Note that the master-slave constraints are not considered in the prediction and update.
 * @tparam TLinearAlgebra The struct containing the linear algebra types
 */
template<class TLinearAlgebra = SerialLinearAlgebraTraits>
class LegacySchemeAdapter : public StaticScheme<TLinearAlgebra>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(LegacySchemeAdapter);

    using BaseType = StaticScheme<TLinearAlgebra>;

    using MatrixType = typename TLinearAlgebra::MatrixType;

    using VectorType = typename TLinearAlgebra::VectorType;

    using TLSType = typename BaseType::TLSType;

    using IndexType = std::size_t;

    using LegacySparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;

    using LegacyLocalSpaceType = UblasSpace<double, Matrix, Vector>;

    using LegacySchemeType = Scheme<LegacySparseSpaceType, LegacyLocalSpaceType>;

    ///@}
    ///@name Life Cycle
    ///@{

    explicit LegacySchemeAdapter(
        ModelPart& rModelPart,
        Parameters ThisParameters,
        typename LegacySchemeType::Pointer pLegacyScheme)
        : BaseType(rModelPart, ThisParameters)
        , mpLegacyScheme(pLegacyScheme)
    {
        KRATOS_ERROR_IF_NOT(mpLegacyScheme) << "No legacy scheme provided." << std::endl;
    }

    virtual ~LegacySchemeAdapter() = default;

    ///@}
    ///@name Operations
    ///@{

    void Initialize(ImplicitStrategyData<TLinearAlgebra>& rImplicitStrategyData) override
    {
        KRATOS_TRY

        // DOFs, equation ids and system arrays: the Future ones
        this->InitializeLinearSystem(rImplicitStrategyData);

        // Scheme, elements and conditions: as the legacy strategies do
        auto& r_model_part = this->GetModelPart();
        if (!mpLegacyScheme->SchemeIsInitialized()) {
            mpLegacyScheme->Initialize(r_model_part);
        }
        if (!mpLegacyScheme->ElementsAreInitialized()) {
            mpLegacyScheme->InitializeElements(r_model_part);
        }
        if (!mpLegacyScheme->ConditionsAreInitialized()) {
            mpLegacyScheme->InitializeConditions(r_model_part);
        }

        KRATOS_CATCH("")
    }

    void InitializeSolutionStep(ImplicitStrategyData<TLinearAlgebra>& rImplicitStrategyData) override
    {
        ResizeLegacyVectors(rImplicitStrategyData);
        mpLegacyScheme->InitializeSolutionStep(this->GetModelPart(), mA, mDx, mb);
    }

    void Predict(ImplicitStrategyData<TLinearAlgebra>& rImplicitStrategyData) override
    {
        ResizeLegacyVectors(rImplicitStrategyData);
        mpLegacyScheme->Predict(this->GetModelPart(), *(rImplicitStrategyData.pGetDofSet()), mA, mDx, mb);
        if (this->GetMoveMesh()) {
            this->MoveMesh();
        }
    }

    void InitializeNonLinIteration(ImplicitStrategyData<TLinearAlgebra>& rImplicitStrategyData) override
    {
        mpLegacyScheme->InitializeNonLinIteration(this->GetModelPart(), mA, mDx, mb);
    }

    void FinalizeNonLinIteration(ImplicitStrategyData<TLinearAlgebra>& rImplicitStrategyData) override
    {
        mpLegacyScheme->FinalizeNonLinIteration(this->GetModelPart(), mA, mDx, mb);
    }

    void FinalizeSolutionStep(ImplicitStrategyData<TLinearAlgebra>& rImplicitStrategyData) override
    {
        mpLegacyScheme->FinalizeSolutionStep(this->GetModelPart(), mA, mDx, mb);
        mpLegacyScheme->Clean();
    }

    /// The increment of the Future linear system (Dx) applied by the legacy scheme, which also updates its time derivatives
    void Update(ImplicitStrategyData<TLinearAlgebra>& rImplicitStrategyData) override
    {
        KRATOS_TRY

        ResizeLegacyVectors(rImplicitStrategyData);
        this->CalculateUpdateVector(rImplicitStrategyData);
        const auto& r_dx = *(rImplicitStrategyData.pGetLinearSystem()->pGetVector(LinearSystemTags::DenseVectorTag::Dx));
        IndexPartition<IndexType>(mDx.size()).for_each([&](IndexType i){
            mDx[i] = r_dx[i];
        });
        mpLegacyScheme->Update(this->GetModelPart(), *(rImplicitStrategyData.pGetDofSet()), mA, mDx, mb);
        if (this->GetMoveMesh()) {
            this->MoveMesh();
        }

        KRATOS_CATCH("")
    }

    void Build(
        MatrixType& rLHS,
        VectorType& rRHS) override
    {
        KRATOS_TRY

        const auto& r_process_info = this->GetModelPart().GetProcessInfo();
        rRHS.BeginAssemble();
        rLHS.BeginAssemble();
        #pragma omp parallel
        {
            TLSType aux_tls;
            const auto assemble = [&](auto& rEntity){
                if (this->CalculateLocalSystemContribution(rEntity, aux_tls, r_process_info)) {
                    rRHS.Assemble(aux_tls.LocalVector, aux_tls.LocalEqIds);
                    rLHS.Assemble(aux_tls.LocalMatrix, aux_tls.LocalEqIds);
                }
            };
            ForEachEntity(this->GetModelPart().Elements(), assemble);
            ForEachEntity(this->GetModelPart().Conditions(), assemble);
        }
        rRHS.FinalizeAssemble();
        rLHS.FinalizeAssemble();

        KRATOS_CATCH("")
    }

    void Build(VectorType& rRHS) override
    {
        KRATOS_TRY

        const auto& r_process_info = this->GetModelPart().GetProcessInfo();
        rRHS.BeginAssemble();
        #pragma omp parallel
        {
            TLSType aux_tls;
            const auto assemble = [&](auto& rEntity){
                if (this->CalculateRightHandSideContribution(rEntity, aux_tls, r_process_info)) {
                    rRHS.Assemble(aux_tls.LocalVector, aux_tls.LocalEqIds);
                }
            };
            ForEachEntity(this->GetModelPart().Elements(), assemble);
            ForEachEntity(this->GetModelPart().Conditions(), assemble);
        }
        rRHS.FinalizeAssemble();

        KRATOS_CATCH("")
    }

    void Build(MatrixType& rLHS) override
    {
        KRATOS_TRY

        const auto& r_process_info = this->GetModelPart().GetProcessInfo();
        rLHS.BeginAssemble();
        #pragma omp parallel
        {
            TLSType aux_tls;
            const auto assemble = [&](auto& rEntity){
                if (this->CalculateLeftHandSideContribution(rEntity, aux_tls, r_process_info)) {
                    rLHS.Assemble(aux_tls.LocalMatrix, aux_tls.LocalEqIds);
                }
            };
            ForEachEntity(this->GetModelPart().Elements(), assemble);
            ForEachEntity(this->GetModelPart().Conditions(), assemble);
        }
        rLHS.FinalizeAssemble();

        KRATOS_CATCH("")
    }

    typename LegacySchemeType::Pointer pGetLegacyScheme()
    {
        return mpLegacyScheme;
    }

    std::string Info() const override
    {
        return "LegacySchemeAdapter";
    }

    ///@}
protected:
    ///@name Protected Operations
    ///@{

    // Contributions of an entity through the legacy scheme (they hide the ones of the implicit scheme, which ask the entity alone, so
    // that the schemes deriving from this one, e.g. the HROM one, get the contributions of the legacy scheme too)
    template<class TEntityType>
    bool CalculateLocalSystemContribution(
        TEntityType& rEntity,
        TLSType& rTLS,
        const ProcessInfo& rProcessInfo)
    {
        if (rEntity.IsActive()) {
            mpLegacyScheme->CalculateSystemContributions(rEntity, rTLS.LocalMatrix, rTLS.LocalVector, rTLS.LocalEqIds, rProcessInfo);
            return true;
        }
        return false;
    }

    template<class TEntityType>
    bool CalculateRightHandSideContribution(
        TEntityType& rEntity,
        TLSType& rTLS,
        const ProcessInfo& rProcessInfo)
    {
        if (rEntity.IsActive()) {
            mpLegacyScheme->CalculateRHSContribution(rEntity, rTLS.LocalVector, rTLS.LocalEqIds, rProcessInfo);
            return true;
        }
        return false;
    }

    template<class TEntityType>
    bool CalculateLeftHandSideContribution(
        TEntityType& rEntity,
        TLSType& rTLS,
        const ProcessInfo& rProcessInfo)
    {
        if (rEntity.IsActive()) {
            mpLegacyScheme->CalculateLHSContribution(rEntity, rTLS.LocalMatrix, rTLS.LocalEqIds, rProcessInfo);
            return true;
        }
        return false;
    }

    ///@}
private:
    ///@name Member Variables
    ///@{

    typename LegacySchemeType::Pointer mpLegacyScheme;

    // Arguments of the legacy interface: the matrix and the RHS are not built (empty matrix, null RHS), Dx holds the last increment
    CompressedMatrix mA;
    Vector mDx;
    Vector mb;

    ///@}
    ///@name Private Operations
    ///@{

    void ResizeLegacyVectors(ImplicitStrategyData<TLinearAlgebra>& rImplicitStrategyData)
    {
        const IndexType system_size = rImplicitStrategyData.pGetDofSet()->size();
        if (mDx.size() != system_size) {
            mDx = ZeroVector(system_size);
            mb = ZeroVector(system_size);
        }
    }

    template<class TEntitiesContainerType, class TFunction>
    static void ForEachEntity(
        TEntitiesContainerType& rEntities,
        const TFunction& rFunction)
    {
        const int n_entities = static_cast<int>(rEntities.size());
        const auto it_begin = rEntities.begin();
        #pragma omp for schedule(guided, 64)
        for (int k = 0; k < n_entities; ++k) {
            rFunction(*(it_begin + k));
        }
    }

    ///@}
}; // Class LegacySchemeAdapter

} // namespace Kratos::Future
