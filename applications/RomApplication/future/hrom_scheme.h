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
#include <utility>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

#include "future/solving_strategies/schemes/implicit_scheme.h"

namespace Kratos::Future
{

///@name Kratos Classes
///@{

/**
 * @class HRomScheme
 * @ingroup RomApplication
 * @brief Implicit scheme that assembles only the entities of a hyper-reduced order model (HROM), scaled by their weights
 * @details This scheme behaves as the scheme it derives from, except in the build of the system: only the HROM elements and
 * conditions are assembled, with their local contributions multiplied by their HROM weight. The rest of the entities, as well as
 * the HROM entities with a null weight, are not computed. Note that the system keeps the full order size and that no row is
 * modified for the DOFs not touched by the HROM entities (the fixed DOFs are imposed through the basis in the projection).
 * @tparam TScheme The implicit scheme the HROM one derives from (e.g. StaticScheme)
 * @author Raul Bravo
 */
template<class TScheme>
class HRomScheme : public TScheme
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of HRomScheme
    KRATOS_CLASS_POINTER_DEFINITION(HRomScheme);

    /// Scheme the HROM one derives from
    using BaseType = TScheme;

    /// Implicit scheme base class, which defines the build helpers (the HROM scheme derives from schemes deriving from it)
    using ImplicitSchemeType = typename TScheme::BaseType;

    /// Matrix type definition
    using MatrixType = typename ImplicitSchemeType::MatrixType;

    /// Vector type definition
    using VectorType = typename ImplicitSchemeType::VectorType;

    /// Thread local storage used by the build helpers
    using TLSType = typename ImplicitSchemeType::TLSType;

    /// Index type definition
    using IndexType = std::size_t;

    ///@}
    ///@name Life Cycle
    ///@{

    /// @brief Constructor
    /// @param rModelPart Model part the scheme refers to
    /// @param ThisParameters Settings of the base scheme
    explicit HRomScheme(
        ModelPart& rModelPart,
        Parameters ThisParameters)
        : BaseType(rModelPart, ThisParameters)
    {
    }

    virtual ~HRomScheme() = default;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Sets the HROM entities and their weights
     * @param rElementIds Ids of the HROM elements
     * @param rElementWeights HROM weights of the elements (same order as the ids)
     * @param rConditionIds Ids of the HROM conditions
     * @param rConditionWeights HROM weights of the conditions (same order as the ids)
     */
    void SetHRomEntities(
        const std::vector<IndexType>& rElementIds,
        const std::vector<double>& rElementWeights,
        const std::vector<IndexType>& rConditionIds,
        const std::vector<double>& rConditionWeights)
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(rElementIds.size() != rElementWeights.size()) << "There are " << rElementIds.size() << " HROM element ids but "
            << rElementWeights.size() << " weights." << std::endl;
        KRATOS_ERROR_IF(rConditionIds.size() != rConditionWeights.size()) << "There are " << rConditionIds.size() << " HROM condition ids but "
            << rConditionWeights.size() << " weights." << std::endl;

        // Note that the entities with a null weight are skipped, as their contribution to the HROM is zero
        auto& r_model_part = this->GetModelPart();
        mHRomElements.clear();
        for (IndexType i = 0; i < rElementIds.size(); ++i) {
            if (rElementWeights[i] != 0.0) {
                mHRomElements.emplace_back(r_model_part.pGetElement(rElementIds[i]), rElementWeights[i]);
            }
        }
        mHRomConditions.clear();
        for (IndexType i = 0; i < rConditionIds.size(); ++i) {
            if (rConditionWeights[i] != 0.0) {
                mHRomConditions.emplace_back(r_model_part.pGetCondition(rConditionIds[i]), rConditionWeights[i]);
            }
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief Builds the LHS and the RHS with the HROM entities, scaled by their weights
     * @param rLHS The LHS matrix
     * @param rRHS The RHS vector
     */
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
            AssembleHRomEntities(mHRomElements, [&](auto& rEntity, const double Weight){
                if (this->CalculateLocalSystemContribution(rEntity, aux_tls, r_process_info)) {
                    aux_tls.LocalMatrix *= Weight;
                    aux_tls.LocalVector *= Weight;
                    rRHS.Assemble(aux_tls.LocalVector, aux_tls.LocalEqIds);
                    rLHS.Assemble(aux_tls.LocalMatrix, aux_tls.LocalEqIds);
                }
            });
            AssembleHRomEntities(mHRomConditions, [&](auto& rEntity, const double Weight){
                if (this->CalculateLocalSystemContribution(rEntity, aux_tls, r_process_info)) {
                    aux_tls.LocalMatrix *= Weight;
                    aux_tls.LocalVector *= Weight;
                    rRHS.Assemble(aux_tls.LocalVector, aux_tls.LocalEqIds);
                    rLHS.Assemble(aux_tls.LocalMatrix, aux_tls.LocalEqIds);
                }
            });
        }
        rRHS.FinalizeAssemble();
        rLHS.FinalizeAssemble();

        KRATOS_CATCH("")
    }

    /**
     * @brief Builds the RHS with the HROM entities, scaled by their weights
     * @param rRHS The RHS vector
     */
    void Build(VectorType& rRHS) override
    {
        KRATOS_TRY

        const auto& r_process_info = this->GetModelPart().GetProcessInfo();
        rRHS.BeginAssemble();
        #pragma omp parallel
        {
            TLSType aux_tls;
            const auto assemble_rhs = [&](auto& rEntity, const double Weight){
                if (this->CalculateRightHandSideContribution(rEntity, aux_tls, r_process_info)) {
                    aux_tls.LocalVector *= Weight;
                    rRHS.Assemble(aux_tls.LocalVector, aux_tls.LocalEqIds);
                }
            };
            AssembleHRomEntities(mHRomElements, assemble_rhs);
            AssembleHRomEntities(mHRomConditions, assemble_rhs);
        }
        rRHS.FinalizeAssemble();

        KRATOS_CATCH("")
    }

    /**
     * @brief Builds the LHS with the HROM entities, scaled by their weights
     * @param rLHS The LHS matrix
     */
    void Build(MatrixType& rLHS) override
    {
        KRATOS_TRY

        const auto& r_process_info = this->GetModelPart().GetProcessInfo();
        rLHS.BeginAssemble();
        #pragma omp parallel
        {
            TLSType aux_tls;
            const auto assemble_lhs = [&](auto& rEntity, const double Weight){
                if (this->CalculateLeftHandSideContribution(rEntity, aux_tls, r_process_info)) {
                    aux_tls.LocalMatrix *= Weight;
                    rLHS.Assemble(aux_tls.LocalMatrix, aux_tls.LocalEqIds);
                }
            };
            AssembleHRomEntities(mHRomElements, assemble_lhs);
            AssembleHRomEntities(mHRomConditions, assemble_lhs);
        }
        rLHS.FinalizeAssemble();

        KRATOS_CATCH("")
    }

    ///@}
    ///@name Access
    ///@{

    /// Number of HROM elements with a non-null weight
    IndexType NumberOfHRomElements() const
    {
        return mHRomElements.size();
    }

    /// Number of HROM conditions with a non-null weight
    IndexType NumberOfHRomConditions() const
    {
        return mHRomConditions.size();
    }

    ///@}
    ///@name Input and output
    ///@{

    std::string Info() const override
    {
        return "HRomScheme";
    }

    ///@}
private:
    ///@name Member Variables
    ///@{

    std::vector<std::pair<Element::Pointer, double>> mHRomElements; /// HROM elements with a non-null weight, and their weights

    std::vector<std::pair<Condition::Pointer, double>> mHRomConditions; /// HROM conditions with a non-null weight, and their weights

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Applies the provided assembly function to the given HROM entities (to be called inside a parallel region)
     * @param rEntities HROM entities and their weights
     * @param rAssemble Function computing and assembling the weighted contribution of an entity
     */
    template<class TEntitiesType, class TAssembleFunction>
    static void AssembleHRomEntities(
        TEntitiesType& rEntities,
        const TAssembleFunction& rAssemble)
    {
        const int n_entities = static_cast<int>(rEntities.size());
        #pragma omp for schedule(guided, 64)
        for (int k = 0; k < n_entities; ++k) {
            rAssemble(*(rEntities[k].first), rEntities[k].second);
        }
    }

    ///@}
}; // Class HRomScheme

///@}

} // namespace Kratos::Future
