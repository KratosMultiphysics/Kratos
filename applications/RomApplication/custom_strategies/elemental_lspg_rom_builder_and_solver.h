//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:   Raul Bravo
//

#pragma once

/* System includes */
#include <vector>

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/lock_object.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "custom_strategies/rom_builder_and_solver.h"
#include "utilities/builtin_timer.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

/* Application includes */
#include "rom_application_variables.h"
#include "custom_utilities/rom_auxiliary_utilities.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/**
 * @class ElementalLeastSquaresPetrovGalerkinROMBuilderAndSolver
 * @ingroup RomApplication
 * @brief Least-Squares Petrov-Galerkin (LSPG) ROM builder and solver that never assembles the global system matrix.
 * @details The reduced system (J*Phi)^T (J*Phi) dq = (J*Phi)^T r is obtained in two loops over the entities:
 *  1. The test functions W = J*Phi are assembled row by row from the elemental products J_e*Phi_e. Row i of W needs
 *     all the entities sharing the DOF i, so in a HROM this loop runs over the selected entities and their nodal
 *     neighbours (the complementary mesh), without weights.
 *  2. The reduced system is the sum of W_e^T J_e Phi_e and W_e^T r_e, where W_e are the rows of W of the DOFs of the
 *     entity. In a HROM this loop runs over the selected entities only, with their weights.
 * It is equivalent to the LeastSquaresPetrovGalerkinROMBuilderAndSolver with the 'normal_equations' solving technique.
 * @tparam TSparseSpace Sparse space template type
 * @tparam TDenseSpace Local space template type
 * @tparam TLinearSolver Linear solver template type
 */
template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
class ElementalLeastSquaresPetrovGalerkinROMBuilderAndSolver : public ROMBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:

    ///@name Type Definitions
    ///@{

    // Class pointer definition
    KRATOS_CLASS_POINTER_DEFINITION(ElementalLeastSquaresPetrovGalerkinROMBuilderAndSolver);

    // The size_t types
    using SizeType = std::size_t;
    using IndexType = std::size_t;

    /// The definition of the current class
    using ClassType = ElementalLeastSquaresPetrovGalerkinROMBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>;

    /// Definitions from the base classes
    using BaseBuilderAndSolverType = BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>;
    using BaseType = ROMBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>;
    using TSchemeType = typename BaseType::TSchemeType;
    using ElementsArrayType = typename BaseType::ElementsArrayType;
    using ConditionsArrayType = typename BaseType::ConditionsArrayType;
    using LocalSystemVectorType = typename BaseType::LocalSystemVectorType;
    using LocalSystemMatrixType = typename BaseType::LocalSystemMatrixType;
    using RomSystemMatrixType = typename BaseType::RomSystemMatrixType;
    using RomSystemVectorType = typename BaseType::RomSystemVectorType;
    using EquationIdVectorType = typename BaseType::EquationIdVectorType;
    using DofsVectorType = typename BaseType::DofsVectorType;

    ///@}
    ///@name Life cycle
    ///@{

    explicit ElementalLeastSquaresPetrovGalerkinROMBuilderAndSolver(
        typename TLinearSolver::Pointer pNewLinearSystemSolver,
        Parameters ThisParameters) : BaseType(pNewLinearSystemSolver)
    {
        // Validate and assign defaults
        Parameters this_parameters_copy = ThisParameters.Clone();
        this_parameters_copy = this->ValidateAndAssignParameters(this_parameters_copy, this->GetDefaultParameters());
        this->AssignSettings(this_parameters_copy);
    }

    ~ElementalLeastSquaresPetrovGalerkinROMBuilderAndSolver() = default;

    ///@}
    ///@name Operations
    ///@{

    typename BaseBuilderAndSolverType::Pointer Create(
        typename TLinearSolver::Pointer pNewLinearSystemSolver,
        Parameters ThisParameters) const override
    {
        return Kratos::make_shared<ClassType>(pNewLinearSystemSolver, ThisParameters);
    }

    void SetUpDofSet(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart) override
    {
        BaseType::SetUpDofSet(pScheme, rModelPart);

        // The test functions of the selected entities need the contributions of their neighbours
        if (this->mHromSimulation) {
            RomAuxiliaryUtilities::GetEntitiesAndNodalNeighbours(rModelPart, this->mSelectedElements, this->mSelectedConditions, mComplementaryElements, mComplementaryConditions);
        }
    }

    /**
     * @brief Computes the test functions of the LSPG ROM (the Jacobian times the right basis) entity by entity.
     * @details The rows of the fixed DOFs are zero. In a HROM only the rows of the DOFs of the selected entities are complete.
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @return The matrix J*Phi (number of DOFs x number of ROM modes)
     */
    Matrix CalculateJPhi(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart)
    {
        KRATOS_ERROR_IF(!pScheme) << "No scheme provided!" << std::endl;
        BuildJPhi(*pScheme, rModelPart);
        return mJPhi;
    }

    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name" : "elemental_lspg_rom_builder_and_solver",
            "nodal_unknowns" : [],
            "number_of_rom_dofs" : 10,
            "rom_bns_settings": {
                "train_petrov_galerkin" : false,
                "solving_technique" : "normal_equations",
                "basis_strategy" : "residuals",
                "monotonicity_preserving" : false
            },
            "weight_vector_index": 0,
            "number_of_hrom_sets": 1
        })");
        default_parameters.AddMissingParameters(BaseType::GetDefaultParameters());

        return default_parameters;
    }

    static std::string Name()
    {
        return "elemental_lspg_rom_builder_and_solver";
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const override
    {
        return "ElementalLeastSquaresPetrovGalerkinROMBuilderAndSolver";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const override
    {
        rOStream << Info();
    }

    ///@}
protected:

    ///@name Protected operations
    ///@{

    void AssignSettings(const Parameters ThisParameters) override
    {
        BaseType::AssignSettings(ThisParameters);

        // Note that the inner settings might be only partially provided
        const auto& r_bns_settings = ThisParameters["rom_bns_settings"];
        const std::string solving_technique = r_bns_settings.Has("solving_technique") ? r_bns_settings["solving_technique"].GetString() : "normal_equations";
        KRATOS_ERROR_IF(solving_technique != "normal_equations") << "The elemental LSPG ROM only supports the 'normal_equations' solving technique. The provided one is '" << solving_technique << "'. Use the 'global' assembling strategy for other techniques." << std::endl;
        const auto is_requested = [&r_bns_settings](const std::string& rKey) {
            return r_bns_settings.Has(rKey) && r_bns_settings[rKey].GetBool();
        };
        KRATOS_ERROR_IF(is_requested("train_petrov_galerkin")) << "The Petrov-Galerkin training is not available with the elemental LSPG ROM. Use the 'global' assembling strategy." << std::endl;
        KRATOS_ERROR_IF(is_requested("monotonicity_preserving")) << "'monotonicity_preserving' is not available with the elemental LSPG ROM. Use the 'global' assembling strategy." << std::endl;
    }

    /**
     * Thread local storage of the assembly
     */
    struct LSPGAssemblyTLS
    {
        LSPGAssemblyTLS(SizeType NRomModes)
            : romA(ZeroMatrix(NRomModes, NRomModes)),
              romB(ZeroVector(NRomModes))
        { }
        LSPGAssemblyTLS() = delete;

        Matrix phiE = {};                // Elemental Phi
        Matrix psiE = {};                // Elemental test functions (rows of J*Phi)
        LocalSystemMatrixType lhs = {};  // Elemental LHS
        LocalSystemVectorType rhs = {};  // Elemental RHS
        EquationIdVectorType eq_id = {}; // Elemental equation ID vector
        DofsVectorType dofs = {};        // Elemental dof vector
        RomSystemMatrixType romA;        // reduced LHS
        RomSystemVectorType romB;        // reduced RHS
        RomSystemMatrixType aux = {};    // Auxiliary: LHS * phi
    };

    /**
     * Builds the reduced system of equations
     */
    void BuildROM(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        RomSystemMatrixType& rA,
        RomSystemVectorType& rb) override
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(!pScheme) << "No scheme provided!" << std::endl;

        const auto assembling_timer = BuiltinTimer();
        const SizeType number_of_rom_modes = this->GetNumberOfROMModes();
        const auto& r_current_process_info = rModelPart.GetProcessInfo();

        rA = ZeroMatrix(number_of_rom_modes, number_of_rom_modes);
        rb = ZeroVector(number_of_rom_modes);

        // First loop: test functions
        BuildJPhi(*pScheme, rModelPart);

        // Second loop: projection of the (weighted) entities onto the test functions
        using SystemSumReducer = CombinedReduction<typename BaseType::template NonTrivialSumReduction<RomSystemMatrixType>, typename BaseType::template NonTrivialSumReduction<RomSystemVectorType>>;
        LSPGAssemblyTLS assembly_tls_container(number_of_rom_modes);

        auto& r_elements = this->mHromSimulation ? this->mSelectedElements : rModelPart.Elements();
        if (!r_elements.empty()) {
            std::tie(rA, rb) = block_for_each<SystemSumReducer>(r_elements, assembly_tls_container,
                [&](Element& rElement, LSPGAssemblyTLS& rTLS)
            {
                return CalculateLocalContribution(rElement, rTLS, *pScheme, r_current_process_info);
            });
        }

        auto& r_conditions = this->mHromSimulation ? this->mSelectedConditions : rModelPart.Conditions();
        if (!r_conditions.empty()) {
            RomSystemMatrixType a_conditions;
            RomSystemVectorType b_conditions;
            std::tie(a_conditions, b_conditions) = block_for_each<SystemSumReducer>(r_conditions, assembly_tls_container,
                [&](Condition& rCondition, LSPGAssemblyTLS& rTLS)
            {
                return CalculateLocalContribution(rCondition, rTLS, *pScheme, r_current_process_info);
            });
            rA += a_conditions;
            rb += b_conditions;
        }

        KRATOS_INFO_IF("ElementalLeastSquaresPetrovGalerkinROMBuilderAndSolver", (this->GetEchoLevel() > 0)) << "Build time: " << assembling_timer.ElapsedSeconds() << std::endl;

        KRATOS_CATCH("")
    }

    ///@}
private:

    ///@name Private member variables
    ///@{

    Matrix mJPhi;                               /// Test functions: the Jacobian times the right basis
    std::vector<LockObject> mRowLocks;          /// One lock per row of the test functions
    ElementsArrayType mComplementaryElements;   /// Selected elements and their nodal neighbours (HROM)
    ConditionsArrayType mComplementaryConditions; /// Selected conditions and their nodal neighbours (HROM)

    ///@}
    ///@name Private operations
    ///@{

    /**
     * Assembles the test functions (J*Phi) from the elemental contributions, without HROM weights
     */
    void BuildJPhi(
        TSchemeType& rScheme,
        ModelPart& rModelPart)
    {
        const SizeType system_size = BaseBuilderAndSolverType::GetEquationSystemSize();
        const SizeType number_of_rom_modes = this->GetNumberOfROMModes();
        const auto& r_current_process_info = rModelPart.GetProcessInfo();

        if (mJPhi.size1() != system_size || mJPhi.size2() != number_of_rom_modes) {
            mJPhi.resize(system_size, number_of_rom_modes, false);
        }
        noalias(mJPhi) = ZeroMatrix(system_size, number_of_rom_modes);
        if (mRowLocks.size() != system_size) {
            mRowLocks = std::vector<LockObject>(system_size);
        }

        LSPGAssemblyTLS assembly_tls_container(number_of_rom_modes);

        auto& r_elements = this->mHromSimulation ? mComplementaryElements : rModelPart.Elements();
        block_for_each(r_elements, assembly_tls_container, [&](Element& rElement, LSPGAssemblyTLS& rTLS)
        {
            AddJPhiContribution(rElement, rTLS, rScheme, r_current_process_info);
        });

        auto& r_conditions = this->mHromSimulation ? mComplementaryConditions : rModelPart.Conditions();
        block_for_each(r_conditions, assembly_tls_container, [&](Condition& rCondition, LSPGAssemblyTLS& rTLS)
        {
            AddJPhiContribution(rCondition, rTLS, rScheme, r_current_process_info);
        });
    }

    /**
     * Adds the contribution of an element or condition to the test functions
     */
    template<typename TEntity>
    void AddJPhiContribution(
        TEntity& rEntity,
        LSPGAssemblyTLS& rTLS,
        TSchemeType& rScheme,
        const ProcessInfo& rCurrentProcessInfo)
    {
        if (rEntity.IsDefined(ACTIVE) && rEntity.IsNot(ACTIVE)) {
            return;
        }

        rScheme.CalculateSystemContributions(rEntity, rTLS.lhs, rTLS.rhs, rTLS.eq_id, rCurrentProcessInfo);
        // Entities without LHS (e.g. loads) do not contribute to the test functions
        if (rTLS.lhs.size1() == 0) {
            return;
        }
        rEntity.GetDofList(rTLS.dofs, rCurrentProcessInfo);

        const SizeType ndofs = rTLS.dofs.size();
        BaseType::ResizeIfNeeded(rTLS.phiE, ndofs, this->GetNumberOfROMModes());
        BaseType::ResizeIfNeeded(rTLS.aux, ndofs, this->GetNumberOfROMModes());

        RomAuxiliaryUtilities::GetPhiElemental(rTLS.phiE, rTLS.dofs, rEntity.GetGeometry(), this->mMapPhi);
        noalias(rTLS.aux) = prod(rTLS.lhs, rTLS.phiE);

        for (IndexType i = 0; i < ndofs; ++i) {
            const auto& r_dof = *rTLS.dofs[i];
            if (!r_dof.IsFixed()) {
                const IndexType equation_id = r_dof.EquationId();
                const std::lock_guard<LockObject> scope_lock(mRowLocks[equation_id]);
                noalias(row(mJPhi, equation_id)) += row(rTLS.aux, i);
            }
        }
    }

    /**
     * Computes the contribution of an element or condition to the reduced system
     */
    template<typename TEntity>
    std::tuple<RomSystemMatrixType, RomSystemVectorType> CalculateLocalContribution(
        TEntity& rEntity,
        LSPGAssemblyTLS& rTLS,
        TSchemeType& rScheme,
        const ProcessInfo& rCurrentProcessInfo)
    {
        const SizeType number_of_rom_modes = this->GetNumberOfROMModes();
        if (rEntity.IsDefined(ACTIVE) && rEntity.IsNot(ACTIVE)) {
            rTLS.romA = ZeroMatrix(number_of_rom_modes, number_of_rom_modes);
            rTLS.romB = ZeroVector(number_of_rom_modes);
            return std::tie(rTLS.romA, rTLS.romB);
        }

        rScheme.CalculateSystemContributions(rEntity, rTLS.lhs, rTLS.rhs, rTLS.eq_id, rCurrentProcessInfo);
        rEntity.GetDofList(rTLS.dofs, rCurrentProcessInfo);

        const SizeType ndofs = rTLS.dofs.size();
        BaseType::ResizeIfNeeded(rTLS.psiE, ndofs, number_of_rom_modes);
        RomAuxiliaryUtilities::GetJPhiElemental(rTLS.psiE, rTLS.dofs, mJPhi);

        const double h_rom_weight = this->mHromSimulation ? rEntity.GetValue(HROM_WEIGHT)[this->mActiveHromSet] : 1.0;

        if (rTLS.lhs.size1() == 0) {
            rTLS.romA = ZeroMatrix(number_of_rom_modes, number_of_rom_modes);
        } else {
            BaseType::ResizeIfNeeded(rTLS.phiE, ndofs, number_of_rom_modes);
            BaseType::ResizeIfNeeded(rTLS.aux, ndofs, number_of_rom_modes);
            RomAuxiliaryUtilities::GetPhiElemental(rTLS.phiE, rTLS.dofs, rEntity.GetGeometry(), this->mMapPhi);
            noalias(rTLS.aux) = prod(rTLS.lhs, rTLS.phiE);
            noalias(rTLS.romA) = prod(trans(rTLS.psiE), rTLS.aux) * h_rom_weight;
        }
        noalias(rTLS.romB) = prod(trans(rTLS.psiE), rTLS.rhs) * h_rom_weight;

        return std::tie(rTLS.romA, rTLS.romB);
    }

    ///@}
}; /* Class ElementalLeastSquaresPetrovGalerkinROMBuilderAndSolver */

///@}

} /* namespace Kratos */
