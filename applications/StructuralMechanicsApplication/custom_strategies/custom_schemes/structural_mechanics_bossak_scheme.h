//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/cfd_variables.h" //TODO: For OSS_SWITCH (think about a better place for this variable)
#include "processes/calculate_nodal_area_process.h"

// Application includes
#include "solving_strategies/schemes/residual_based_bossak_displacement_scheme.hpp"

namespace Kratos
{

///@name Kratos Classes
///@{

/**
 * @brief Structural mechanics Bossak scheme
 * This scheme extends the implementation in the ResidualBasedBossakDisplacementScheme to do some structural mechanics specifics
 * @tparam TSparseSpace The sparse space considered
 * @tparam TDenseSpace The dense space considered
 */
template<class TSparseSpace,  class TDenseSpace >
class StructuralMechanicsBossakScheme
    : public ResidualBasedBossakDisplacementScheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(StructuralMechanicsBossakScheme);

    /// Base type for the scheme
    using BaseType = Scheme<TSparseSpace, TDenseSpace>;

    /// Implicit base type for the scheme
    using BossakBaseType = ResidualBasedBossakDisplacementScheme<TSparseSpace, TDenseSpace>;

    /// Class type for the scheme
    using ClassType = StructuralMechanicsBossakScheme<TSparseSpace, TDenseSpace>;

    /// Type for the system matrix within BossakBaseType
    using TSystemMatrixType = typename BossakBaseType::TSystemMatrixType;

    /// Type for the system vector within BossakBaseType
    using TSystemVectorType = typename BossakBaseType::TSystemVectorType;

    /// Pointer type for the BaseType
    using BaseTypePointer = typename BaseType::Pointer;

    ///@}
    ///@name Life Cycle
    ///@{

    explicit StructuralMechanicsBossakScheme(Parameters ThisParameters)
        : BossakBaseType()
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);

        // For pure Newmark Scheme
        this->mNewmark.gamma = 0.5;

        this->AuxiliarInitializeBossak();
    }

    explicit StructuralMechanicsBossakScheme(StructuralMechanicsBossakScheme& rOther)
        :BossakBaseType(rOther)
    {}

    ~StructuralMechanicsBossakScheme() override
    {}

    ///@}
    ///@name Operations
    ///@{

    typename BaseType::Pointer Create(Parameters ThisParameters) const override
    {
        return Kratos::make_shared<ClassType>(ThisParameters);
    }

    BaseTypePointer Clone() override
    {
        return BaseTypePointer(new StructuralMechanicsBossakScheme(*this));
    }

    void Initialize(ModelPart& rModelPart) override
    {
        // Call the base Initialize method
        BossakBaseType::Initialize(rModelPart);

        // Allocate the OSS projection variables
        const auto &r_process_info = rModelPart.GetProcessInfo();
        const bool oss_switch = r_process_info.Has(OSS_SWITCH) ? r_process_info[OSS_SWITCH] : false;
        if (oss_switch) {
            const array_1d<double,3> aux_zero = ZeroVector(3);
            block_for_each(rModelPart.Nodes(), aux_zero, [this](Node& rNode, array_1d<double,3>& rAuxZero){
                for (auto& rp_var : mProjectionScalarVariablesList) {
                    rNode.SetValue(*rp_var, 0.0);
                }
                for (auto& rp_var : mProjectionVectorVariablesList) {
                    rNode.SetValue(*rp_var, rAuxZero);
                }
            });
        }
    }

    void InitializeSolutionStep(
        ModelPart &rModelPart,
        TSystemMatrixType &A,
        TSystemVectorType &Dx,
        TSystemVectorType &b) override
    {
        // Call the base InitializeSolutionStep method
        BossakBaseType::InitializeSolutionStep(rModelPart, A, Dx, b);

        // Update the NODAL_AREA (note that this strictly required only in the updated Lagrangian case)
        const auto& r_current_process_info = rModelPart.GetProcessInfo();
        const bool oss_switch = r_current_process_info.Has(OSS_SWITCH) ? r_current_process_info[OSS_SWITCH] : false;
        if (oss_switch && mUpdateNodalArea) {
            CalculateNodalAreaProcess<false>(rModelPart).Execute();
        }
    }

    void FinalizeNonLinIteration(
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb) override
    {
        KRATOS_TRY

        // Call base class FinalizeNonLinIteration
        BossakBaseType::FinalizeNonLinIteration(rModelPart, rA, rDx, rb);

        // Check if the Orthogonal SubScales (OSS) are active
        const auto& r_process_info = rModelPart.GetProcessInfo();
        const bool oss_switch = r_process_info.Has(OSS_SWITCH) ? r_process_info[OSS_SWITCH] : false;

        // Calculate the OSS projections
        if (oss_switch) {
            // Initialize the projection values
            const array_1d<double,3> aux_zero = ZeroVector(3);
            block_for_each(rModelPart.Nodes(), aux_zero, [this](Node& rNode, array_1d<double,3>& rAuxZero){
                for (auto& rp_var : mProjectionScalarVariablesList) {
                    rNode.FastGetSolutionStepValue(*rp_var) = 0.0;
                }
                for (auto& rp_var : mProjectionVectorVariablesList) {
                    rNode.FastGetSolutionStepValue(*rp_var) = rAuxZero;
                }
            });

            // Calculate the element residuals projection
            std::tuple<double, array_1d<double,3>> oss_proj_tls;
            block_for_each(rModelPart.Elements(), oss_proj_tls, [&, this](Element& rElement, std::tuple<double, array_1d<double,3>>& rOssProjTLS){
                double& r_scalar_proj = std::get<0>(rOssProjTLS);
                array_1d<double,3>& r_vector_proj = std::get<1>(rOssProjTLS);
                for (auto& rp_var : mProjectionScalarVariablesList) {
                    rElement.Calculate(*rp_var, r_scalar_proj, r_process_info);
                }
                for (auto& rp_var : mProjectionVectorVariablesList) {
                    rElement.Calculate(*rp_var, r_vector_proj, r_process_info);
                }
            });

            // Do the nodal weighting
            //TODO: We the weighted L2 projection with the density should be done in the multimaterial case
            //TODO: Note that this would require customizing the weights calculation in here (standard NODAL_AREA cannot be used)
            block_for_each(rModelPart.Nodes(), [this](Node& rNode){
                const double nodal_area = rNode.GetValue(NODAL_AREA);
                for (auto& rp_var : mProjectionScalarVariablesList) {
                    rNode.FastGetSolutionStepValue(*rp_var) /= nodal_area;
                }
                for (auto& rp_var : mProjectionVectorVariablesList) {
                    rNode.FastGetSolutionStepValue(*rp_var) /= nodal_area;
                }
            });
        }

        KRATOS_CATCH("")
    }

    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"({
            "name" : "structural_mechanics_bossak_scheme",
            "update_nodal_area" : true,
            "projection_variables_list" : []
        })");

        // Getting base class default parameters
        const Parameters base_default_parameters = BossakBaseType::GetDefaultParameters();
        default_parameters.RecursivelyAddMissingParameters(base_default_parameters);
        return default_parameters;
    }

    /// @brief Return the name of the class as used in the settings (@a snake_case).
    static std::string Name()
    {
        return "structural_mechanics_bossak_scheme";
    }

    ///@}
    ///@name Input and output
    ///@{

    /// @brief Return information as a string.
    std::string Info() const override
    {
        return "StructuralMechanicsBossakScheme";
    }

    /// @brief Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// @brief Print the instance's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    ///@}
protected:
    ///@name Protected Operations
    ///@{

    void AssignSettings(const Parameters ThisParameters) override
    {
        BossakBaseType::AssignSettings(ThisParameters);

        mUpdateNodalArea = ThisParameters["update_nodal_area"].GetBool();
        const auto& r_proj_var_list = ThisParameters["projection_variables_list"].GetStringArray();
        for (const std::string& r_var_name : r_proj_var_list) {
            if (KratosComponents<Variable<double>>::Has(r_var_name)) {
                mProjectionScalarVariablesList.push_back(&KratosComponents<Variable<double>>::Get(r_var_name));
            } else if (KratosComponents<Variable<array_1d<double,3>>>::Has(r_var_name)) {
                mProjectionVectorVariablesList.push_back(&KratosComponents<Variable<array_1d<double, 3>>>::Get(r_var_name));
            } else {
                KRATOS_ERROR << "Wrong projection variable '" << r_var_name << "'." << std::endl;
            }
        }
    }

    ///@}
private:
    ///@name Private Member Variables
    ///@{

    bool mUpdateNodalArea;

    std::vector<const Variable<double>*> mProjectionScalarVariablesList;

    std::vector<const Variable<array_1d<double,3>>*> mProjectionVectorVariablesList;

    ///@}
}; // Class StructuralMechanicsBossakScheme

///@}

} // namespace Kratos
