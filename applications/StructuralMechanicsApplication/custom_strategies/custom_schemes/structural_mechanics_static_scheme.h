//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/variables.h"
#include "processes/calculate_nodal_area_process.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"

// Application includes

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
 * @class StructuralMechanicsStaticScheme
 * @ingroup KratosCore
 * @brief This class provides the implementation of a static scheme
 * @details The only operation done in this  scheme is the update of the database, no predict is done
 * @tparam TSparseSpace The sparse space considered
 * @tparam TDenseSpace The dense space considered
 * @see Scheme
 * @author Riccardo Rossi
 */
template<class TSparseSpace, class TDenseSpace>
class StructuralMechanicsStaticScheme
    : public ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace>
{

public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of StructuralMechanicsStaticScheme
    KRATOS_CLASS_POINTER_DEFINITION( StructuralMechanicsStaticScheme);

    /// Base class definition
    using BaseType = Scheme<TSparseSpace, TDenseSpace>;

    /// Base class definition
    using StaticBaseType = ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace, TDenseSpace>;

    // The current class definition
    using ClassType = StructuralMechanicsStaticScheme<TSparseSpace, TDenseSpace>;

    /// Matrix type definition
    using TSystemMatrixType = typename StaticBaseType::TSystemMatrixType;

    /// Vector type definition
    using TSystemVectorType = typename StaticBaseType::TSystemVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    explicit StructuralMechanicsStaticScheme(Parameters ThisParameters)
        : StaticBaseType()
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);
    }

    explicit StructuralMechanicsStaticScheme()
        : StaticBaseType()
    {}

    explicit StructuralMechanicsStaticScheme(StructuralMechanicsStaticScheme& rOther)
        :StaticBaseType(rOther)
    {}

    ~StructuralMechanicsStaticScheme() override
    {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    typename BaseType::Pointer Create(Parameters ThisParameters) const override
    {
        return Kratos::make_shared<ClassType>(ThisParameters);
    }

    void Initialize(ModelPart& rModelPart) override
    {
        // Call the base Initialize method
        StaticBaseType::Initialize(rModelPart);

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

        // Calculate the NODAL_AREA
        if (oss_switch && mCalculateNodalArea) {
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
        StaticBaseType::FinalizeNonLinIteration(rModelPart, rA, rDx, rb);

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
            "name" : "structural_mechanics_static_scheme",
            "calculate_nodal_area" : true,
            "projection_variables_list" : []
        })");

        // Getting base class default parameters
        const Parameters base_default_parameters = StaticBaseType::GetDefaultParameters();
        default_parameters.RecursivelyAddMissingParameters(base_default_parameters);
        return default_parameters;
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "structural_mechanics_static_scheme";
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "StructuralMechanicsStaticScheme";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}
protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    void AssignSettings(const Parameters ThisParameters) override
    {
        StaticBaseType::AssignSettings(ThisParameters);

        mCalculateNodalArea = ThisParameters["calculate_nodal_area"].GetBool();
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
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}
private:
    ///@name Private static Member Variables
    ///@{

    ///@}
    ///@name Private member Variables
    ///@{

    bool mCalculateNodalArea;

    std::vector<const Variable<double>*> mProjectionScalarVariablesList;

    std::vector<const Variable<array_1d<double,3>>*> mProjectionVectorVariablesList;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Private LifeCycle
    ///@{

    ///@}
}; // Class StructuralMechanicsStaticScheme

}  // namespace Kratos
