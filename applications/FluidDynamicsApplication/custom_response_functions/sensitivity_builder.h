//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    
//

#if !defined(KRATOS_SENSITIVITY_BUILDER_H_INCLUDED)
#define KRATOS_SENSITIVITY_BUILDER_H_INCLUDED

// System includes
#include <vector>
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "utilities/variable_utils.h"
#include "response_functions/adjoint_response_function.h"
#include "solving_strategies/response_functions/nodal_sensitivity_builder.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Classes
///@{

class SensitivityBuilder
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(SensitivityBuilder);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    SensitivityBuilder(Parameters& rParameters, ModelPart& rModelPart, AdjointResponseFunction::Pointer pResponseFunction)
    : mrModelPart(rModelPart), mpResponseFunction(pResponseFunction)
    {
        KRATOS_TRY;

        Parameters default_settings(R"(
        {
            "sensitivity_model_part_name": "PLEASE_SPECIFY_SENSITIVITY_MODEL_PART",
            "nodal_sensitivity_variables": ["SHAPE_SENSITIVITY"],
            "integrate_in_time": true
        })");

        Parameters custom_settings = rParameters;
        custom_settings.ValidateAndAssignDefaults(default_settings);

        mSensitivityModelPartName =
            custom_settings["sensitivity_model_part_name"].GetString();

        Parameters nodal_sensitivity_variables = custom_settings["nodal_sensitivity_variables"];
        mNodalSensitivityVariables.resize(nodal_sensitivity_variables.size());
        for (unsigned int i = 0; i < nodal_sensitivity_variables.size(); ++i)
            mNodalSensitivityVariables[i] = nodal_sensitivity_variables[i].GetString();

        mIntegrateInTime = custom_settings["integrate_in_time"].GetBool();

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void Initialize()
    {
        KRATOS_TRY;

        Check();

        for (const std::string& r_label : mNodalSensitivityVariables)
            SetNodalSensitivityVariableToZero(r_label);

        VariableUtils().SetNonHistoricalVariable(UPDATE_SENSITIVITIES, false,
                                                 mrModelPart.Nodes());
        VariableUtils().SetNonHistoricalVariable(
            UPDATE_SENSITIVITIES, true,
            mrModelPart.GetSubModelPart(mSensitivityModelPartName).Nodes());

        KRATOS_CATCH("");
    }

    void UpdateSensitivities()
    {
        KRATOS_TRY;

        double delta_time;
        if (mIntegrateInTime)
            delta_time = -mrModelPart.GetProcessInfo()[DELTA_TIME];
        else
            delta_time = 1.0;

        NodalSensitivityBuilder nodal_sensitivity_builder(mrModelPart, mpResponseFunction);
        for (const std::string& r_label : mNodalSensitivityVariables)
            nodal_sensitivity_builder.BuildNodalSolutionStepSensitivities(r_label, delta_time);

        KRATOS_CATCH("");
    }

    ///@}

protected:
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}

private:
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;
    AdjointResponseFunction::Pointer mpResponseFunction;
    std::string mSensitivityModelPartName;
    std::vector<std::string> mNodalSensitivityVariables;
    bool mIntegrateInTime;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void Check()
    {
        KRATOS_TRY;

        if (mrModelPart.HasSubModelPart(mSensitivityModelPartName) == false)
            KRATOS_ERROR << "No sub model part \"" << mSensitivityModelPartName
                         << "\"" << std::endl;

        KRATOS_CATCH("");
    }

    void SetNodalSensitivityVariableToZero(std::string const& rVariableName)
    {
        KRATOS_TRY;

        if (KratosComponents<Variable<double>>::Has(rVariableName) == true)
        {
            const Variable<double>& r_variable =
                KratosComponents<Variable<double>>::Get(rVariableName);

            VariableUtils().SetScalarVar(r_variable, r_variable.Zero(), mrModelPart.Nodes());
        }
        else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(rVariableName) == true)
        {
            const Variable<array_1d<double, 3>>& r_variable =
                KratosComponents<Variable<array_1d<double, 3>>>::Get(rVariableName);

            VariableUtils().SetVectorVar(r_variable, r_variable.Zero(), mrModelPart.Nodes());
        }
        else
            KRATOS_ERROR << "Unsupported variable: " << rVariableName << "." << std::endl;

        KRATOS_CATCH("");
    }


    ///@}
};

///@} // Kratos Classes

///@} // FluidDynamicsApplication group

} /* namespace Kratos.*/

#endif /* KRATOS_SENSITIVITY_BUILDER_H_INCLUDED defined */
