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
#include "utilities/nodal_sensitivity_builder.h"
#include "utilities/element_data_sensitivity_builder.h"

namespace Kratos
{
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
    SensitivityBuilder(Parameters Settings,
                       ModelPart& rModelPart,
                       AdjointResponseFunction::Pointer pResponseFunction)
        : mrModelPart(rModelPart), mpResponseFunction(pResponseFunction)
    {
        KRATOS_TRY;

        Parameters default_settings(R"(
        {
            "sensitivity_model_part_name": "PLEASE_SPECIFY_SENSITIVITY_MODEL_PART",
            "nodal_solution_step_sensitivity_variables": [],
            "element_data_sensitivity_variables" : [],
            "build_mode": "integrate"
        })");

        Settings.ValidateAndAssignDefaults(default_settings);

        auto sensitivity_model_part_name =
            Settings["sensitivity_model_part_name"].GetString();
        if (sensitivity_model_part_name !=
            "PLEASE_SPECIFY_SENSITIVITY_MODEL_PART")
        {
            mpSensitivityModelPart =
                &mrModelPart.GetSubModelPart(sensitivity_model_part_name);
        }
        else
        {
            mpSensitivityModelPart = &mrModelPart;
        }

        Parameters nodal_sensitivity_variables =
            Settings["nodal_solution_step_sensitivity_variables"];
        mNodalSolutionStepSensitivityVariables.resize(nodal_sensitivity_variables.size());
        for (unsigned int i = 0; i < nodal_sensitivity_variables.size(); ++i)
        {
            mNodalSolutionStepSensitivityVariables[i] = nodal_sensitivity_variables[i].GetString();
        }

        Parameters element_sensitivity_variables =
            Settings["element_data_sensitivity_variables"];
        mElementDataSensitivityVariables.resize(element_sensitivity_variables.size());
        for (unsigned int i = 0; i < element_sensitivity_variables.size(); ++i)
        {
            mElementDataSensitivityVariables[i] =
                element_sensitivity_variables[i].GetString();
        }

        mBuildMode = Settings["build_mode"].GetString();

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

        SetAllSensitivityVariablesToZero();

        VariableUtils().SetNonHistoricalVariable(UPDATE_SENSITIVITIES, false,
                                                 mrModelPart.Nodes());
        VariableUtils().SetNonHistoricalVariable(
            UPDATE_SENSITIVITIES, true, mpSensitivityModelPart->Nodes());

        KRATOS_CATCH("");
    }

    void UpdateSensitivities()
    {
        KRATOS_TRY;

        double weight;
        if (mBuildMode == "integrate")
        {
            // integrate in time.
            weight = -mrModelPart.GetProcessInfo()[DELTA_TIME];
        }
        else if (mBuildMode == "sum")
        {
            // sum in time.
            weight = 1.0;
        }
        else if (mBuildMode == "static")
        {
            // overwrite existing.
            SetAllSensitivityVariablesToZero();
            weight = 1.0;
        }
        else
        {
            KRATOS_ERROR << "Unsupported \"build_mode\": " << mBuildMode << std::endl;
        }

        NodalSensitivityBuilder nodal_sensitivity_builder(mrModelPart, mpResponseFunction);
        for (const std::string& r_label : mNodalSolutionStepSensitivityVariables)
        {
            nodal_sensitivity_builder.BuildNodalSolutionStepSensitivities(r_label, weight);
        }
        
        ElementSensitivityBuilder element_sensitivity_builder(mrModelPart, mpResponseFunction);
        for (const std::string& r_label : mElementDataSensitivityVariables)
        {
            element_sensitivity_builder.BuildElementSensitivities(r_label, weight);
        }

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
    ModelPart* mpSensitivityModelPart = nullptr;
    AdjointResponseFunction::Pointer mpResponseFunction;
    std::vector<std::string> mNodalSolutionStepSensitivityVariables;
    std::vector<std::string> mElementDataSensitivityVariables;
    std::string mBuildMode;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void SetAllSensitivityVariablesToZero()
    {
        KRATOS_TRY;

        for (const std::string& r_label : mNodalSolutionStepSensitivityVariables)
        {
            SetNodalSolutionStepSensitivityVariableToZero(r_label);
        }

        for (const std::string& r_label : mElementDataSensitivityVariables)
        {
            SetElementDataSensitivityVariableToZero(r_label);
        }

        KRATOS_CATCH("");
    }

    void SetNodalSolutionStepSensitivityVariableToZero(std::string const& rVariableName)
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

    void SetElementDataSensitivityVariableToZero(std::string const& rVariableName)
    {
        KRATOS_TRY;

        auto& r_elements = mrModelPart.Elements();
        if (KratosComponents<Variable<double>>::Has(rVariableName) == true)
        {
            const Variable<double>& r_variable =
                KratosComponents<Variable<double>>::Get(rVariableName);
            #pragma omp parallel for
            for (int k = 0; k< static_cast<int> (r_elements.size()); ++k)
            {
                auto it = r_elements.begin() + k;
                it->SetValue(r_variable, r_variable.Zero());
            }
        }
        else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(rVariableName) == true)
        {
            const Variable<array_1d<double, 3>>& r_variable =
                KratosComponents<Variable<array_1d<double, 3>>>::Get(rVariableName);
            #pragma omp parallel for
            for (int k = 0; k< static_cast<int> (r_elements.size()); ++k)
            {
                auto it = r_elements.begin() + k;
                it->SetValue(r_variable, r_variable.Zero());
            }
        }
        else
        {
            KRATOS_ERROR << "Unsupported variable: " << rVariableName << "." << std::endl;
        }

        KRATOS_CATCH("");
    }


    ///@}
};

///@} // Kratos Classes

} /* namespace Kratos.*/

#endif /* KRATOS_SENSITIVITY_BUILDER_H_INCLUDED defined */
