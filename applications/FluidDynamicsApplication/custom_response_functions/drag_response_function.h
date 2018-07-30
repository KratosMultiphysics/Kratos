//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Michael Andre, https://github.com/msandre
//

#if !defined(KRATOS_DRAG_RESPONSE_FUNCTION_H_INCLUDED)
#define KRATOS_DRAG_RESPONSE_FUNCTION_H_INCLUDED

// System includes
#include <vector>
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "includes/kratos_parameters.h"
#include "utilities/variable_utils.h"
#include "solving_strategies/response_functions/response_function.h"
#include "solving_strategies/response_functions/response_function_sensitivity_builder_utility.h"

namespace Kratos
{
///@addtogroup AdjointFluidApplication
///@{

///@name Kratos Classes
///@{

class SensitivityBuilder : protected ResponseFunctionSensitivityBuilderUtility
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(SensitivityBuilder);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    SensitivityBuilder(Parameters& rParameters, ModelPart& rModelPart, ResponseFunction::Pointer pResponseFunction)
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
            SetNodalSensitivityVariableToZero(r_label, mrModelPart.Nodes());

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
        for (const std::string& r_label : mNodalSensitivityVariables)
            BuildNodalSolutionStepSensitivities(r_label, mrModelPart, delta_time);

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

    void CalculatePartialSensitivity(Variable<array_1d<double, 3>> const& rVariable,
                                     Element const& rElement,
                                     Matrix const& rSensitivityMatrix,
                                     Vector& rPartialSensitivity,
                                     ProcessInfo const& rProcessInfo) const override
    {
        mpResponseFunction->CalculatePartialSensitivity2(
            rVariable, rElement, rSensitivityMatrix, rPartialSensitivity, rProcessInfo);
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;
    ResponseFunction::Pointer mpResponseFunction;
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

    void SetNodalSensitivityVariableToZero(std::string const& rVariableName,
                                           NodesContainerType& rNodes)
    {
        KRATOS_TRY;

        if (KratosComponents<Variable<double>>::Has(rVariableName) == true)
        {
            const Variable<double>& r_variable =
                KratosComponents<Variable<double>>::Get(rVariableName);

            VariableUtils().SetScalarVar(r_variable, r_variable.Zero(), rNodes);
        }
        else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(rVariableName) == true)
        {
            const Variable<array_1d<double, 3>>& r_variable =
                KratosComponents<Variable<array_1d<double, 3>>>::Get(rVariableName);

            VariableUtils().SetVectorVar(r_variable, r_variable.Zero(), rNodes);
        }
        else
            KRATOS_ERROR << "Unsupported variable: " << rVariableName << "." << std::endl;

        KRATOS_CATCH("");
    }


    ///@}
};

/// A response function for drag.
/**
 * The response function is defined as:
 *
 * \f[
 * \bar{D} = \Sigma_{n=1}^N D^n \Delta t
 * \f]
 *
 * if "integrate_in_time" is true.
 */
template <unsigned int TDim>
class DragResponseFunction : public ResponseFunction
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(DragResponseFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    DragResponseFunction(Parameters& rParameters, ModelPart& rModelPart)
    : mrModelPart(rModelPart)
    {
        KRATOS_TRY;

        Parameters default_settings(R"(
        {
            "structure_model_part_name": "PLEASE_SPECIFY_STRUCTURE_MODEL_PART",
            "drag_direction": [1.0, 0.0, 0.0]
        })");

        Parameters custom_settings = rParameters["custom_settings"];
        custom_settings.ValidateAndAssignDefaults(default_settings);

        mStructureModelPartName =
            custom_settings["structure_model_part_name"].GetString();

        if (custom_settings["drag_direction"].IsArray() == false ||
            custom_settings["drag_direction"].size() != 3)
        {
            KRATOS_ERROR << "Invalid \"drag_direction\"." << std::endl;
        }

        for (unsigned int d = 0; d < TDim; ++d)
            mDragDirection[d] = custom_settings["drag_direction"][d].GetDouble();

        if (std::abs(norm_2(mDragDirection) - 1.0) > 1e-3)
        {
            const double magnitude = norm_2(mDragDirection);
            if (magnitude == 0.0)
                KRATOS_ERROR << "\"drag_direction\" is zero." << std::endl;

            std::cout << "WARNING: Non unit magnitude in \"drag_direction\"." << std::endl;
            std::cout << "WARNING: Normalizing ..." << std::endl;

            for (unsigned int d = 0; d < TDim; d++)
                mDragDirection[d] /= magnitude;
        }

        KRATOS_CATCH("");
    }

    /// Destructor.
    ~DragResponseFunction() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void Initialize() override
    {
        KRATOS_TRY;

        Check();

        VariableUtils().SetFlag(STRUCTURE, false, mrModelPart.Nodes());
        VariableUtils().SetFlag(
            STRUCTURE, true, mrModelPart.GetSubModelPart(mStructureModelPartName).Nodes());

        KRATOS_CATCH("");
    }

    void CalculateGradient(Element const& rElement,
                           Matrix const& rAdjointMatrix,
                           Vector& rResponseGradient,
                           ProcessInfo const& rProcessInfo) const override
    {
        CalculateDragContribution(
            rAdjointMatrix, rElement.GetGeometry().Points(), rResponseGradient);
    }

    void CalculateFirstDerivativesGradient(Element const& rElement,
                                           Matrix const& rAdjointMatrix,
                                           Vector& rResponseGradient,
                                           ProcessInfo const& rProcessInfo) const override
    {
        CalculateDragContribution(
            rAdjointMatrix, rElement.GetGeometry().Points(), rResponseGradient);
    }

    void CalculateSecondDerivativesGradient(Element const& rElement,
                                            Matrix const& rAdjointMatrix,
                                            Vector& rResponseGradient,
                                            ProcessInfo const& rProcessInfo) const override
    {
        CalculateDragContribution(
            rAdjointMatrix, rElement.GetGeometry().Points(), rResponseGradient);
    }

    // This is a temporary crutch to upgrade the response function base class
    // without completely breaking the tests.
    virtual void CalculatePartialSensitivity2(Variable<array_1d<double, 3>> const& rVariable,
                                             Element const& rElement,
                                             Matrix const& rSensitivityMatrix,
                                             Vector& rPartialSensitivity,
                                             ProcessInfo const& rProcessInfo) const override
    {
        KRATOS_TRY;

        CalculateDragContribution(
            rSensitivityMatrix, rElement.GetGeometry().Points(), rPartialSensitivity);

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
    std::string mStructureModelPartName;
    array_1d<double, TDim> mDragDirection;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void Check()
    {
        KRATOS_TRY;

        if (mrModelPart.HasSubModelPart(mStructureModelPartName) == false)
            KRATOS_ERROR << "No sub model part \"" << mStructureModelPartName
                         << "\"" << std::endl;

        KRATOS_CATCH("");
    }

    void CalculateDragContribution(const Matrix& rDerivativesOfResidual,
                                   const Element::NodesArrayType& rNodes,
                                   Vector& rDerivativesOfDrag) const
    {
        constexpr std::size_t max_size = 50;
        BoundedVector<double, max_size> drag_flag_vector(rDerivativesOfResidual.size2());

        const unsigned num_nodes = rNodes.size();
        unsigned local_index = 0;
        for (unsigned i_node = 0; i_node < num_nodes; ++i_node)
        {
            if (rNodes[i_node].Is(STRUCTURE))
            {
                for (unsigned d = 0; d < TDim; ++d)
                    drag_flag_vector[local_index++] = mDragDirection[d];
            }
            else
            {
                for (unsigned int d = 0; d < TDim; ++d)
                    drag_flag_vector[local_index++] = 0.0;
            }

            drag_flag_vector[local_index++] = 0.0; // pressure dof
        }

        if (rDerivativesOfDrag.size() != rDerivativesOfResidual.size1())
            rDerivativesOfDrag.resize(rDerivativesOfResidual.size1(), false);

        noalias(rDerivativesOfDrag) = prod(rDerivativesOfResidual, drag_flag_vector);
    }

    ///@}
};

///@} // Kratos Classes

///@} // Adjoint Fluid Application group

} /* namespace Kratos.*/

#endif /* KRATOS_DRAG_RESPONSE_FUNCTION_H_INCLUDED defined */
