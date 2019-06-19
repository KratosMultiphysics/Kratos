// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Jordi Cotela
//

#ifndef KRATOS_POINT_TEMPERATURE_RESPONSE_FUNCTION_H_INCLUDED
#define KRATOS_POINT_TEMPERATURE_RESPONSE_FUNCTION_H_INCLUDED

#include "includes/kratos_flags.h"
#include "includes/model_part.h"
#include "response_functions/adjoint_response_function.h"

namespace Kratos {
///@addtogroup ConvectionDiffusionApplication
///@{

///@name Kratos Classes
///@{

class PointTemperatureResponseFunction: public AdjointResponseFunction
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(PointTemperatureResponseFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    PointTemperatureResponseFunction(Parameters Settings, ModelPart& rModelPart)
    : mrModelPart(rModelPart)
    {
        KRATOS_TRY;

        KRATOS_CATCH("");
    }

    /// Destructor.
    ~PointTemperatureResponseFunction() override
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

        KRATOS_CATCH("");
    }

    void CalculateGradient(const Element& rAdjointElement,
                           const Matrix& rResidualGradient,
                           Vector& rResponseGradient,
                           const ProcessInfo& rProcessInfo) override
    {
        ComputePointTemperatureSensitivityContribution(rResidualGradient, rAdjointElement.GetGeometry().Points(),rResponseGradient);
    }

    void CalculateGradient(const Condition& rAdjointCondition,
                           const Matrix& rResidualGradient,
                           Vector& rResponseGradient,
                           const ProcessInfo& rProcessInfo) override
    {
        noalias(rResponseGradient) = ZeroVector(rResidualGradient.size2());
    }

    void CalculateFirstDerivativesGradient(const Element& rAdjointElement,
                                           const Matrix& rResidualGradient,
                                           Vector& rResponseGradient,
                                           const ProcessInfo& rProcessInfo) override
    {
        ComputePointTemperatureSensitivityContribution(rResidualGradient, rAdjointElement.GetGeometry().Points(),rResponseGradient);
    }

    void CalculateSecondDerivativesGradient(const Element& rAdjointElement,
                                            const Matrix& rResidualGradient,
                                            Vector& rResponseGradient,
                                            const ProcessInfo& rProcessInfo) override
    {
        ComputePointTemperatureSensitivityContribution(rResidualGradient, rAdjointElement.GetGeometry().Points(),rResponseGradient);
    }

    void CalculatePartialSensitivity(Element& rAdjointElement,
                                     const Variable<array_1d<double, 3>>& rVariable,
                                     const Matrix& rSensitivityMatrix,
                                     Vector& rSensitivityGradient,
                                     const ProcessInfo& rProcessInfo) override
    {
        ComputePointTemperatureSensitivityContribution(rSensitivityMatrix, rAdjointElement.GetGeometry().Points(),rSensitivityGradient);
    }

    double CalculateValue(ModelPart& rModelPart) override
    {
        KRATOS_TRY;
        KRATOS_ERROR
            << "PointTemperature::CalculateValue(ModelPart& rModelPart) is not implemented!!!\n";
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

    ///@}
    ///@name Private Operators
    ///@{

    void ComputePointTemperatureSensitivityContribution(
        const Matrix& rDerivativesOfResidual,
        const Element::NodesArrayType& rNodes,
        Vector& rLocalSensitivityContribution) const
    {
        constexpr std::size_t max_size = 50;
        BoundedVector<double, max_size> nodal_flag_vector(rDerivativesOfResidual.size2());

        const unsigned int num_nodes = rNodes.size();
        for (unsigned int i_node = 0; i_node < num_nodes; ++i_node)
        {
            nodal_flag_vector[i_node] = rNodes[i_node].Is(STRUCTURE) ? 1.0 : 0.0;
        }

        if (rLocalSensitivityContribution.size() != rDerivativesOfResidual.size1())
            rLocalSensitivityContribution.resize(rDerivativesOfResidual.size1(), false);

        noalias(rLocalSensitivityContribution) = prod(rDerivativesOfResidual, nodal_flag_vector);
    }

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
};

///@} // Kratos Classes

///@} // ConvectionDiffusionApplication group

}

#endif // KRATOS_POINT_TEMPERATURE_RESPONSE_FUNCTION_H_INCLUDED