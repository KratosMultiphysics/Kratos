//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    
//

#if !defined(KRATOS_NODE_DISPLACEMENT_RESPONSE_FUNCTION_H_INCLUDED)
#define KRATOS_NODE_DISPLACEMENT_RESPONSE_FUNCTION_H_INCLUDED

// System includes
#include <algorithm>

// External includes

// Project includes
#include "includes/kratos_parameters.h"
#include "response_functions/adjoint_response_function.h"

namespace Kratos
{
///@name Kratos Classes
///@{

class NodeDisplacementResponseFunction : public AdjointResponseFunction
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(NodeDisplacementResponseFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    explicit NodeDisplacementResponseFunction(Parameters Settings, ModelPart& rModelPart)
    : mrModelPart(rModelPart)
    {
        KRATOS_TRY;
        Parameters default_params(R"(
        {
            "node_id" : 0,
            "direction" : [1.0, 0.0, 0.0]
        })");
        Settings.RecursivelyValidateAndAssignDefaults(default_params);
        mNodeId = Settings["node_id"].GetInt();
        mDirection.resize(3);
        for (unsigned i = 0; i < 3; ++i)
            mDirection[i] = Settings["direction"][i].GetDouble();
        KRATOS_CATCH("");
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    double CalculateValue() override
    {
        KRATOS_TRY;
        const std::size_t domain_size = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];
        const auto& r_node = mrModelPart.GetNode(mNodeId);
        const array_1d<double, 3>& r_displacement = r_node.FastGetSolutionStepValue(DISPLACEMENT);
        double value{0.0};
        for (std::size_t d = 0; d < domain_size; ++d)
            value += r_displacement[d] * mDirection[d];
        return value;
        KRATOS_CATCH("");
    }

    void CalculateGradient(const Element& rAdjointElement,
                           const Matrix& rResidualGradient,
                           Vector& rResponseGradient,
                           const ProcessInfo& rProcessInfo) override
    {
        KRATOS_TRY;
        rResponseGradient = ZeroVector(rResidualGradient.size1());
        const auto& r_geom = rAdjointElement.GetGeometry();
        auto it = std::find_if(r_geom.begin(), r_geom.end(),
                               [this](const Element::NodeType& r_node) {
                                   return (r_node.Id() == mNodeId);
                               });

        if (it != r_geom.end()) // Found node in current element.
        {
            std::size_t ws_dim = rAdjointElement.WorkingSpaceDimension();
            const std::size_t pos = it - r_geom.begin();
            for (std::size_t d = 0; d < ws_dim; ++d) // TODO: Divide by number of neighbour elements.
                rResponseGradient[pos * ws_dim + d] = mDirection[d];
        }
        KRATOS_CATCH("");
    }

    void CalculateGradient(const Condition& rAdjointCondition,
                           const Matrix& rResidualGradient,
                           Vector& rResponseGradient,
                           const ProcessInfo& rProcessInfo) override
    {
        KRATOS_TRY;
        rResponseGradient = ZeroVector(rResidualGradient.size1());
        KRATOS_CATCH("");
    }

    void CalculatePartialSensitivity(Element& rAdjointElement,
                                     const Variable<array_1d<double, 3>>& rVariable,
                                     const Matrix& rSensitivityMatrix,
                                     Vector& rSensitivityGradient,
                                     const ProcessInfo& rProcessInfo) override
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    void CalculatePartialSensitivity(Condition& rAdjointCondition,
                                     const Variable<array_1d<double, 3>>& rVariable,
                                     const Matrix& rSensitivityMatrix,
                                     Vector& rSensitivityGradient,
                                     const ProcessInfo& rProcessInfo) override
    {
        KRATOS_TRY;
        rSensitivityGradient = ZeroVector(rSensitivityMatrix.size1());
        KRATOS_CATCH("");
    }

    ///@}
    private:
    ///@name Member Variables
    ///@{
        ModelPart& mrModelPart;
        std::size_t mNodeId = -1;
        Vector mDirection;
    ///@}
    ///@name Private Operations
    ///@{
        
    ///@}
};

///@} // Kratos Classes
} /* namespace Kratos.*/

#endif /* KRATOS_NODE_DISPLACEMENT_RESPONSE_FUNCTION_H_INCLUDED defined */
