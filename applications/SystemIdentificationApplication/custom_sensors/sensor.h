//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         SystemIdentificationApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <variant>

// External includes

// Project includes
#include "containers/array_1d.h"
#include "containers/data_value_container.h"
#include "expression/container_expression.h"
#include "geometries/point.h"
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "response_functions/adjoint_response_function.h"

// Application includes

namespace Kratos {
///@name Kratos Classes
///@{

class KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) Sensor: public AdjointResponseFunction
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(Sensor);

    using ContainerExpressionType = std::variant<
                                        ContainerExpression<ModelPart::NodesContainerType>::Pointer,
                                        ContainerExpression<ModelPart::ConditionsContainerType>::Pointer,
                                        ContainerExpression<ModelPart::ElementsContainerType>::Pointer
                                        >;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    Sensor(
        const std::string& rName,
        Node::Pointer pNode,
        const double Weight);

    /// Destructor.
    ~Sensor() = default;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Get the Sensor parameters.
     *
     * This method is used retrieve sensor parameters used in this sensor.
     * These parameters are used to reconstruct the sensor from scratch
     * if required.
     *
     * @return const Parameters
     */
    virtual Parameters GetSensorParameters() const;

    /**
     * @brief Calculate the current value of the sensor using the given model part.
     *
     * @param rModelPart        Model part to calculate the sensor value.
     * @return double           Calculated sensor value.
     */
    double CalculateValue(ModelPart& rModelPart) override
    {
        KRATOS_ERROR << "Calling base class Sensor::CalculateValue. Please implement it in the derrived class.";
    }

    /**
     * @brief Get the Type of the sensor
     *
     * @return std::string
     */
    std::string GetName() const;

    /**
     * @brief Get the sensor location.
     *
     * @return Point
     */
    Node::Pointer GetNode() const;

    /**
     * @brief Get the Weight of the sensor.
     *
     * @return double       The weight of the sensor .
     */
    double GetWeight() const;

    /**
     * @brief Get the Sensor value
     *
     * @return double       Value of the sensor.
     */
    double GetSensorValue() const;

    /**
     * @brief Set the Sensor Value
     *
     * @param Value         Value to be set.
     */
    void SetSensorValue(const double Value);

    /**
     * @brief Adds given container expression
     *
     * @throws If an expression is found already under the given name.
     *
     * @param rExpressionName       Name of the container expression.
     * @param pContainerExpression  Container expression pointer to be added.
     */
    void AddContainerExpression(
        const std::string& rExpressionName,
        ContainerExpressionType pContainerExpression);

    /**
     * @brief Get the Nodal Expression for specified expression name
     *
     * @throws If the @ref rExpressionName is not found in the map of container expressions.
     *
     * @param rExpressionName            Expression name
     * @return ContainerExpressionType   Container expression
     */
    ContainerExpressionType GetContainerExpression(const std::string& rExpressionName) const;

    /**
     * @brief Get the Container Expressions map
     *
     * @return std::unordered_map<std::string, ContainerExpressionType>   Container expressions map
     */
    std::unordered_map<std::string, ContainerExpressionType> GetContainerExpressionsMap() const;

    /**
     * @brief Clear container expressions from the sensors.
     *
     */
    void ClearContainerExpressions();

    ///@}
    ///@name Input and output
    ///@{

    virtual std::string Info() const;

    virtual void PrintInfo(std::ostream& rOStream) const;

    virtual void PrintData(std::ostream& rOStream) const;

    ///@}

private:
    ///@name Private member variables
    ///@{

    const std::string mName;

    Node::Pointer mpNode;

    const double mWeight;

    double mSensorValue;

    std::unordered_map<std::string, ContainerExpressionType> mContainerExpressions;

    ///@}
};

///@} // Kratos Classes

/// output stream functions
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const Sensor& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} /* namespace Kratos.*/