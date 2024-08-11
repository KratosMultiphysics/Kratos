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
#include "containers/data_value_container.h"
#include "expression/container_expression.h"
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
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

    using ContainerExpressionPointerType = std::variant<
                                            ContainerExpression<ModelPart::NodesContainerType>::Pointer,
                                            ContainerExpression<ModelPart::ConditionsContainerType>::Pointer,
                                            ContainerExpression<ModelPart::ElementsContainerType>::Pointer>;

    KRATOS_CLASS_POINTER_DEFINITION(Sensor);

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
     * @brief Get the Node representing the sensor
     *
     * @return Node representing the sensor
     */
    Node::Pointer GetNode() const;

    /**
     * @brief Get the Sensor parameters.
     *
     * This method is used retrieve sensor parameters used in this sensor.
     * These parameters are used to reconstruct the sensor from scratch
     * if required.
     *
     * @return const Parameters
     */
    virtual const Parameters GetSensorParameters() const
    {
        KRATOS_ERROR << "Calling base class Sensor::GetSensorParameters. Please implement it in the derived class.";
        return Parameters(R"({})" );
    }

    /**
     * @brief Calculate the current value of the sensor using the given model part.
     *
     * @param rModelPart        Model part to calculate the sensor value.
     * @return double           Calculated sensor value.
     */
    double CalculateValue(ModelPart& rModelPart) override
    {
        KRATOS_ERROR << "Calling base class Sensor::CalculateValue. Please implement it in the derived class.";
    }

    /**
     * @brief Get the Type of the sensor
     *
     * @return std::string
     */
    std::string GetName() const;

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
     * @throws If an expression is found already under the given name in the respective TContainerType container.
     *
     * @tparam TContainerType       Type of the container expression.
     * @param rExpressionName       Name of the container expression.
     */
    void AddContainerExpression(
        const std::string& rExpressionName,
        ContainerExpressionPointerType pContainerExpression);

    /**
     * @brief Get the Container Expression for specified expression name
     *
     * @throws If the @ref rExpressionName is not found in the map of nodal expressions.
     *
     * @param rExpressionName                   Expression name
     * @return ContainerExpressionPointerType   Container expression
     */
    ContainerExpressionPointerType GetContainerExpression(const std::string& rExpressionName) const;

    /**
     * @brief Get the Container Expressions map
     *
     * @return std::unordered_map<std::string, ContainerExpressionPointerType> Container expressions map
     */
    std::unordered_map<std::string, ContainerExpressionPointerType> GetContainerExpressionsMap() const;

    /**
     * @brief Clear expressions from the sensor.
     *
     */
    void Clear();

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

    std::unordered_map<std::string, ContainerExpressionPointerType> mContainerExpressionsMap;

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