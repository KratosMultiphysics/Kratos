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

// External includes

// Project includes
#include "containers/array_1d.h"
#include "containers/data_value_container.h"
#include "expression/container_expression.h"
#include "geometries/point.h"
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "response_functions/adjoint_response_function.h"

// Application includes

namespace Kratos {
///@name Kratos Classes
///@{

class Sensor: public AdjointResponseFunction, public DataValueContainer
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(Sensor);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    Sensor(
        const std::string& rName,
        const Point& rLocation,
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
    virtual const Parameters GetSensorParameters() const
    {
        KRATOS_ERROR << "Calling base class Sensor::GetSensorParameters. Please implement it in the derrived class.";
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
    Point GetLocation() const;

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
     * @param pContainerExpression  Container expression to be added.
     */
    template<class TContainerType>
    void AddContainerExpression(
        const std::string& rExpressionName,
        typename ContainerExpression<TContainerType>::Pointer pContainerExpression);

    /**
     * @brief Adds a given nodal expression under the given name.
     * @throws If an expression is found already under the given name.
     *
     * @param rExpressionName       Name of the nodal expression.
     * @param pNodalExpression      Nodal expression to be added.
     */
    void AddNodalExpression(
        const std::string& rExpressionName,
        ContainerExpression<ModelPart::NodesContainerType>::Pointer pNodalExpression);

    /**
     * @brief Get the Nodal Expression for specified expression name
     *
     * @throws If the @ref rExpressionName is not found in the map of nodal expressions.
     *
     * @param rExpressionName                                                Expression name
     * @return ContainerExpression<ModelPart::NodesContainerType>::Pointer   Nodal container expression
     */
    ContainerExpression<ModelPart::NodesContainerType>::Pointer GetNodalExpression(const std::string& rExpressionName) const;

    /**
     * @brief Get the Nodal Expressions map
     *
     * @return std::unordered_map<std::string, ContainerExpression<ModelPart::NodesContainerType>::Pointer>   Nodal container expressions map
     */
    std::unordered_map<std::string, ContainerExpression<ModelPart::NodesContainerType>::Pointer> GetNodalExpressionsMap() const;

    /**
     * @brief Adds a given condition expression under the given name.
     * @throws If an expression is found already under the given name.
     *
     * @param rExpressionName           Name of the condition expression.
     * @param pConditionExpression      Condition expression to be added.
     */
    void AddConditionExpression(
        const std::string& rExpressionName,
        ContainerExpression<ModelPart::ConditionsContainerType>::Pointer pConditionExpression);

    /**
     * @brief Get the Condition Expression for specified expression name
     *
     * @throws If the @ref rExpressionName is not found in the map of condition expressions.
     *
     * @param rExpressionName                                                     Expression name
     * @return ContainerExpression<ModelPart::ConditionsContainerType>::Pointer   Condition container expression
     */
    ContainerExpression<ModelPart::ConditionsContainerType>::Pointer GetConditionExpression(const std::string& rExpressionName) const;

    /**
     * @brief Get the condition Expressions map
     *
     * @return std::unordered_map<std::string, ContainerExpression<ModelPart::ConditionsContainerType>::Pointer>   Condition container expressions map
     */
    std::unordered_map<std::string, ContainerExpression<ModelPart::ConditionsContainerType>::Pointer> GetConditionExpressionsMap() const;

    /**
     * @brief Adds a given element expression under the given name.
     * @throws If an expression is found already under the given name.
     *
     * @param rExpressionName           Name of the element expression.
     * @param pElementExpression      Element expression to be added.
     */
    void AddElementExpression(
        const std::string& rExpressionName,
        ContainerExpression<ModelPart::ElementsContainerType>::Pointer pElementExpression);

    /**
     * @brief Get the Element Expression for specified expression name
     *
     * @throws If the @ref rExpressionName is not found in the map of element expressions.
     *
     * @param rExpressionName                                                    Expression name
     * @return ContainerExpression<ModelPart::ElementsContainerType>::Pointer    Element container expression
     */
    ContainerExpression<ModelPart::ElementsContainerType>::Pointer GetElementExpression(const std::string& rExpressionName) const;

    /**
     * @brief Get the element expressions map
     *
     * @return std::unordered_map<std::string, ContainerExpression<ModelPart::ElementsContainerType>::Pointer>   Element container expressions map
     */
    std::unordered_map<std::string, ContainerExpression<ModelPart::ElementsContainerType>::Pointer> GetElementExpressionsMap() const;

    /**
     * @brief Get the list of variable names in the specification.
     *
     * @return std::vector<std::string>     List of variable names in the specification.
     */
    std::vector<std::string> GetDataVariableNames() const;

    /**
     * @brief Clear nodal expressions from the sensos.
     *
     */
    void ClearNodalExpressions();

    /**
     * @brief Clear condition expressions from the sensor.
     *
     */
    void ClearConditionExpressions();

    /**
     * @brief Clear element expressions from the sensor.
     *
     */
    void ClearElementExpressions();

    ///@}
    ///@name Input and output
    ///@{

    std::string Info() const override;

    void PrintInfo(std::ostream& rOStream) const override;

    void PrintData(std::ostream& rOStream) const override;

    ///@}

private:
    ///@name Private member variables
    ///@{

    const std::string mName;

    const Point mLocation;

    const double mWeight;

    double mSensorValue;

    std::unordered_map<std::string, ContainerExpression<ModelPart::NodesContainerType>::Pointer> mNodalExpressions;

    std::unordered_map<std::string, ContainerExpression<ModelPart::ConditionsContainerType>::Pointer> mConditionExpressions;

    std::unordered_map<std::string, ContainerExpression<ModelPart::ElementsContainerType>::Pointer> mElementExpressions;

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