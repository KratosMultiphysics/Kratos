//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: DigitalTwinApplication/license.txt
//
//  Main authors:    Suneth Wranakulasuriya
//

#pragma once

// System includes
#include <string>
#include <unordered_map>

// External includes

// Project includes
#include "containers/data_value_container.h"
#include "expression/container_expression.h"
#include "geometries/point.h"
#include "includes/define.h"
#include "includes/indexed_object.h"
#include "includes/model_part.h"

// Application includes

namespace Kratos {
///@name Kratos Classes
///@{

class SensorSpecification: public IndexedObject, public DataValueContainer
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = IndexedObject::IndexType;

    KRATOS_CLASS_POINTER_DEFINITION(SensorSpecification);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    SensorSpecification(
        const std::string &rName,
        const IndexType NewId);

    /// Destructor.
    ~SensorSpecification() override = default;

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief Get the Name of the sensor
     *
     * @return std::string
     */
    std::string GetName() const;

    /**
     * @brief Set the Location of the sensor
     *
     * @param rLocation     New location to set.
     */
    void SetLocation(const array_1d<double, 3>& rLocation);

    /**
     * @brief Get the sensor location.
     *
     * @return Point
     */
    array_1d<double, 3> GetLocation() const;

    /**
     * @brief Set the Sensor Value
     *
     * @param SensorValue   Sensor value to be set.
     */
    void SetSensorValue(const double SensorValue);

    /**
     * @brief Get the Sensor value
     *
     * @return double
     */
    double GetSensorValue() const;

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

    void ClearNodalExpressions();

    void ClearConditionExpressions();

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

    double mSensorValue;

    array_1d<double, 3> mLocation;

    std::unordered_map<std::string, ContainerExpression<ModelPart::NodesContainerType>::Pointer> mNodalExpressions;

    std::unordered_map<std::string, ContainerExpression<ModelPart::ConditionsContainerType>::Pointer> mConditionExpressions;

    std::unordered_map<std::string, ContainerExpression<ModelPart::ElementsContainerType>::Pointer> mElementExpressions;

    ///@}
};

/// output stream functions
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const SensorSpecification& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

///@} // Kratos Classes

} /* namespace Kratos.*/