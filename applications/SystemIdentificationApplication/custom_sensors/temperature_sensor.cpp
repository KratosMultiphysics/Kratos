//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         SystemIdentificationApplication/license.txt
//
//  Main authors:    Talhah Ansari
//                   Suneth Warnakulasuriya
//

// System includes
#include <cmath>
#include <limits>

// External includes

// Project includes
#include "utilities/brute_force_point_locator.h"

// Application includes
#include "custom_utilities/sensor_utils.h"
#include "system_identification_application_variables.h"

// Include base h
#include "temperature_sensor.h"

namespace Kratos {

/// Constructor.
TemperatureSensor::TemperatureSensor(
    const std::string& rName,
    Node::Pointer pNode,
    const Element& rElement,
    const double Weight)
    : BaseType(rName, pNode, Weight),
      mElementId(rElement.Id())
{
    const auto& r_geometry = rElement.GetGeometry();
    const auto& current_sensor_location = *(this->GetNode());

    Point local_point;
    if (r_geometry.IsInside(current_sensor_location, local_point)) {
        // point is within the geometry. Use shape function evaluations
        // from the element geometry to get the shape function values.
        r_geometry.ShapeFunctionsValues(mNs, local_point);
    } else {
        // resize the shape functions
        if (mNs.size() != r_geometry.size()) {
            mNs.resize(r_geometry.size(), false);
        }

        // then first clear the shape functions
        mNs.clear();

        // point can be either on a node or on an edge
        bool is_point_on_a_node{false}, is_point_on_an_edge{false};

        // check if the point is on a node.
        for (IndexType i_node = 0; i_node < r_geometry.size(); ++i_node) {
            if (norm_2(r_geometry[i_node] - current_sensor_location) < std::numeric_limits<double>::epsilon()) {
                mNs[i_node] = 1.0;
                is_point_on_a_node = true;
                break;
            }
        }

        // point is not on a node, then check if it is on an edge.
        if (!is_point_on_a_node) {
            for (IndexType i_node = 0; i_node < r_geometry.size(); ++i_node) {
                const auto& r_node_1 = r_geometry[i_node];
                const auto& r_node_2 = r_geometry[(i_node + 1) % r_geometry.size()];

                // compute the area of the triangle from the 3 points
                const double area = r_node_1.X() * (r_node_2.Y() - current_sensor_location.Y()) + r_node_2.X() * (current_sensor_location.Y() - r_node_1.Y()) + current_sensor_location.X() * (r_node_1.Y() - r_node_2.Y());
                if (std::abs(area) <= std::numeric_limits<double>::epsilon()) {
                    const double r = norm_2(current_sensor_location - r_node_1) / norm_2(r_node_2 - r_node_1);
                    mNs[i_node] = 1 - r;
                    mNs[(i_node + 1) % r_geometry.size()] = r;
                    is_point_on_an_edge = true;
                    break;
                }
            }
        }

        KRATOS_ERROR_IF_NOT(is_point_on_a_node || is_point_on_an_edge)
                << "The point " << this->GetNode()->Coordinates() << " is not inside or on the boundary of the geometry of element with id "
                << mElementId << ".";
    }

    this->GetNode()->SetValue(SENSOR_ELEMENT_ID, static_cast<int>(mElementId));
}

Sensor::Pointer TemperatureSensor::Create(
    ModelPart& rDomainModelPart,
    ModelPart& rSensorModelPart,
    const IndexType Id,
    Parameters SensorParameters)
{
    KRATOS_TRY

    SensorParameters.ValidateAndAssignDefaults(TemperatureSensor::GetDefaultParameters());

    const auto& location = SensorParameters["location"].GetVector();
    KRATOS_ERROR_IF_NOT(location.size() == 3)
        << "Location of the sensor \"" << SensorParameters["name"].GetString()
        << "\" should have 3 components. [ location = " << location << " ].\n";

    Point loc(location[0], location[1], location[2]);

    Vector dummy_shape_functions;

    const auto element_id = BruteForcePointLocator(rDomainModelPart).FindElement(loc, dummy_shape_functions);
    const auto& r_element = rDomainModelPart.GetElement(element_id);

    auto p_node = rSensorModelPart.CreateNewNode(Id, location[0], location[1], location[2]);

    auto p_sensor = Kratos::make_shared<TemperatureSensor>(
        SensorParameters["name"].GetString(),
        p_node,
        r_element,
        SensorParameters["weight"].GetDouble()
    );

    SensorUtils::ReadVariableData(p_sensor->GetNode()->GetData(), SensorParameters["variable_data"]);

    return p_sensor;

    KRATOS_CATCH("");
}
Parameters TemperatureSensor::GetSensorParameters() const
{
    auto parameters = BaseType::GetSensorParameters();
    parameters.AddString("type", "temperature_sensor");
    return parameters;
}

Parameters TemperatureSensor::GetDefaultParameters()
{
    return Parameters(R"(
    {
        "type"         : "temperature_sensor",
        "name"         : "",
        "value"        : 0,
        "location"     : [0.0, 0.0, 0.0],
        "weight"       : 1.0,
        "variable_data": {}
    })" );
}

double TemperatureSensor::CalculateValue(ModelPart& rModelPart)
{
    double temperature = 0.0;
    if (rModelPart.HasElement(mElementId)) {
        const auto& r_element = rModelPart.GetElement(mElementId);
        const auto& r_geometry = r_element.GetGeometry();

        for (IndexType i = 0; i < r_geometry.size(); ++i) {
            temperature += r_geometry[i].FastGetSolutionStepValue(TEMPERATURE) * mNs[i];
        }
    }
    return rModelPart.GetCommunicator().GetDataCommunicator().SumAll(temperature);

}

void TemperatureSensor::CalculateGradient(
    const Element& rAdjointElement,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    // temperature sensor does not affect the adjoint system
    SetVectorToZero(rResponseGradient, rResidualGradient.size1());
}

void TemperatureSensor::CalculateGradient(
    const Condition& rAdjointCondition,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    SetVectorToZero(rResponseGradient, rResidualGradient.size1());
}

void TemperatureSensor::CalculateFirstDerivativesGradient(
    const Element& rAdjointElement,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    SetVectorToZero(rResponseGradient, rResidualGradient.size1());
}

void TemperatureSensor::CalculateFirstDerivativesGradient(
    const Condition& rAdjointCondition,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    SetVectorToZero(rResponseGradient, rResidualGradient.size1());
}

void TemperatureSensor::CalculateSecondDerivativesGradient(
    const Element& rAdjointElement,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    SetVectorToZero(rResponseGradient, rResidualGradient.size1());
}

void TemperatureSensor::CalculateSecondDerivativesGradient(
    const Condition& rAdjointCondition,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    SetVectorToZero(rResponseGradient, rResidualGradient.size1());
}

void TemperatureSensor::CalculatePartialSensitivity(
    Element& rAdjointElement,
    const Variable<double>& rVariable,
    const Matrix& rSensitivityMatrix,
    Vector& rSensitivityGradient,
    const ProcessInfo& rProcessInfo)
{
    // temperature is part of the cost function
    KRATOS_TRY;

    if (rSensitivityGradient.size() != rSensitivityMatrix.size1()) {
        rSensitivityGradient.resize(rSensitivityMatrix.size1(), false);
    }

    rSensitivityGradient.clear();

    if (rVariable == TEMPERATURE) {
        if (rAdjointElement.Id() == mElementId) {
            const auto& r_geometry = rAdjointElement.GetGeometry();
            const IndexType block_size = rSensitivityGradient.size() / r_geometry.size();

            for (IndexType i = 0; i < r_geometry.size(); ++i) {
                // -1 because in structures the residual is opposite to that of fluids and requires a negative sign
                // TODO: Fix this negative in structural responses when redesigning the adjoint framework
                rSensitivityGradient[i * block_size] = -1.0 * mNs[i];
            }
        }
    }

    KRATOS_CATCH("");
}

void TemperatureSensor::CalculatePartialSensitivity(
    Condition& rAdjointCondition,
    const Variable<double>& rVariable,
    const Matrix& rSensitivityMatrix,
    Vector& rSensitivityGradient,
    const ProcessInfo& rProcessInfo)
{
    SetVectorToZero(rSensitivityGradient, rSensitivityMatrix.size1());
}

void TemperatureSensor::CalculatePartialSensitivity(
    Element& rAdjointElement,
    const Variable<array_1d<double, 3>>& rVariable,
    const Matrix& rSensitivityMatrix,
    Vector& rSensitivityGradient,
    const ProcessInfo& rProcessInfo)
{
    SetVectorToZero(rSensitivityGradient, rSensitivityMatrix.size1());
}

void TemperatureSensor::CalculatePartialSensitivity(
    Condition& rAdjointCondition,
    const Variable<array_1d<double, 3>>& rVariable,
    const Matrix& rSensitivityMatrix,
    Vector& rSensitivityGradient,
    const ProcessInfo& rProcessInfo)
{
    SetVectorToZero(rSensitivityGradient, rSensitivityMatrix.size1());
}

std::string TemperatureSensor::Info() const
{
    std::stringstream msg;
    msg << "TemperatureSensor " << this->GetName();
    return msg.str();
}

void TemperatureSensor::PrintInfo(std::ostream& rOStream) const
{
    rOStream << Info() << std::endl;
}

void TemperatureSensor::PrintData(std::ostream& rOStream) const
{
    rOStream << "    Element Id: " << mElementId << std::endl;
    Sensor::PrintData(rOStream);
}

void TemperatureSensor::SetVectorToZero(
    Vector& rVector,
    const IndexType Size)
{
    if (rVector.size() != Size) {
        rVector.resize(Size, false);
    }

    rVector.clear();
}

}; // namespace Kratos
