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

// System includes

// External includes

// Project includes
#include "includes/kratos_components.h"

// Application includes

// Include base h
#include "sensor_utils.h"

namespace Kratos {

bool SensorUtils::IsPointInGeometry(
    const Point& rPoint,
    const Geometry<Node>& rGeometry)
{
    Point result;
    return rGeometry.IsInside(rPoint, result);
}

SensorUtils::SensorViewType SensorUtils::CreateSensorView(
    Sensor::Pointer pSensor,
    const std::string& rExpressionName)
{
    KRATOS_TRY

    auto p_expression = pSensor->GetContainerExpression(rExpressionName);
    if (std::holds_alternative<ContainerExpression<ModelPart::NodesContainerType>::Pointer>(p_expression)) {
        return std::make_shared<SensorView<ModelPart::NodesContainerType>>(pSensor, rExpressionName);
    } else if (std::holds_alternative<ContainerExpression<ModelPart::ConditionsContainerType>::Pointer>(p_expression)) {
        return std::make_shared<SensorView<ModelPart::ConditionsContainerType>>(pSensor, rExpressionName);
    } else if (std::holds_alternative<ContainerExpression<ModelPart::ElementsContainerType>::Pointer>(p_expression)) {
        return std::make_shared<SensorView<ModelPart::ElementsContainerType>>(pSensor, rExpressionName);
    } else {
        KRATOS_ERROR << "Unsupported expression type.";
        return std::make_shared<SensorView<ModelPart::ElementsContainerType>>(pSensor, rExpressionName);
    }

    KRATOS_CATCH("");
}

void SensorUtils::ReadVariableData(
    DataValueContainer& rDataValueContainer,
    Parameters VariableDataParameters)
{
    KRATOS_TRY

    for (auto var_data_itr = VariableDataParameters.begin(); var_data_itr != VariableDataParameters.end(); ++var_data_itr ) {
        const auto& var_name = var_data_itr.name();
        if (KratosComponents<Variable<bool>>::Has(var_name)) {
            rDataValueContainer.SetValue(KratosComponents<Variable<bool>>::Get(var_name), var_data_itr->GetBool());
        } else if (KratosComponents<Variable<int>>::Has(var_name)) {
            rDataValueContainer.SetValue(KratosComponents<Variable<int>>::Get(var_name), var_data_itr->GetInt());
        } else if (KratosComponents<Variable<double>>::Has(var_name)) {
            rDataValueContainer.SetValue(KratosComponents<Variable<double>>::Get(var_name), var_data_itr->GetDouble());
        } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(var_name)) {
            const auto& r_value = var_data_itr->GetVector();

            KRATOS_ERROR_IF(r_value.size() != 3)
                << "The variable \"" << var_name
                << "\" is of type Array3D, but the given vector size = "
                << r_value.size() << " which does not match with required 3 components vector. [ value = "
                << r_value << " ].\n";

            rDataValueContainer.SetValue(KratosComponents<Variable<array_1d<double, 3>>>::Get(var_name),
                                         array_1d<double, 3>{r_value[0], r_value[1], r_value[2]});
        } else if (KratosComponents<Variable<array_1d<double, 4>>>::Has(var_name)) {
            const auto& r_value = var_data_itr->GetVector();

            KRATOS_ERROR_IF(r_value.size() != 4)
                << "The variable \"" << var_name
                << "\" is of type Array4D, but the given vector size = "
                << r_value.size() << " which does not match with required 4 components vector. [ value = "
                << r_value << " ].\n";

            rDataValueContainer.SetValue(KratosComponents<Variable<array_1d<double, 4>>>::Get(var_name),
                                         array_1d<double, 4>{r_value[0], r_value[1], r_value[2], r_value[3]});
        } else if (KratosComponents<Variable<array_1d<double, 4>>>::Has(var_name)) {
            const auto& r_value = var_data_itr->GetVector();

            KRATOS_ERROR_IF(r_value.size() != 4)
                << "The variable \"" << var_name
                << "\" is of type Array4D, but the given vector size = "
                << r_value.size() << " which does not match with required 4 components vector. [ value = "
                << r_value << " ].\n";

            rDataValueContainer.SetValue(KratosComponents<Variable<array_1d<double, 4>>>::Get(var_name),
                                         array_1d<double, 4>{r_value[0], r_value[1], r_value[2], r_value[3]});
        } else if (KratosComponents<Variable<array_1d<double, 4>>>::Has(var_name)) {
            const auto& r_value = var_data_itr->GetVector();

            KRATOS_ERROR_IF(r_value.size() != 4)
                << "The variable \"" << var_name
                << "\" is of type Array4D, but the given vector size = "
                << r_value.size() << " which does not match with required 4 components vector. [ value = "
                << r_value << " ].\n";

            rDataValueContainer.SetValue(KratosComponents<Variable<array_1d<double, 4>>>::Get(var_name),
                                         array_1d<double, 4>{r_value[0], r_value[1], r_value[2], r_value[3]});
        } else if (KratosComponents<Variable<array_1d<double, 4>>>::Has(var_name)) {
            const auto& r_value = var_data_itr->GetVector();

            KRATOS_ERROR_IF(r_value.size() != 4)
                << "The variable \"" << var_name
                << "\" is of type Array4D, but the given vector size = "
                << r_value.size() << " which does not match with required 4 components vector. [ value = "
                << r_value << " ].\n";

            rDataValueContainer.SetValue(KratosComponents<Variable<array_1d<double, 4>>>::Get(var_name),
                                         array_1d<double, 4>{r_value[0], r_value[1], r_value[2], r_value[3]});
        } else if (KratosComponents<Variable<array_1d<double, 6>>>::Has(var_name)) {
            const auto& r_value = var_data_itr->GetVector();

            KRATOS_ERROR_IF(r_value.size() != 6)
                << "The variable \"" << var_name
                << "\" is of type Array6D, but the given vector size = "
                << r_value.size() << " which does not match with required 6 components vector. [ value = "
                << r_value << " ].\n";

            rDataValueContainer.SetValue(KratosComponents<Variable<array_1d<double, 6>>>::Get(var_name),
                                         array_1d<double, 6>{r_value[0], r_value[1], r_value[2], r_value[3], r_value[4], r_value[5]});
        } else if (KratosComponents<Variable<array_1d<double, 6>>>::Has(var_name)) {
            const auto& r_value = var_data_itr->GetVector();

            KRATOS_ERROR_IF(r_value.size() != 6)
                << "The variable \"" << var_name
                << "\" is of type Array6D, but the given vector size = "
                << r_value.size() << " which does not match with required 6 components vector. [ value = "
                << r_value << " ].\n";

            rDataValueContainer.SetValue(KratosComponents<Variable<array_1d<double, 6>>>::Get(var_name),
                                         array_1d<double, 6>{r_value[0], r_value[1], r_value[2], r_value[3], r_value[4], r_value[5]});
        } else if (KratosComponents<Variable<array_1d<double, 9>>>::Has(var_name)) {
            const auto& r_value = var_data_itr->GetVector();

            KRATOS_ERROR_IF(r_value.size() != 9)
                << "The variable \"" << var_name
                << "\" is of type Array9D, but the given vector size = "
                << r_value.size() << " which does not match with required 9 components vector. [ value = "
                << r_value << " ].\n";

            rDataValueContainer.SetValue(KratosComponents<Variable<array_1d<double, 9>>>::Get(var_name),
                                        array_1d<double, 9>{r_value[0], r_value[1], r_value[2], r_value[3], r_value[4], r_value[5], r_value[6], r_value[7], r_value[8]});
        } else if (KratosComponents<Variable<Vector>>::Has(var_name)) {
            rDataValueContainer.SetValue(KratosComponents<Variable<Vector>>::Get(var_name), var_data_itr->GetVector());
        } else if (KratosComponents<Variable<Matrix>>::Has(var_name)) {
            rDataValueContainer.SetValue(KratosComponents<Variable<Matrix>>::Get(var_name), var_data_itr->GetMatrix());
        } else if (KratosComponents<Variable<std::string>>::Has(var_name)) {
            rDataValueContainer.SetValue(KratosComponents<Variable<std::string>>::Get(var_name), var_data_itr->GetString());
        } else {
            KRATOS_ERROR << "Unsupported variable type for variable name = \"" << var_name << "\"."
                         << "Supported variable types are:"
                         << "\n\tbool"
                         << "\n\tint"
                         << "\n\tdouble"
                         << "\n\tarray_1d<double, 3>"
                         << "\n\tarray_1d<double, 4>"
                         << "\n\tarray_1d<double, 6>"
                         << "\n\tarray_1d<double, 9>"
                         << "\n\tVector"
                         << "\n\tMatrix"
                         << "\n\tstd::string"
                         << "\nProvided parameter data = " << *var_data_itr << "\n";
        }
    }

    KRATOS_CATCH("");
}

} /* namespace Kratos.*/