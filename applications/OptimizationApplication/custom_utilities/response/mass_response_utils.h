//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Reza Najian Asl,
//                   Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <vector>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

// Application includes
#include "custom_utilities/optimization_utils.h"

namespace Kratos
{

///@name Kratos Classes
///@{

class KRATOS_API(OPTIMIZATION_APPLICATION) MassResponseUtils
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    using GeometryType = ModelPart::ElementType::GeometryType;

    using SensitivityFieldVariableTypes = OptimizationUtils::SensitivityFieldVariableTypes;

    using SensitivityModelPartVariablesListMap = OptimizationUtils::SensitivityModelPartVariablesListMap;

    ///@}
    ///@name Static operations
    ///@{

    static void Check(const std::vector<ModelPart const*>& rModelParts);

    static double CalculateValue(const std::vector<ModelPart const*>& rModelParts);

    static void CalculateSensitivity(
        const std::vector<ModelPart*>& rEvaluatedModelParts,
        const SensitivityModelPartVariablesListMap& rSensitivityModelPartVariableInfo);

    ///@}
private:
    ///@name Private operations
    ///@{

    static void CheckModelPart(const ModelPart& rModelPart);

    static double CalculateModelPartValue(const ModelPart& rModelPart);

    static bool HasVariableInProperties(
        const ModelPart& rModelPart,
        const Variable<double>& rVariable);

    static void CalculateMassGeometricalPropertySensitivity(
        ModelPart& rModelPart,
        const Variable<double>& rGeometricalPropertySensitivityVariable,
        const Variable<double>& rGeometricalCoflictingPropertySensitivityVariable,
        const Variable<double>& rOutputSensitivityVariable);

    static void CalculateMassShapeSensitivity(
        ModelPart& rModelPart,
        const Variable<array_1d<double, 3>>& rOutputSensitivityVariable);

    static void CalculateMassDensitySensitivity(
        ModelPart& rModelPart,
        const Variable<double>& rOutputSensitivityVariable);

    static void CalculateMassThicknessSensitivity(
        ModelPart& rModelPart,
        const Variable<double>& rOutputSensitivityVariable);

    static void CalculateMassCrossAreaSensitivity(
        ModelPart& rModelPart,
        const Variable<double>& rOutputSensitivityVariable);

    ///@}
};

///@}
}