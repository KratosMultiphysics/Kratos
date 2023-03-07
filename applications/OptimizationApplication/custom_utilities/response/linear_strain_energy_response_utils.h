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

class KRATOS_API(OPTIMIZATION_APPLICATION) LinearStrainEnergyResponseUtils
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

    static double CalculateValue(const std::vector<ModelPart*>& rModelParts);

    static void CalculateSensitivity(
        const std::vector<ModelPart*>& rEvaluatedModelParts,
        const SensitivityModelPartVariablesListMap& rSensitivityModelPartVariableInfo,
        const double PerturbationSize);

    ///@}
private:
    ///@name Private static operations
    ///@{

    static double CalculateModelPartValue(ModelPart& rModelPart);

    static void CalculateStrainEnergyShapeSensitivity(
        ModelPart& rModelPart,
        const double Delta,
        const Variable<array_1d<double, 3>>& rOutputSensitivityVariable);

    static void CalculateStrainEnergyLinearlyDependentPropertySensitivity(
        ModelPart& rModelPart,
        const Variable<double>& rPrimalVariable,
        const Variable<double>& rOutputSensitivityVariable);

    static void CalculateStrainEnergyFiniteDifferencePropertySensitivity(
        ModelPart& rModelPart,
        const double Delta,
        const Variable<double>& rPrimalVariable,
        const Variable<double>& rOutputSensitivityVariable);

    ///@}
};

///@}
}