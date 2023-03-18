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
#include <unordered_map>
#include <variant>
#include <vector>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

// Application includes

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

    using SensitivityFieldVariableTypes = std::variant<const Variable<double>*, const Variable<array_1d<double, 3>>*>;

    using SensitivityVariableModelPartsListMap = std::unordered_map<SensitivityFieldVariableTypes, std::vector<ModelPart*>>;

    ///@}
    ///@name Static operations
    ///@{

    static double CalculateValue(const std::vector<ModelPart*>& rModelParts);

    static void CalculateSensitivity(
        ModelPart& rAnalysisModelPart,
        const std::vector<ModelPart*>& rEvaluatedModelParts,
        const SensitivityVariableModelPartsListMap& rSensitivityVariableModelPartInfo,
        const double PerturbationSize);

    ///@}
private:
    ///@name Private static operations
    ///@{

    template<class TEntityType>
    static double CalculateEntityStrainEnergy(
        TEntityType& rEntity,
        Matrix& rLHS,
        Vector& rRHS,
        Vector& rX,
        const ProcessInfo& rProcessInfo);

    static double CalculateModelPartValue(ModelPart& rModelPart);

    template<class TEntityType>
    static void CalculateStrainEnergyEntitySemiAnalyticShapeSensitivity(
        TEntityType& rEntity,
        Vector& rX,
        Vector& rRefRHS,
        Vector& rPerturbedRHS,
        typename TEntityType::Pointer& pThreadLocalEntity,
        ModelPart& rModelPart,
        std::vector<std::string>& rModelPartNames,
        const double Delta,
        const IndexType MaxNodeId,
        const Variable<array_1d<double, 3>>& rOutputSensitivityVariable);

    static void CalculateStrainEnergySemiAnalyticShapeSensitivity(
        ModelPart& rModelPart,
        const double Delta,
        const Variable<array_1d<double, 3>>& rOutputSensitivityVariable);

    static void CalculateStrainEnergyLinearlyDependentPropertySensitivity(
        ModelPart& rModelPart,
        const Variable<double>& rPrimalVariable,
        const Variable<double>& rOutputSensitivityVariable);

    static void CalculateStrainEnergySemiAnalyticPropertySensitivity(
        ModelPart& rModelPart,
        const double Delta,
        const Variable<double>& rPrimalVariable,
        const Variable<double>& rOutputSensitivityVariable);

    ///@}
};

///@}
}