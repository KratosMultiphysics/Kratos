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

class KRATOS_API(OPTIMIZATION_APPLICATION) ResponseUtils
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    using SensitivityFieldVariableTypes = OptimizationUtils::SensitivityFieldVariableTypes;

    using SensitivityModelPartVariablesListMap = OptimizationUtils::SensitivityModelPartVariablesListMap;

    ///@}
    ///@name Static operations
    ///@{

    static void CheckAndPrepareModelPartsForSensitivityComputation(
        const std::vector<ModelPart*>& rEvaluatedModelParts,
        const SensitivityModelPartVariablesListMap& rSensitivityModelPartVariableInfo,
        const Flags& rFlag,
        const std::vector<SensitivityFieldVariableTypes>& rUsedNodalSensitivityVariables);

    ///@}
};

///@}
}