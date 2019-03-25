//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors: Martin Fusseder, https://github.com/MFusseder
//                Michael Andre, https://github.com/msandre
//
//

#if !defined(KRATOS_SENSITIVITY_BUILDER_H_INCLUDED)
#define KRATOS_SENSITIVITY_BUILDER_H_INCLUDED

// System includes
#include <string>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "response_functions/adjoint_response_function.h"

namespace Kratos
{
///@name Kratos Classes
///@{

class KRATOS_API(KRATOS_CORE) SensitivityBuilder
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(SensitivityBuilder);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    SensitivityBuilder(Parameters Settings,
                       ModelPart& rModelPart,
                       AdjointResponseFunction::Pointer pResponseFunction);

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    static void CalculateNodalSolutionStepSensitivities(const std::vector<std::string>& rVariables,
                                                        ModelPart& rModelPart,
                                                        AdjointResponseFunction& rResponseFunction,
                                                        double ScalingFactor);

    static void CalculateNonHistoricalSensitivities(const std::vector<std::string>& rVariables,
                                                    ModelPart::ElementsContainerType& rElements,
                                                    AdjointResponseFunction& rResponseFunction,
                                                    const ProcessInfo& rProcessInfo,
                                                    double ScalingFactor);

    static void CalculateNonHistoricalSensitivities(const std::vector<std::string>& rVariables,
                                                    ModelPart::ConditionsContainerType& rConditions,
                                                    AdjointResponseFunction& rResponseFunction,
                                                    const ProcessInfo& rProcessInfo,
                                                    double ScalingFactor);

    void Initialize();

    void UpdateSensitivities();

    void Clear();

    /// Clear the flags which are indicating the membership of a node, element or condition in the sensitivity model part
    void ClearFlags();

    /// Clear sensitivities in historical and non-historical database
    void ClearSensitivities();

    ///@}

private:
    ///@name Member Variables
    ///@{

    ModelPart* mpModelPart = nullptr;
    ModelPart* mpSensitivityModelPart = nullptr;
    AdjointResponseFunction::Pointer mpResponseFunction;
    std::vector<std::string> mNodalSolutionStepSensitivityVariables;
    std::vector<std::string> mElementDataValueSensitivityVariables;
    std::vector<std::string> mConditionDataValueSensitivityVariables;
    std::string mBuildMode = "static";
    bool mNodalSolutionStepSensitivityCalculationIsThreadSafe = false;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
};

///@} // Kratos Classes

} /* namespace Kratos.*/

#endif /* KRATOS_SENSITIVITY_BUILDER_H_INCLUDED defined */
