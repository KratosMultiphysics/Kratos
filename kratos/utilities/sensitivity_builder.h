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
//                Suneth Warnakulasuriya, https://github.com/sunethwarna
//

#pragma once

// System includes
#include <string>
#include <vector>
#include <unordered_map>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "response_functions/adjoint_response_function.h"
#include "solving_strategies/schemes/sensitivity_builder_scheme.h"

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

    /**
     * @brief Contains the sensitivity design and output variables
     *
     * The design variable is passed to CalculateSensitivityMatrix()
     * and CalculatePartialSensitivity() when calculating the local
     * sensitivity contributions. The local sensitivities are assembled
     * to the output variable.
     *
     * Example 1:
     * SensitivityVariables<double> vars("THICKNESS");
     * vars.pDesignVariable->Name(); // "THICKNESS"
     * vars.pOutputVariable->Name(); // "THICKNESS_SENSITIVITY"
     *
     * Example 2:
     * SensitivityVariables<double> vars("THICKNESS_SENSITIVITY");
     * vars.pDesignVariable->Name(); // "THICKNESS_SENSITIVITY"
     * vars.pOutputVariable->Name(); // "THICKNESS_SENSITIVITY"
     */
    template <class TDataType>
    struct SensitivityVariables
    {
        const Variable<TDataType>* pDesignVariable = nullptr;
        const Variable<TDataType>* pOutputVariable = nullptr;

        explicit SensitivityVariables(const std::string& rName)
        {
            KRATOS_TRY;
            const std::string output_suffix = "_SENSITIVITY";
            pDesignVariable = &KratosComponents<Variable<TDataType>>::Get(rName);
            if (rName.size() > output_suffix.size() &&
                std::equal(output_suffix.rbegin(), output_suffix.rend(), rName.rbegin()))
            {
                pOutputVariable = pDesignVariable;
            }
            else
            {
                pOutputVariable =
                    &KratosComponents<Variable<TDataType>>::Get(rName + output_suffix);
            }
            KRATOS_CATCH("");
        }
    };


    /**
     * @brief This type holds list of variables for the sensitivity analysis
     *
     *      1. Design variable
     *      2. Output variable
     *
     * In here, both design variable and output variable should be of the same type,
     * hence they are homogeneous variables.
     *
     * @tparam TDataType
     */
    template <class TDataType>
    using THomogeneousSensitivityVariables = std::vector<SensitivityVariables<TDataType>>;

    /**
     * @brief This type holds list of homogeneous variable lists.
     *
     */
    using TSensitivityVariables = std::tuple<
                                        THomogeneousSensitivityVariables<double>,
                                        THomogeneousSensitivityVariables<array_1d<double, 3>>
                                        >;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    SensitivityBuilder(Parameters Settings,
                       ModelPart& rModelPart,
                       AdjointResponseFunction::Pointer pResponseFunction);

    /// Constructor with SensitivityBuilderScheme.
    SensitivityBuilder(Parameters Settings,
                       ModelPart& rModelPart,
                       AdjointResponseFunction::Pointer pResponseFunction,
                       SensitivityBuilderScheme::Pointer pSensitivityBuilderScheme);

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Set the Response Function
     *
     * This sets the response function used in the sensitivity builder. This
     * is useful in cases where the LHS of the adjoint problem does not change,
     * but the RHS changes due to change in the the response function. In these
     * cases, this allows re-use of the already constructed LHS with different
     * RHSs.
     *
     * @param pResponseFunction         New Response function to be set.
     */
    void SetResponseFunction(AdjointResponseFunction::Pointer pResponseFunction);

    void Initialize();

    void InitializeSolutionStep();

    void UpdateSensitivities();

    void FinalizeSolutionStep();

    void Finalize();

    void Clear();

    /// Clear the flags which are indicating the membership of a node, element or condition in the sensitivity model part
    void ClearFlags();

    /// Clear sensitivities in historical and non-historical database
    void ClearSensitivities();

    ///@}
    ///@name static Operations
    ///@{

    static void CalculateNodalSolutionStepSensitivities(
        const std::vector<std::string>& rVariables,
        ModelPart& rModelPart,
        AdjointResponseFunction& rResponseFunction,
        double ScalingFactor);

    static void CalculateNodalSolutionStepSensitivities(
        const TSensitivityVariables& rVariables,
        ModelPart& rModelPart,
        AdjointResponseFunction& rResponseFunction,
        SensitivityBuilderScheme& rSensitivityBuilderScheme,
        double ScalingFactor);

    static void CalculateNonHistoricalSensitivities(
        const std::vector<std::string>& rVariables,
        ModelPart::ElementsContainerType& rElements,
        AdjointResponseFunction& rResponseFunction,
        const ProcessInfo& rProcessInfo,
        double ScalingFactor);

    static void CalculateNonHistoricalSensitivities(
        const TSensitivityVariables& rVariables,
        ModelPart::ElementsContainerType& rElements,
        AdjointResponseFunction& rResponseFunction,
        SensitivityBuilderScheme& rSensitivityBuilderScheme,
        const ProcessInfo& rProcessInfo,
        double ScalingFactor);

    static void CalculateNonHistoricalSensitivities(
        const std::vector<std::string>& rVariables,
        ModelPart::ConditionsContainerType& rConditions,
        AdjointResponseFunction& rResponseFunction,
        const ProcessInfo& rProcessInfo,
        double ScalingFactor);

    static void CalculateNonHistoricalSensitivities(
        const TSensitivityVariables& rVariables,
        ModelPart::ConditionsContainerType& rConditions,
        AdjointResponseFunction& rResponseFunction,
        SensitivityBuilderScheme& rSensitivityBuilderScheme,
        const ProcessInfo& rProcessInfo,
        double ScalingFactor);

private:
    ///@name Member Variables
    ///@{

    ModelPart* mpModelPart = nullptr;
    ModelPart* mpSensitivityModelPart = nullptr;
    AdjointResponseFunction::Pointer mpResponseFunction;
    SensitivityBuilderScheme::Pointer mpSensitivityBuilderScheme;

    TSensitivityVariables mNodalSolutionStepSensitivityVariablesList;
    TSensitivityVariables mElementDataValueSensitivityVariablesList;
    TSensitivityVariables mConditionDataValueSensitivityVariablesList;

    std::string mBuildMode = "static";
    bool mNodalSolutionStepSensitivityCalculationIsThreadSafe = false;

    ///@}
    ///@name Private Operations
    ///@{

    static TSensitivityVariables GetVariableLists(const std::vector<std::string>& rVariableNames);

    ///@}
};

///@} // Kratos Classes

} /* namespace Kratos.*/