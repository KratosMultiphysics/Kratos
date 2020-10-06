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

#if !defined(KRATOS_SENSITIVITY_BUILDER_H_INCLUDED)
#define KRATOS_SENSITIVITY_BUILDER_H_INCLUDED

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
#include "utilities/assemble_utilities.h"

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


    template <class TDataType>
    using THomogeneousSensitivityVariables = std::vector<SensitivityVariables<TDataType>>;

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
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
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

    /**
     * @brief Assigns entity derivatives to nodes based on a given constant weighting
     *
     * This method assigns derivatives of a value calculated on given entity to nodes. The derivatives given by entity
     * should have the following order in the matrix
     *
     * \[
     *      M\left(a, i\right) = \frac{\partial y^i}{\partial x^c_k}
     * \]
     *
     * Where $a = c * DerivativeDimension + k$, and $i$ is the value dimension, $c$ is derivative node and $k$ is derivative node's direction
     * When these derivatives are stored on nodes, following order is assumed.
     *      Firstly, derivatives w.r.t. itself
     *      Then in the order of rNeighbourNodeIdsMap[current_node.Id] vector
     *
     * rNeighbourNodeIdsMap can be generated using FindGlobalNodalNeighboursProcess
     *
     * @see FindGlobalNodalNeighboursProcess
     *
     * @tparam TContainerType           Container type
     * @param rModelPart                Model part to use entities and nodes
     * @param DerivativeDimension       Dimensionality of derivative variable
     * @param rDerivativeVariable       Matrix type derivative variable which holds derivatives in an entity.
     * @param rNeighbourNodeIdsMap      Neighbour node ids map for all the nodes in rModelPart
     * @param Weight                    Constant weighting
     * @param rFlag                     Flag to check in entities whether entity derivatives should be distributed to nodes or not
     * @param CheckValue                Flag check value
     */
    template <class TContainerType>
    static void AssignEntityDerivativesToNodes(
        ModelPart& rModelPart,
        const int DerivativeDimension,
        const Variable<Matrix>& rDerivativeVariable,
        const std::unordered_map<int, std::vector<int>>& rNeighbourNodeIdsMap,
        const double Weight,
        const Flags& rFlag,
        const bool CheckValue = true);

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
    SensitivityBuilderScheme::Pointer mpSensitivityBuilderScheme;

    TSensitivityVariables mNodalSolutionStepSensitivityVariablesList;
    TSensitivityVariables mElementDataValueSensitivityVariablesList;
    TSensitivityVariables mConditionDataValueSensitivityVariablesList;

    std::string mBuildMode = "static";
    bool mNodalSolutionStepSensitivityCalculationIsThreadSafe = false;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    static TSensitivityVariables GetVariableLists(const std::vector<std::string>& rVariableNames);

    ///@}
};

///@} // Kratos Classes

} /* namespace Kratos.*/

#endif /* KRATOS_SENSITIVITY_BUILDER_H_INCLUDED defined */
