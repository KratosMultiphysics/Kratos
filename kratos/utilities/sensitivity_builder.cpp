#include "utilities/sensitivity_builder.h"

#include "includes/kratos_parameters.h"
#include "includes/ublas_interface.h"
#include "input_output/logger.h"
#include "utilities/openmp_utils.h"
#include "utilities/variable_utils.h"
#include <algorithm>
#include <utility>

namespace
{
namespace sensitivity_builder_cpp // cotire guard
{
using namespace Kratos;

bool RequiresUpdateSensitivities(const Geometry<Node<3>>& rGeom)
{
    for (auto& r_node : rGeom)
    {
        if (r_node.GetValue(UPDATE_SENSITIVITIES))
        {
            return true;
        }
    }
    return false;
}

void AssembleOnDataValueContainer(const Variable<double>& rVariable,
                                  const Vector& rSrc,
                                  DataValueContainer& rDest)
{
    KRATOS_ERROR_IF(rSrc.size() != 1) << "Variable: " << rVariable.Name()
                                      << ", rSrc.size() = " << rSrc.size() << std::endl;
    rDest[rVariable] += rSrc[0];
}

void AssembleOnDataValueContainer(const Variable<array_1d<double, 3>>& rVariable,
                                  const Vector& rSrc,
                                  DataValueContainer& rDest)
{
    KRATOS_ERROR_IF(rSrc.size() > 3) << "Variable: " << rVariable.Name()
                                     << ", rSrc.size() = " << rSrc.size() << std::endl;
    array_1d<double, 3>& r_dest = rDest[rVariable];
    for (std::size_t d = 0; d < rSrc.size(); ++d)
    {
        r_dest[d] += rSrc[d];
    }
}

void AssembleNodalSolutionStepValues(const Variable<double>& rVariable,
                                     const Vector& rValues,
                                     Geometry<Node<3>>& rGeom)
{
    KRATOS_TRY;
    KRATOS_ERROR_IF(rGeom.size() != rValues.size())
        << "Geometry size: " << rGeom.size()
        << " is incompatible with vector size: " << rValues.size() << std::endl;
    std::size_t index = 0;
    for (auto& r_node : rGeom)
    {
        if (r_node.GetValue(UPDATE_SENSITIVITIES))
        {
            double& r_dest = r_node.FastGetSolutionStepValue(rVariable);
            const double rhs = rValues[index];
#pragma omp atomic
            r_dest += rhs;
        }
        ++index;
    }
    KRATOS_CATCH("");
}

void AssembleNodalSolutionStepValues(const Variable<array_1d<double, 3>>& rVariable,
                                     const Vector& rValues,
                                     Geometry<Node<3>>& rGeom)
{
    KRATOS_TRY;
    if (rGeom.size() * rGeom.WorkingSpaceDimension() != rValues.size())
    {
        KRATOS_ERROR << "Geometry size: " << rGeom.size()
                     << " and working space dimension: " << rGeom.WorkingSpaceDimension()
                     << " are incompatible with incompatible with vector size: "
                     << rValues.size() << std::endl;
    }
    std::size_t index = 0;
    for (auto& r_node : rGeom)
    {
        if (r_node.GetValue(UPDATE_SENSITIVITIES))
        {
            array_1d<double, 3>& r_value = r_node.FastGetSolutionStepValue(rVariable);
            for (std::size_t d = 0; d < rGeom.WorkingSpaceDimension(); ++d)
            {
                double& r_dest = r_value[d];
                const double rhs = rValues[index + d];
#pragma omp atomic
                r_dest += rhs;
            }
        }
        index += rGeom.WorkingSpaceDimension();
    }
    KRATOS_CATCH("");
}

struct LocalSensitivityBuilder
{
    Vector LocalSensitivity; /**< Local contribution to the total sensitivity */
    Vector PartialSensitivity; /**< Local contribution to the partial derivative of the response function */
    Vector AdjointVector; /**< Local adjoint vector */
    Matrix SensitivityMatrix; /**< Local contribution to the partial derivative of the residual */

    /// Calculate the local contributions of the sensitivity
    template <typename TDataType, typename TElement>
    LocalSensitivityBuilder& CalculateLocalSensitivity(const Variable<TDataType>& rVariable,
                                                       TElement& rElement,
                                                       AdjointResponseFunction& rResponseFunction,
                                                       const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY;
        rElement.CalculateSensitivityMatrix(rVariable, SensitivityMatrix, rProcessInfo);
        rElement.GetValuesVector(AdjointVector);
        KRATOS_ERROR_IF(AdjointVector.size() != SensitivityMatrix.size2())
            << "AdjointVector.size(): " << AdjointVector.size()
            << " incompatible with SensitivityMatrix.size1(): "
            << SensitivityMatrix.size1() << ". Variable: " << rVariable << std::endl;
        if (LocalSensitivity.size() != SensitivityMatrix.size1())
        {
            LocalSensitivity.resize(SensitivityMatrix.size1(), false);
        }
        rResponseFunction.CalculatePartialSensitivity(
            rElement, rVariable, SensitivityMatrix, PartialSensitivity, rProcessInfo);
        KRATOS_ERROR_IF(PartialSensitivity.size() != SensitivityMatrix.size1())
            << "PartialSensitivity.size(): " << PartialSensitivity.size()
            << " incompatible with SensitivityMatrix.size1(): "
            << SensitivityMatrix.size1() << ". Variable: " << rVariable << std::endl;
        noalias(LocalSensitivity) = prod(SensitivityMatrix, AdjointVector) + PartialSensitivity;
        return *this;
        KRATOS_CATCH("");
    }
};

/**
 * @brief Call a function on each container element in a parallel loop
 *
 * Any data associated with the function (e.g., functors or capture by value)
 * is copied to each thread.
 */
template <typename TContainer, typename TCallable>
void ParallelForEach(TContainer& rContainer, TCallable&& Fun)
{
#pragma omp parallel
    {
        auto local_fun(Fun);
#pragma omp for
        for (int i = 0; i < static_cast<int>(rContainer.size()); ++i)
        {
            std::forward<TCallable>(local_fun)(*(rContainer.begin() + i));
        }
    }
}

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

template <typename TDataType, typename TContainer>
void AssembleNodalSolutionStepContainerContributions(const SensitivityVariables<TDataType>& rVariables,
                                                     TContainer& rContainer,
                                                     AdjointResponseFunction& rResponseFunction,
                                                     const ProcessInfo& rProcessInfo,
                                                     double ScalingFactor)
{
    KRATOS_TRY;
    LocalSensitivityBuilder builder;
    ParallelForEach(rContainer, [&, builder](typename TContainer::data_type& rElement) mutable {
        auto& r_geom = rElement.GetGeometry();
        if (RequiresUpdateSensitivities(r_geom) == false)
            return;
        builder.CalculateLocalSensitivity(*rVariables.pDesignVariable, rElement,
                                          rResponseFunction, rProcessInfo);
        // if rElement does not contribute to local sensitivity, skip assembly
        if (builder.LocalSensitivity.size() != 0)
        {
            builder.LocalSensitivity *= ScalingFactor;
            AssembleNodalSolutionStepValues(*rVariables.pOutputVariable,
                                        builder.LocalSensitivity, r_geom);
        }
    });
    KRATOS_CATCH("");
}

template <typename TDataType>
void CalculateNodalSolutionStepSensitivities(const SensitivityVariables<TDataType>& rVariables,
                                             ModelPart& rModelPart,
                                             AdjointResponseFunction& rResponseFunction,
                                             double ScalingFactor)
{
    KRATOS_TRY;
    auto& r_comm = rModelPart.GetCommunicator();
    const auto& r_output_variable = *rVariables.pOutputVariable;
    if (r_comm.TotalProcesses() > 1)
    {
        // Make sure we only add the old sensitivity once when we assemble.
        ParallelForEach(rModelPart.Nodes(), [&r_output_variable, &r_comm](Node<3>& rNode) {
            if (rNode.FastGetSolutionStepValue(PARTITION_INDEX) != r_comm.MyPID())
                rNode.FastGetSolutionStepValue(r_output_variable) =
                    r_output_variable.Zero();
        });
    }
    AssembleNodalSolutionStepContainerContributions(
        rVariables, rModelPart.Elements(), rResponseFunction,
        rModelPart.GetProcessInfo(), ScalingFactor);
    AssembleNodalSolutionStepContainerContributions(
        rVariables, rModelPart.Conditions(), rResponseFunction,
        rModelPart.GetProcessInfo(), ScalingFactor);
    r_comm.AssembleCurrentData(r_output_variable);
    KRATOS_CATCH("");
}

void CalculateNodalSolutionStepSensitivities(const std::string& rVariable,
                                             ModelPart& rModelPart,
                                             AdjointResponseFunction& rResponseFunction,
                                             double ScalingFactor)
{
    KRATOS_TRY;
    if (KratosComponents<Variable<double>>::Has(rVariable) == true)
    {
        CalculateNodalSolutionStepSensitivities(
            SensitivityVariables<double>{rVariable}, rModelPart,
            rResponseFunction, ScalingFactor);
    }
    else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(rVariable) == true)
    {
        CalculateNodalSolutionStepSensitivities(
            SensitivityVariables<array_1d<double, 3>>{rVariable}, rModelPart,
            rResponseFunction, ScalingFactor);
    }
    else
    {
        KRATOS_ERROR << "Unsupported variable: " << rVariable << "." << std::endl;
    }
    KRATOS_CATCH("");
}

template <typename TDataType, typename TContainer>
void CalculateNonHistoricalSensitivities(const SensitivityVariables<TDataType>& rVariables,
                                         TContainer& rContainer,
                                         AdjointResponseFunction& rResponseFunction,
                                         const ProcessInfo& rProcessInfo,
                                         double ScalingFactor)
{
    KRATOS_TRY;
    LocalSensitivityBuilder builder;
    ParallelForEach(rContainer, [&, builder](typename TContainer::data_type& rElement) mutable {
        if (rElement.GetValue(UPDATE_SENSITIVITIES) == false)
            return;
        builder.CalculateLocalSensitivity(*rVariables.pDesignVariable, rElement,
                                          rResponseFunction, rProcessInfo);
        // if rElement does not contribute to local sensitivity, skip assembly
        if (builder.LocalSensitivity.size() != 0)
        {
            builder.LocalSensitivity *= ScalingFactor;
            AssembleOnDataValueContainer(*rVariables.pOutputVariable,
                                     builder.LocalSensitivity, rElement.Data());
        }
    });
    KRATOS_CATCH("");
}

template <typename TContainer>
void CalculateNonHistoricalSensitivities(const std::string& rVariable,
                                         TContainer& rContainer,
                                         AdjointResponseFunction& rResponseFunction,
                                         const ProcessInfo& rProcessInfo,
                                         double ScalingFactor)
{
    KRATOS_TRY;
    if (KratosComponents<Variable<double>>::Has(rVariable) == true)
    {
        CalculateNonHistoricalSensitivities(SensitivityVariables<double>{rVariable},
                                            rContainer, rResponseFunction,
                                            rProcessInfo, ScalingFactor);
    }
    else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(rVariable) == true)
    {
        CalculateNonHistoricalSensitivities(
            SensitivityVariables<array_1d<double, 3>>{rVariable}, rContainer,
            rResponseFunction, rProcessInfo, ScalingFactor);
    }
    else
    {
        KRATOS_ERROR << "Unsupported variable: " << rVariable << "." << std::endl;
    }
    KRATOS_CATCH("");
}

template <typename TContainer>
void CalculateNonHistoricalSensitivities(const std::vector<std::string>& rVariables,
                                         TContainer& rContainer,
                                         AdjointResponseFunction& rResponseFunction,
                                         const ProcessInfo& rProcessInfo,
                                         double ScalingFactor)
{
    KRATOS_TRY;
    for (auto& r_variable : rVariables)
    {
        CalculateNonHistoricalSensitivities(
            r_variable, rContainer, rResponseFunction, rProcessInfo, ScalingFactor);
    }
    KRATOS_CATCH("");
}

void SetNodalSolutionStepSensitivityVariableToZero(std::string const& rVariable,
                                                   ModelPart::NodesContainerType& rNodes)
{
    KRATOS_TRY;
    if (KratosComponents<Variable<double>>::Has(rVariable) == true)
    {
        auto sensitivity_variables = SensitivityVariables<double>{rVariable};
        VariableUtils().SetVariable(*sensitivity_variables.pOutputVariable,
                                    sensitivity_variables.pOutputVariable->Zero(), rNodes);
    }
    else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(rVariable) == true)
    {
        auto sensitivity_variables = SensitivityVariables<array_1d<double, 3>>{rVariable};
        VariableUtils().SetVariable(*sensitivity_variables.pOutputVariable,
                                    sensitivity_variables.pOutputVariable->Zero(), rNodes);
    }
    else
    {
        KRATOS_ERROR << "Unsupported variable: " << rVariable << "." << std::endl;
    }
    KRATOS_CATCH("");
}

void SetNodalSolutionStepSensitivityVariablesToZero(const std::vector<std::string>& rVariables,
                                                    ModelPart::NodesContainerType& rNodes)
{
    KRATOS_TRY;
    for (auto& r_variable : rVariables)
    {
        SetNodalSolutionStepSensitivityVariableToZero(r_variable, rNodes);
    }
    KRATOS_CATCH("");
}

template <typename TContainer>
void SetNonHistoricalSensitivityVariableToZero(std::string const& rVariable, TContainer& rElements)
{
    KRATOS_TRY;
    if (KratosComponents<Variable<double>>::Has(rVariable) == true)
    {
        auto sensitivity_variables = SensitivityVariables<double>{rVariable};
        VariableUtils().SetNonHistoricalVariable(
            *sensitivity_variables.pOutputVariable,
            sensitivity_variables.pOutputVariable->Zero(), rElements);
    }
    else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(rVariable) == true)
    {
        auto sensitivity_variables = SensitivityVariables<array_1d<double, 3>>{rVariable};
        VariableUtils().SetNonHistoricalVariable(
            *sensitivity_variables.pOutputVariable,
            sensitivity_variables.pOutputVariable->Zero(), rElements);
    }
    else
    {
        KRATOS_ERROR << "Unsupported variable: " << rVariable << "." << std::endl;
    }
    KRATOS_CATCH("");
}

template <typename TContainer>
void SetNonHistoricalSensitivityVariablesToZero(const std::vector<std::string>& rVariables,
                                                TContainer& rElements)
{
    KRATOS_TRY;
    for (auto& r_variable : rVariables)
    {
        SetNonHistoricalSensitivityVariableToZero(r_variable, rElements);
    }
    KRATOS_CATCH("");
}

void ReplaceDeprecatedNameIfExists(Parameters& rSettings,
                                   const std::string& rDeprecatedName,
                                   const std::string& rNewName)
{
    KRATOS_TRY;
    if (rSettings.Has(rDeprecatedName))
    {
        KRATOS_WARNING("SensitivityBuilder")
            << "Setting \"" << rDeprecatedName
            << "\" is deprecated and will be removed in the future." << std::endl;
        Parameters value = rSettings[rDeprecatedName].Clone();
        rSettings.RemoveValue(rDeprecatedName);
        rSettings.AddValue(rNewName, value);
    }
    KRATOS_CATCH("");
}

std::vector<std::string> ArrayOfStringToVector(Parameters Settings)
{
    KRATOS_DEBUG_ERROR_IF_NOT(Settings.IsArray())
        << "Settings must be an array." << std::endl;
    std::vector<std::string> result(Settings.size());
    for (std::size_t i = 0; i < result.size(); ++i)
    {
        result[i] = Settings.GetArrayItem(i).GetString();
    }
    return result;
}

} // namespace sensitivity_builder_cpp
} // namespace

namespace Kratos
{
SensitivityBuilder::SensitivityBuilder(Parameters Settings,
                                       ModelPart& rModelPart,
                                       AdjointResponseFunction::Pointer pResponseFunction)
    : mpModelPart(&rModelPart), mpResponseFunction(pResponseFunction)
{
    KRATOS_TRY;
    using sensitivity_builder_cpp::ReplaceDeprecatedNameIfExists;
    ReplaceDeprecatedNameIfExists(Settings, "nodal_sensitivity_variables",
                                  "nodal_solution_step_sensitivity_variables");
    ReplaceDeprecatedNameIfExists(Settings, "element_sensitivity_variables",
                                  "element_data_value_sensitivity_variables");
    ReplaceDeprecatedNameIfExists(Settings, "element_data_sensitivity_variables",
                                  "element_data_value_sensitivity_variables");
    ReplaceDeprecatedNameIfExists(Settings,
                                  "condition_data_sensitivity_variables",
                                  "condition_data_value_sensitivity_variables");
    ReplaceDeprecatedNameIfExists(Settings, "condition_sensitivity_variables",
                                  "condition_data_value_sensitivity_variables");
    Parameters default_settings(R"(
        {
            "sensitivity_model_part_name": "PLEASE_SPECIFY_SENSITIVITY_MODEL_PART",
            "nodal_solution_step_sensitivity_variables": [],
            "element_data_value_sensitivity_variables" : [],
            "condition_data_value_sensitivity_variables" : [],
            "build_mode": "static",
            "nodal_solution_step_sensitivity_calculation_is_thread_safe" : false
        })");
    Settings.ValidateAndAssignDefaults(default_settings);

    auto sensitivity_model_part_name =
        Settings["sensitivity_model_part_name"].GetString();
    if (sensitivity_model_part_name != "PLEASE_SPECIFY_SENSITIVITY_MODEL_PART")
    {
        mpSensitivityModelPart = mpModelPart->pGetSubModelPart(sensitivity_model_part_name);
    }
    else
    {
        mpSensitivityModelPart = mpModelPart;
    }
    using sensitivity_builder_cpp::ArrayOfStringToVector;
    mNodalSolutionStepSensitivityVariables = ArrayOfStringToVector(
        Settings["nodal_solution_step_sensitivity_variables"]);
    mElementDataValueSensitivityVariables = ArrayOfStringToVector(
        Settings["element_data_value_sensitivity_variables"]);
    mConditionDataValueSensitivityVariables = ArrayOfStringToVector(
        Settings["condition_data_value_sensitivity_variables"]);
    mBuildMode = Settings["build_mode"].GetString();
    mNodalSolutionStepSensitivityCalculationIsThreadSafe =
        Settings["nodal_solution_step_sensitivity_calculation_is_thread_safe"].GetBool();
    KRATOS_CATCH("");
}

void SensitivityBuilder::CalculateNodalSolutionStepSensitivities(
    const std::vector<std::string>& rVariables,
    ModelPart& rModelPart,
    AdjointResponseFunction& rResponseFunction,
    double ScalingFactor)
{
    KRATOS_TRY;
    using sensitivity_builder_cpp::CalculateNodalSolutionStepSensitivities;
    for (auto& r_variable : rVariables)
    {
        CalculateNodalSolutionStepSensitivities(
            r_variable, rModelPart, rResponseFunction, ScalingFactor);
    }
    KRATOS_CATCH("");
}

void SensitivityBuilder::CalculateNonHistoricalSensitivities(
    const std::vector<std::string>& rVariables,
    ModelPart::ElementsContainerType& rElements,
    AdjointResponseFunction& rResponseFunction,
    const ProcessInfo& rProcessInfo,
    double ScalingFactor)
{
    using sensitivity_builder_cpp::CalculateNonHistoricalSensitivities;
    CalculateNonHistoricalSensitivities(
        rVariables, rElements, rResponseFunction, rProcessInfo, ScalingFactor);
}

void SensitivityBuilder::CalculateNonHistoricalSensitivities(
    const std::vector<std::string>& rVariables,
    ModelPart::ConditionsContainerType& rConditions,
    AdjointResponseFunction& rResponseFunction,
    const ProcessInfo& rProcessInfo,
    double ScalingFactor)
{
    using sensitivity_builder_cpp::CalculateNonHistoricalSensitivities;
    CalculateNonHistoricalSensitivities(
        rVariables, rConditions, rResponseFunction, rProcessInfo, ScalingFactor);
}

void SensitivityBuilder::Initialize()
{
    KRATOS_TRY;
    Clear();
    VariableUtils().SetNonHistoricalVariable(UPDATE_SENSITIVITIES, true,
                                             mpSensitivityModelPart->Nodes());
    VariableUtils().SetNonHistoricalVariable(
        UPDATE_SENSITIVITIES, true, mpSensitivityModelPart->Elements());
    VariableUtils().SetNonHistoricalVariable(
        UPDATE_SENSITIVITIES, true, mpSensitivityModelPart->Conditions());
    KRATOS_CATCH("");
}

void SensitivityBuilder::UpdateSensitivities()
{
    KRATOS_TRY;
    double scaling_factor{};
    if (mBuildMode == "integrate")
    {
        // integrate in time
        scaling_factor = -mpModelPart->GetProcessInfo()[DELTA_TIME];
    }
    else if (mBuildMode == "sum")
    {
        scaling_factor = 1.0;
    }
    else if (mBuildMode == "static")
    {
        scaling_factor = 1.0;
        ClearSensitivities();
    }
    else
    {
        KRATOS_ERROR << "Unsupported \"build_mode\": " << mBuildMode << std::endl;
    }
    if (mNodalSolutionStepSensitivityCalculationIsThreadSafe)
    {
        CalculateNodalSolutionStepSensitivities(
            mNodalSolutionStepSensitivityVariables, *mpModelPart,
            *mpResponseFunction, scaling_factor);
    }
    else
    {
#ifdef _OPENMP
        const int max_threads = omp_get_max_threads();
        omp_set_num_threads(1);
#endif
        CalculateNodalSolutionStepSensitivities(
            mNodalSolutionStepSensitivityVariables, *mpModelPart,
            *mpResponseFunction, scaling_factor);
#ifdef _OPENMP
        omp_set_num_threads(max_threads);
#endif
    }
    CalculateNonHistoricalSensitivities(
        mElementDataValueSensitivityVariables, mpModelPart->Elements(),
        *mpResponseFunction, mpModelPart->GetProcessInfo(), scaling_factor);
    CalculateNonHistoricalSensitivities(
        mConditionDataValueSensitivityVariables, mpModelPart->Conditions(),
        *mpResponseFunction, mpModelPart->GetProcessInfo(), scaling_factor);
    KRATOS_CATCH("");
}

void SensitivityBuilder::Clear()
{
    KRATOS_TRY;

    ClearFlags();
    ClearSensitivities();

    KRATOS_CATCH("");
}

void SensitivityBuilder::ClearFlags()
{
    KRATOS_TRY;
    VariableUtils().SetNonHistoricalVariable(UPDATE_SENSITIVITIES, false,
                                             mpModelPart->Nodes());
    VariableUtils().SetNonHistoricalVariable(UPDATE_SENSITIVITIES, false,
                                             mpModelPart->Elements());
    VariableUtils().SetNonHistoricalVariable(UPDATE_SENSITIVITIES, false,
                                             mpModelPart->Conditions());
    KRATOS_CATCH("");
}

void SensitivityBuilder::ClearSensitivities()
{
    KRATOS_TRY;
    using sensitivity_builder_cpp::SetNodalSolutionStepSensitivityVariablesToZero;
    using sensitivity_builder_cpp::SetNonHistoricalSensitivityVariablesToZero;
    SetNodalSolutionStepSensitivityVariablesToZero(
        mNodalSolutionStepSensitivityVariables, mpModelPart->Nodes());
    SetNonHistoricalSensitivityVariablesToZero(
        mElementDataValueSensitivityVariables, mpModelPart->Elements());
    SetNonHistoricalSensitivityVariablesToZero(
        mConditionDataValueSensitivityVariables, mpModelPart->Conditions());
    KRATOS_CATCH("");
}

} // namespace Kratos
