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

// System includes

// External includes

// Project includes
#include "containers/global_pointers_vector.h"
#include "includes/kratos_parameters.h"
#include "includes/parallel_environment.h"
#include "includes/ublas_interface.h"
#include "input_output/logger.h"
#include "solving_strategies/schemes/sensitivity_builder_scheme.h"
#include "utilities/openmp_utils.h"
#include "utilities/parallel_utilities.h"
#include "utilities/pointer_map_communicator.h"
#include "utilities/variable_utils.h"

// Include base h
#include "utilities/sensitivity_builder.h"

namespace sensitivity_builder_cpp // unity build guard
{
using namespace Kratos;

template <template <class T> class TFunctor, class... TArgs>
void ExecuteFunctor(
    const SensitivityBuilder::TSensitivityVariables& rVariables,
    TArgs&... rArgs)
{
    for (const auto& r_variable : std::get<0>(rVariables)) {
        TFunctor<double>()(r_variable, rArgs...);
    }
    for (const auto& r_variable : std::get<1>(rVariables)) {
        TFunctor<array_1d<double, 3>>()(r_variable, rArgs...);
    }
}

template <template <class T> class TFunctor, class TContainer>
void ExecuteFunctorInContainer(
    const SensitivityBuilder::TSensitivityVariables& rVariables,
    TContainer& rContainer)
{
    block_for_each(rContainer, [&](typename TContainer::value_type& rEntity) {
        ExecuteFunctor<TFunctor, typename TContainer::value_type>(
            rVariables, rEntity);
    });
}

template <class TDataType>
class SetNonHistoricalValueToZeroFunctor
{
public:
    template <class TEntityType>
    void operator()(
        const SensitivityBuilder::SensitivityVariables<TDataType>& rVariable,
        TEntityType& rEntity)
    {
        rEntity.SetValue(*rVariable.pOutputVariable, rVariable.pOutputVariable->Zero());
    }
};

template <class TDataType>
class SetHistoricalValueToZeroFunctor
{
public:
    void operator()(
        const SensitivityBuilder::SensitivityVariables<TDataType>& rVariable,
        ModelPart::NodeType& rNode)
    {
        rNode.FastGetSolutionStepValue(*rVariable.pOutputVariable) =
            rVariable.pOutputVariable->Zero();
    }
};

template<class TDerivativeEntityType>
bool HasSensitivityContributions(
    const Geometry<ModelPart::NodeType>& rGeometry)
{
    return rGeometry.GetValue(UPDATE_SENSITIVITIES);
}

template<>
bool HasSensitivityContributions<ModelPart::NodeType>(
    const Geometry<ModelPart::NodeType>& rGeometry)
{
    for (auto& r_node : rGeometry) {
        if (r_node.GetValue(UPDATE_SENSITIVITIES)) {
            return true;
        }
    }
    return false;
}

template <class TDataType>
void GetDataFromVector(
    const Vector& rValues,
    const IndexType Position,
    const IndexType OutputSize,
    TDataType& rOutput);

template<>
void GetDataFromVector(
    const Vector& rValues,
    const IndexType Position,
    const IndexType OutputSize,
    double& rOutput)
{
    KRATOS_TRY

    KRATOS_DEBUG_ERROR_IF(rValues.size() <= Position)
        << "rValues vector size mismatch in double type [ rValues.size() = " << rValues.size()
        << ", required size > " << Position << " ].\n";

    KRATOS_DEBUG_ERROR_IF(OutputSize != 1)
        << "OutputSize should be 1 for double values. [ OutputSize = " << OutputSize
        << " ].\n";

    rOutput = rValues[Position];

    KRATOS_CATCH("");
}

template<>
void GetDataFromVector(
    const Vector& rValues,
    const IndexType Position,
    const IndexType OutputSize,
    array_1d<double, 3>& rOutput)
{
    KRATOS_TRY

    KRATOS_DEBUG_ERROR_IF(rValues.size() < Position + OutputSize)
        << "rValues vector size mismatch in array_1d<double, 3> type [ "
           "rValues.size() = "
        << rValues.size() << ", required size >= " << (Position + OutputSize)
        << ", requested position = " << Position
        << ", requested OutputSize = " << OutputSize << " ].\n";

    KRATOS_DEBUG_ERROR_IF(OutputSize > 3)
        << "OutputSize should be less than or equal to 3 [ OutputSize = " << OutputSize
        << " ].\n";

    rOutput.clear();

    for (IndexType i = 0; i < OutputSize; ++i) {
        rOutput[i] = rValues[Position + i];
    }

    KRATOS_CATCH("");
}

template<class TContainerType, class TDataType, class TProxyType>
void AssembleContainerContributions(
    TContainerType& rContainer,
    AdjointResponseFunction& rResponseFunction,
    SensitivityBuilderScheme& rSensitivityBuilderScheme,
    TProxyType& rProxy,
    const Variable<TDataType>& rDesignVariable,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY

    using TDerivativeEntityType = typename TProxyType::GlobalPointerType::element_type;

    using TThreadLocalStorageType = std::tuple<Vector, GlobalPointersVector<TDerivativeEntityType>, TDataType>;

    block_for_each(
        rContainer, TThreadLocalStorageType(),
        [&](typename TContainerType::value_type& rEntity, TThreadLocalStorageType& rTLS) {
            auto& r_geometry = rEntity.GetGeometry();

            auto& sensitivities = std::get<0>(rTLS);
            auto& gp_vector = std::get<1>(rTLS);
            auto& data = std::get<2>(rTLS);

            if (HasSensitivityContributions<TDerivativeEntityType>(r_geometry)) {
                rSensitivityBuilderScheme.CalculateSensitivity(
                    rEntity, rResponseFunction, sensitivities, gp_vector,
                    rDesignVariable, rProcessInfo);

                KRATOS_DEBUG_ERROR_IF(sensitivities.size() != gp_vector.size() &&
                                      sensitivities.size() == 0)
                    << "gp_vector and sensitivities size mismatch [ "
                       "gp_vector.size() = "
                    << gp_vector.size()
                    << ", sensitivities.size() = " << sensitivities.size() << " ].\n";

                // assign sensitivities to correct entities
                if (gp_vector.size() != 0) {
                    const IndexType data_size = sensitivities.size() / gp_vector.size();

                    for (IndexType i = 0; i < gp_vector.size(); ++i) {
                        GetDataFromVector(sensitivities, i * data_size, data_size, data);
                        rProxy.Assign(gp_vector(i), data);
                    }
                }
            }
        });

    KRATOS_CATCH("");
}

template <class TDataType>
class CalculateNodalSolutionStepSensitivityFunctor
{
public:
    void operator()(
        const SensitivityBuilder::SensitivityVariables<TDataType>& rVariable,
        ModelPart& rModelPart,
        AdjointResponseFunction& rResponseFunction,
        SensitivityBuilderScheme& rSensitivityBuilderScheme,
        const double& ScalingFactor)
    {
        KRATOS_TRY

        GlobalPointerMapCommunicator<Node, TDataType> pointer_map_communicator(
            rModelPart.GetCommunicator().GetDataCommunicator());

        auto apply_sensitivities_proxy = pointer_map_communicator.GetApplyProxy(
            [&](Node& rNode, const TDataType& NewValue) {
                if (rNode.Has(UPDATE_SENSITIVITIES) && rNode.GetValue(UPDATE_SENSITIVITIES)) {
                    rNode.SetLock();
                    rNode.FastGetSolutionStepValue(*(rVariable.pOutputVariable)) += NewValue * ScalingFactor;
                    rNode.UnSetLock();
                }
            });

        AssembleContainerContributions(
            rModelPart.Elements(), rResponseFunction, rSensitivityBuilderScheme,
            apply_sensitivities_proxy, *rVariable.pDesignVariable,
            rModelPart.GetProcessInfo());

        AssembleContainerContributions(
            rModelPart.Conditions(), rResponseFunction, rSensitivityBuilderScheme,
            apply_sensitivities_proxy, *rVariable.pDesignVariable,
            rModelPart.GetProcessInfo());

        apply_sensitivities_proxy.SendAndApplyRemotely();

        // synchronize to populate ghost mesh properly
        rModelPart.GetCommunicator().SynchronizeVariable(*(rVariable.pOutputVariable));

        KRATOS_CATCH("");
    }
};

template <class TDataType>
class CalculateNonHistoricalSensitivitiesFunctor
{
public:
    template<class TContainerType>
    void operator()(
        const SensitivityBuilder::SensitivityVariables<TDataType>& rVariable,
        TContainerType& rContainer,
        AdjointResponseFunction& rResponseFunction,
        SensitivityBuilderScheme& rSensitivityBuilderScheme,
        const DataCommunicator& rDataCommunicator,
        const ProcessInfo& rProcessInfo,
        const double& ScalingFactor)
    {
        KRATOS_TRY

        GlobalPointerMapCommunicator<typename TContainerType::value_type, TDataType> pointer_map_communicator(
            rDataCommunicator);

        auto apply_sensitivities_proxy = pointer_map_communicator.GetApplyProxy(
            [&](typename TContainerType::value_type& rEntity, const TDataType& NewValue) {
                if (rEntity.Has(UPDATE_SENSITIVITIES) && rEntity.GetValue(UPDATE_SENSITIVITIES)) {
                    rEntity.GetValue(*(rVariable.pOutputVariable)) += NewValue * ScalingFactor;
                }
            });

        AssembleContainerContributions(
            rContainer, rResponseFunction, rSensitivityBuilderScheme,
            apply_sensitivities_proxy, *rVariable.pDesignVariable,
            rProcessInfo);

        apply_sensitivities_proxy.SendAndApplyRemotely();

        KRATOS_CATCH("");
    }
};

void ReplaceDeprecatedNameIfExists(
    Parameters& rSettings,
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

} // namespace sensitivity_builder_cpp

namespace Kratos
{
SensitivityBuilder::SensitivityBuilder(
    Parameters Settings,
    ModelPart& rModelPart,
    AdjointResponseFunction::Pointer pResponseFunction)
    : SensitivityBuilder(
          Settings,
          rModelPart,
          pResponseFunction,
          Kratos::make_unique<SensitivityBuilderScheme>())
{
}

SensitivityBuilder::SensitivityBuilder(
    Parameters Settings,
    ModelPart& rModelPart,
    AdjointResponseFunction::Pointer pResponseFunction,
    SensitivityBuilderScheme::Pointer pSensitivityBuilderScheme)
    : mpModelPart(&rModelPart),
      mpResponseFunction(pResponseFunction),
      mpSensitivityBuilderScheme(pSensitivityBuilderScheme)
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
    if (sensitivity_model_part_name !=
        "PLEASE_SPECIFY_SENSITIVITY_MODEL_PART") {
        mpSensitivityModelPart = mpModelPart->pGetSubModelPart(sensitivity_model_part_name);
    } else {
        mpSensitivityModelPart = mpModelPart;
    }

    mBuildMode = Settings["build_mode"].GetString();
    mNodalSolutionStepSensitivityCalculationIsThreadSafe =
        Settings["nodal_solution_step_sensitivity_calculation_is_thread_safe"].GetBool();

    mNodalSolutionStepSensitivityVariablesList = SensitivityBuilder::GetVariableLists(
        Settings["nodal_solution_step_sensitivity_variables"].GetStringArray());
    mElementDataValueSensitivityVariablesList = SensitivityBuilder::GetVariableLists(
        Settings["element_data_value_sensitivity_variables"].GetStringArray());
    mConditionDataValueSensitivityVariablesList = SensitivityBuilder::GetVariableLists(
        Settings["condition_data_value_sensitivity_variables"].GetStringArray());

    KRATOS_CATCH("");
}

void SensitivityBuilder::SetResponseFunction(AdjointResponseFunction::Pointer pResponseFunction)
{
    mpResponseFunction = pResponseFunction;
}

void SensitivityBuilder::Initialize()
{
    KRATOS_TRY;

    ClearFlags();
    ClearSensitivities();
    VariableUtils().SetNonHistoricalVariable(UPDATE_SENSITIVITIES, true,
                                             mpSensitivityModelPart->Nodes());
    VariableUtils().SetNonHistoricalVariable(
        UPDATE_SENSITIVITIES, true, mpSensitivityModelPart->Elements());
    VariableUtils().SetNonHistoricalVariable(
        UPDATE_SENSITIVITIES, true, mpSensitivityModelPart->Conditions());

    mpSensitivityBuilderScheme->Initialize(
        *mpModelPart, *mpSensitivityModelPart, *mpResponseFunction);

    KRATOS_CATCH("");
}

void SensitivityBuilder::InitializeSolutionStep()
{
    mpSensitivityBuilderScheme->InitializeSolutionStep(
        *mpModelPart, *mpSensitivityModelPart, *mpResponseFunction);
}

void SensitivityBuilder::UpdateSensitivities()
{
    KRATOS_TRY;
    double scaling_factor{};
    if (mBuildMode == "integrate") {
        // integrate in time
        scaling_factor = -mpModelPart->GetProcessInfo()[DELTA_TIME];
    } else if (mBuildMode == "sum") {
        scaling_factor = 1.0;
    } else if (mBuildMode == "static") {
        scaling_factor = 1.0;
        ClearSensitivities();
    } else {
        KRATOS_ERROR << "Unsupported \"build_mode\": " << mBuildMode << std::endl;
    }
    if (mNodalSolutionStepSensitivityCalculationIsThreadSafe) {
        CalculateNodalSolutionStepSensitivities(
            mNodalSolutionStepSensitivityVariablesList, *mpModelPart,
            *mpResponseFunction, *mpSensitivityBuilderScheme, scaling_factor);
    } else {
#ifdef _OPENMP
        const int max_threads = omp_get_max_threads();
        omp_set_num_threads(1);
#endif
        CalculateNodalSolutionStepSensitivities(
            mNodalSolutionStepSensitivityVariablesList, *mpModelPart,
            *mpResponseFunction, *mpSensitivityBuilderScheme, scaling_factor);
#ifdef _OPENMP
        omp_set_num_threads(max_threads);
#endif
    }
    CalculateNonHistoricalSensitivities(
        mElementDataValueSensitivityVariablesList, mpModelPart->Elements(),
        *mpResponseFunction, *mpSensitivityBuilderScheme, mpModelPart->GetProcessInfo(), scaling_factor);
    CalculateNonHistoricalSensitivities(
        mConditionDataValueSensitivityVariablesList, mpModelPart->Conditions(),
        *mpResponseFunction, *mpSensitivityBuilderScheme, mpModelPart->GetProcessInfo(), scaling_factor);

    mpSensitivityBuilderScheme->Update(
        *mpModelPart, *mpSensitivityModelPart, *mpResponseFunction);

    KRATOS_CATCH("");
}

void SensitivityBuilder::FinalizeSolutionStep()
{
    mpSensitivityBuilderScheme->FinalizeSolutionStep(
        *mpModelPart, *mpSensitivityModelPart, *mpResponseFunction);
}

void SensitivityBuilder::Finalize()
{
    mpSensitivityBuilderScheme->Finalize(
        *mpModelPart, *mpSensitivityModelPart, *mpResponseFunction);
}

void SensitivityBuilder::Clear()
{
    KRATOS_TRY;

    ClearFlags();
    ClearSensitivities();
    mpSensitivityBuilderScheme->Clear();

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

    using namespace sensitivity_builder_cpp;

    ExecuteFunctorInContainer<SetNonHistoricalValueToZeroFunctor>(
        mElementDataValueSensitivityVariablesList, mpModelPart->Elements());
    ExecuteFunctorInContainer<SetNonHistoricalValueToZeroFunctor>(
        mConditionDataValueSensitivityVariablesList, mpModelPart->Conditions());
    ExecuteFunctorInContainer<SetHistoricalValueToZeroFunctor>(
        mNodalSolutionStepSensitivityVariablesList, mpModelPart->Nodes());

    KRATOS_CATCH("");
}

void SensitivityBuilder::CalculateNodalSolutionStepSensitivities(
    const std::vector<std::string>& rVariables,
    ModelPart& rModelPart,
    AdjointResponseFunction& rResponseFunction,
    double ScalingFactor)
{
    KRATOS_TRY;

    using namespace sensitivity_builder_cpp;

    const auto& r_variables_list = GetVariableLists(rVariables);
    SensitivityBuilderScheme scheme;

    scheme.Initialize(rModelPart, rModelPart, rResponseFunction);

    CalculateNodalSolutionStepSensitivities(
        r_variables_list, rModelPart, rResponseFunction, scheme, ScalingFactor);

    KRATOS_CATCH("");
}

void SensitivityBuilder::CalculateNodalSolutionStepSensitivities(
    const TSensitivityVariables& rVariables,
    ModelPart& rModelPart,
    AdjointResponseFunction& rResponseFunction,
    SensitivityBuilderScheme& rSensitivityBuilderScheme,
    double ScalingFactor)
{
    KRATOS_TRY;

    using namespace sensitivity_builder_cpp;

    ExecuteFunctor<CalculateNodalSolutionStepSensitivityFunctor>(
        rVariables, rModelPart, rResponseFunction, rSensitivityBuilderScheme, ScalingFactor);

    KRATOS_CATCH("");
}

void SensitivityBuilder::CalculateNonHistoricalSensitivities(
    const std::vector<std::string>& rVariables,
    ModelPart::ElementsContainerType& rElements,
    AdjointResponseFunction& rResponseFunction,
    const ProcessInfo& rProcessInfo,
    double ScalingFactor)
{
    using namespace sensitivity_builder_cpp;

    const auto& r_variables_list = GetVariableLists(rVariables);
    SensitivityBuilderScheme scheme;

    CalculateNonHistoricalSensitivities(r_variables_list, rElements, rResponseFunction,
                                        scheme, rProcessInfo, ScalingFactor);
}

void SensitivityBuilder::CalculateNonHistoricalSensitivities(
    const TSensitivityVariables& rVariables,
    ModelPart::ElementsContainerType& rElements,
    AdjointResponseFunction& rResponseFunction,
    SensitivityBuilderScheme& rSensitivityBuilderScheme,
    const ProcessInfo& rProcessInfo,
    double ScalingFactor)
{
    using namespace sensitivity_builder_cpp;

    ExecuteFunctor<CalculateNonHistoricalSensitivitiesFunctor>(
        rVariables, rElements, rResponseFunction, rSensitivityBuilderScheme,
        ParallelEnvironment::GetDefaultDataCommunicator(), rProcessInfo, ScalingFactor);
}

void SensitivityBuilder::CalculateNonHistoricalSensitivities(
    const std::vector<std::string>& rVariables,
    ModelPart::ConditionsContainerType& rConditions,
    AdjointResponseFunction& rResponseFunction,
    const ProcessInfo& rProcessInfo,
    double ScalingFactor)
{
    using namespace sensitivity_builder_cpp;

    const auto& r_variables_list = GetVariableLists(rVariables);
    SensitivityBuilderScheme scheme;

    CalculateNonHistoricalSensitivities(r_variables_list, rConditions, rResponseFunction,
                                        scheme, rProcessInfo, ScalingFactor);
}

void SensitivityBuilder::CalculateNonHistoricalSensitivities(
    const TSensitivityVariables& rVariables,
    ModelPart::ConditionsContainerType& rConditions,
    AdjointResponseFunction& rResponseFunction,
    SensitivityBuilderScheme& rSensitivityBuilderScheme,
    const ProcessInfo& rProcessInfo,
    double ScalingFactor)
{
    using namespace sensitivity_builder_cpp;

    ExecuteFunctor<CalculateNonHistoricalSensitivitiesFunctor>(
        rVariables, rConditions, rResponseFunction, rSensitivityBuilderScheme,
        ParallelEnvironment::GetDefaultDataCommunicator(), rProcessInfo, ScalingFactor);
}

SensitivityBuilder::TSensitivityVariables SensitivityBuilder::GetVariableLists(
        const std::vector<std::string>& rVariableNames)
{
    KRATOS_TRY;

    THomogeneousSensitivityVariables<double> double_variables;
    THomogeneousSensitivityVariables<array_1d<double, 3>> array_3d_variables;

    for (const auto& r_variable : rVariableNames) {
        if (KratosComponents<Variable<double>>::Has(r_variable)) {
            double_variables.push_back(SensitivityVariables<double>{r_variable});
        } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(r_variable)) {
            array_3d_variables.push_back(SensitivityVariables<array_1d<double, 3>>{r_variable});
        } else {
            KRATOS_ERROR << "Unsupported variable: " << r_variable << "." << std::endl;
        }
    }

    return std::make_tuple(double_variables, array_3d_variables);

    KRATOS_CATCH("");
}

} // namespace Kratos
