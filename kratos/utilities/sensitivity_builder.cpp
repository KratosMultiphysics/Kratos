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
#include <algorithm>
#include <utility>

// External includes

// Project includes
#include "containers/global_pointers_vector.h"
#include "includes/kratos_parameters.h"
#include "includes/parallel_environment.h"
#include "includes/ublas_interface.h"
#include "input_output/logger.h"
#include "solving_strategies/schemes/sensitivity_builder_scheme.h"
#include "utilities/communication_coloring_utilities.h"
#include "utilities/openmp_utils.h"
#include "utilities/parallel_utilities.h"
#include "utilities/pointer_communicator.h"
#include "utilities/pointer_map_communicator.h"
#include "utilities/variable_utils.h"

// Include base h
#include "utilities/sensitivity_builder.h"

namespace sensitivity_builder_cpp // cotire guard
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

template <template <class T> class TFunctor, class TContainer, class... TArgs>
void ExecuteFunctorInContainer(
    const SensitivityBuilder::TSensitivityVariables& rVariables,
    TContainer& rContainer,
    TArgs&... rArgs)
{
    block_for_each(rContainer, [&](typename TContainer::value_type& rEntity) {
        ExecuteFunctor<TFunctor, typename TContainer::value_type, TArgs...>(
            rVariables, rEntity, rArgs...);
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

template<class TContainerType>
TContainerType& GetContainer(ModelPart& rModelPart);

template <>
ModelPart::ElementsContainerType& GetContainer<ModelPart::ElementsContainerType>(
    ModelPart& rModelPart)
{
    return rModelPart.Elements();
}

template <>
ModelPart::ConditionsContainerType& GetContainer<ModelPart::ConditionsContainerType>(
    ModelPart& rModelPart)
{
    return rModelPart.Conditions();
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

                // assign sensitivities to correct entities

                const IndexType data_size = sensitivities.size() / gp_vector.size();

                for (IndexType i = 0; i < gp_vector.size(); ++i) {
                    GetDataFromVector(sensitivities, i * data_size, data_size, data);
                    rProxy.Assign(gp_vector(i), data);
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

        GlobalPointerMapCommunicator<Node<3>, TDataType> pointer_map_communicator(
            rModelPart.GetCommunicator().GetDataCommunicator());

        auto apply_sensitivities_proxy = pointer_map_communicator.GetApplyProxy(
            [&](Node<3>& rNode, const TDataType& NewValue) {
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

void AddMatrixSubBlock(
    Matrix& rOutput,
    const Matrix& rInput,
    const int RowOffset,
    const int ColOffset)
{
    KRATOS_TRY

    const int rows = rInput.size1();
    const int cols = rInput.size2();

    KRATOS_DEBUG_ERROR_IF(static_cast<int>(rOutput.size1()) < RowOffset + rows)
        << "Output matrix number of rows is smaller than input matrix "
           "number of rows with offset. [ Output Matrix size = ( "
        << rOutput.size1() << ", " << rOutput.size2()
        << " ), Input matrix size = ( " << rows << ", " << cols
        << " ), Offsets = ( " << RowOffset << ", " << ColOffset << " ) ].";

    KRATOS_DEBUG_ERROR_IF(static_cast<int>(rOutput.size2()) < ColOffset + cols)
        << "Output matrix number of cols is smaller than input matrix "
           "number of cols with offset. [ Output Matrix size = ( "
        << rOutput.size1() << ", " << rOutput.size2()
        << " ), Input matrix size = ( " << rows << ", " << cols
        << " ), Offsets = ( " << RowOffset << ", " << ColOffset << " ) ].";

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            rOutput(i + RowOffset, j + ColOffset) += rInput(i, j);
        }
    }

    KRATOS_CATCH("");
}

void GetMatrixSubBlock(
    Matrix& rOutput,
    const Matrix& rInput,
    const int RowOffset,
    const int Rows,
    const int ColOffset,
    const int Cols)
{
    KRATOS_TRY

    KRATOS_DEBUG_ERROR_IF(static_cast<int>(rInput.size1()) < RowOffset + Rows)
        << "rInput matrix number of rows is smaller than number of rows "
           "with offset. [ Input Matrix size = ( "
        << rInput.size1() << ", " << rInput.size2() << " ), SubBlock size = ( "
        << Rows << ", " << Cols << " ), Offsets = ( " << RowOffset << ", "
        << ColOffset << " ) ].";

    KRATOS_DEBUG_ERROR_IF(static_cast<int>(rInput.size2()) < ColOffset + Cols)
        << "rInput matrix number of cols is smaller than number of cols "
           "with offset. [ Input Matrix size = ( "
        << rInput.size1() << ", " << rInput.size2() << " ), SubBlock size = ( "
        << Rows << ", " << Cols << " ), Offsets = ( " << RowOffset << ", "
        << ColOffset << " ) ].";

    KRATOS_DEBUG_ERROR_IF(static_cast<int>(rOutput.size1()) != Rows || static_cast<int>(rOutput.size2()) != Cols)
        << "Output matrix type mismatch. [ Output matrix size = ( "
        << rOutput.size1() << ", " << rOutput.size2() << " ), SubBlockSize = ( "
        << Rows << ", " << Cols << " ) ].";

    for (int i = 0; i < Rows; ++i) {
        for (int j = 0; j < Cols; ++j) {
            rOutput(i, j) = rInput(RowOffset + i, ColOffset + j);
        }
    }

    KRATOS_CATCH("");
}

void ComputeEntityGeometryNeighbourNodeMap(
    std::unordered_map<int, std::unordered_map<int, int>>& rDerivativeNodesMap,
    const std::unordered_map<int, std::vector<int>>& rNeighbourNodeIdsMap,
    const Geometry<ModelPart::NodeType>& rEntityGeometry,
    const Flags& rFlag,
    const bool CheckValue)
{
    KRATOS_TRY

    const int number_of_nodes = rEntityGeometry.PointsNumber();

    for (int i_base_node = 0; i_base_node < number_of_nodes; ++i_base_node) {
        const auto& r_base_node = rEntityGeometry[i_base_node];

        if (r_base_node.Is(rFlag) == CheckValue) {
            const int base_node_id = r_base_node.Id();

            std::unordered_map<int, int> derivative_node_map;
            for (int i_deriv_node = 0; i_deriv_node < number_of_nodes; ++i_deriv_node) {
                const int deriv_node_id = rEntityGeometry[i_deriv_node].Id();

                if (base_node_id == deriv_node_id) {
                    derivative_node_map[i_deriv_node] = 0;
                } else {
                    const auto p_itr = rNeighbourNodeIdsMap.find(base_node_id);
                    const auto& r_neighbour_node_indices = p_itr->second;

                    derivative_node_map[i_deriv_node] =
                        r_neighbour_node_indices.size() + 2;
                    for (int j = 0;
                         j < static_cast<int>(r_neighbour_node_indices.size()); ++j) {
                        if (r_neighbour_node_indices[j] == deriv_node_id) {
                            derivative_node_map[i_deriv_node] = j + 1;
                            break;
                        }
                    }

                    KRATOS_ERROR_IF(derivative_node_map[i_deriv_node] ==
                                    static_cast<int>(r_neighbour_node_indices.size() + 2))
                        << "Derivative node id " << deriv_node_id
                        << " not found in neighbour nodes list in node with id "
                        << base_node_id << ".";
                }
            }

            rDerivativeNodesMap[i_base_node] = derivative_node_map;
        }
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

    mpSensitivityBuilderScheme->Initialize(
        *mpModelPart, *mpSensitivityModelPart, *mpResponseFunction);

    KRATOS_CATCH("");
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

template <class TContainerType>
void SensitivityBuilder::AssignEntityDerivativesToNodes(
    ModelPart& rModelPart,
    const int DerivativeDimension,
    const Variable<Matrix>& rDerivativeVariable,
    const std::unordered_map<int, std::vector<int>>& rNeighbourNodeIdsMap,
    const double Weight,
    const Flags& rFlag,
    const bool CheckValue)
{
    KRATOS_TRY

    using sensitivity_builder_cpp::AddMatrixSubBlock;
    using sensitivity_builder_cpp::GetContainer;
    using sensitivity_builder_cpp::GetMatrixSubBlock;
    using sensitivity_builder_cpp::ComputeEntityGeometryNeighbourNodeMap;


    auto& entity_container = GetContainer<TContainerType>(rModelPart);

    if (entity_container.size() != 0) {
        auto& r_nodes = rModelPart.Nodes();

        const int value_dimension =
            entity_container.begin()->GetValue(rDerivativeVariable).size2();

        KRATOS_ERROR_IF(value_dimension == 0)
            << "Column dimension (representing dimensionality of the value "
               "where derivatives are calculated) of the matrix values are "
               "zero in "
            << rModelPart.Name() << ". Please assign proper matrix values for "
            << rDerivativeVariable.Name() << ".";

        VariableUtils().SetFlag(VISITED, false, r_nodes);

        // identify nodes where neighbours are required
        block_for_each(entity_container, [&](typename TContainerType::value_type& rEntity) {
            if (rEntity.Is(rFlag) == CheckValue) {

                KRATOS_ERROR_IF(!rEntity.Has(rDerivativeVariable))
                    << rDerivativeVariable.Name() << " not found in data value container of "
                    << rEntity.Info() << ".";

                const Matrix& r_value = rEntity.GetValue(rDerivativeVariable);

                KRATOS_ERROR_IF(value_dimension != static_cast<int>(r_value.size2()))
                    << rDerivativeVariable.Name()
                    << " matrix value second dimension is not consistent at "
                    << rEntity.Info() << " [ required dimension size: " << value_dimension
                    << ", obtained matrix: " << r_value << " ].\n";

                auto& r_geometry = rEntity.GetGeometry();
                for (auto& r_node : r_geometry) {
                    if (r_node.Is(rFlag) == CheckValue) {
                        r_node.SetLock();
                        r_node.Set(VISITED, true);
                        r_node.UnSetLock();
                    }
                }
            }
        });

        rModelPart.GetCommunicator().SynchronizeOrNodalFlags(VISITED);

        // resizing matrices
        block_for_each(r_nodes, [&](ModelPart::NodeType& rNode) {
            if (rNode.Is(VISITED)) {
                const int node_id = rNode.Id();

                const auto p_itr = rNeighbourNodeIdsMap.find(node_id);
                KRATOS_ERROR_IF(p_itr == rNeighbourNodeIdsMap.end())
                    << node_id << " is not found in nodal neighbours map.";

                const int number_of_neighbour_nodes = p_itr->second.size();

                // neighbour nodes only contain neighbours, not self. But derivatives are there for neighbour
                // nodes and self node. So we reserve space for number_of_neighbour_nodes + 1 nodes
                // first block always represent self node, rest in the order of rNeighbourNodeIdsMap vector
                rNode.SetValue(rDerivativeVariable,
                               Matrix((number_of_neighbour_nodes + 1) * DerivativeDimension,
                                      value_dimension, 0.0));
            }
            else {
                // initializing these unused node matrices are required to do non-historical assembly
                // otherwise, this method will fail with TContainerType = ModelPart::ConditionsContainerType
                // at TransferDistributedValues in MPICommunicator::AssembleDynamicMatrixValues
                rNode.SetValue(rDerivativeVariable, Matrix(1, 1, 0.0));
            }
        });

        block_for_each(entity_container, [&](typename TContainerType::value_type& rEntity) {
            if (rEntity.Is(rFlag) == CheckValue) {
                auto& r_geometry = rEntity.GetGeometry();
                const int number_of_nodes = r_geometry.PointsNumber();

                std::unordered_map<int, std::unordered_map<int, int>> derivative_nodes_map;

                // calculate the node mapping
                ComputeEntityGeometryNeighbourNodeMap(
                    derivative_nodes_map, rNeighbourNodeIdsMap, r_geometry, rFlag, CheckValue);

                const Matrix& r_entity_derivatives =
                    rEntity.GetValue(rDerivativeVariable) * Weight;

                // move this variable also to TLS storage
                Matrix nodal_derivative(DerivativeDimension, value_dimension);

                // placing derivatives correctly
                for (int i_base_node = 0; i_base_node < number_of_nodes; ++i_base_node) {
                    auto& r_base_node = r_geometry[i_base_node];

                    if (r_base_node.Is(rFlag) == CheckValue) {
                        const auto& r_derivative_nodes_map =
                            derivative_nodes_map.find(i_base_node)->second;

                        for (int i_deriv_node = 0; i_deriv_node < number_of_nodes; ++i_deriv_node) {
                            GetMatrixSubBlock(nodal_derivative, r_entity_derivatives,
                                              i_deriv_node * DerivativeDimension,
                                              DerivativeDimension, 0, value_dimension);

                            r_base_node.SetLock();
                            AddMatrixSubBlock(
                                r_base_node.GetValue(rDerivativeVariable), nodal_derivative,
                                r_derivative_nodes_map.find(i_deriv_node)->second * DerivativeDimension,
                                0);
                            r_base_node.UnSetLock();
                        }
                    }
                }
            }
        });

        rModelPart.GetCommunicator().AssembleNonHistoricalData(rDerivativeVariable);
    }

    KRATOS_CATCH("");
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

// template instantiations
template KRATOS_API(KRATOS_CORE) void SensitivityBuilder::AssignEntityDerivativesToNodes<ModelPart::ElementsContainerType>
    (
        ModelPart&,
        const int,
        const Variable<Matrix>&,
        const std::unordered_map<int, std::vector<int>>&,
        const double,
        const Flags&,
        const bool
    );

template KRATOS_API(KRATOS_CORE) void SensitivityBuilder::AssignEntityDerivativesToNodes<ModelPart::ConditionsContainerType>
    (
        ModelPart&,
        const int,
        const Variable<Matrix>&,
        const std::unordered_map<int, std::vector<int>>&,
        const double,
        const Flags&,
        const bool
    );


} // namespace Kratos
