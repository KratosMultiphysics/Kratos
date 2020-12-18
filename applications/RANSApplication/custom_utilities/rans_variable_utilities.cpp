//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <cmath>
#include <limits>

// External includes

// Project includes
#include "includes/define.h"
#include "utilities/variable_utils.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// Application includes
#include "rans_application_variables.h"

// Include base h
#include "rans_variable_utilities.h"

namespace Kratos
{
namespace RansVariableUtilities
{
std::tuple<unsigned int, unsigned int> ClipScalarVariable(
    const double MinimumValue,
    const double MaximumValue,
    const Variable<double>& rVariable,
    ModelPart& rModelPart)
{
    KRATOS_TRY

    auto& r_communicator = rModelPart.GetCommunicator();
    auto& r_nodes = r_communicator.LocalMesh().Nodes();

    unsigned int number_of_nodes_below_minimum, number_of_nodes_above_maximum;
    std::tie(number_of_nodes_below_minimum, number_of_nodes_above_maximum) =
        block_for_each<CombinedReduction<SumReduction<unsigned int>, SumReduction<unsigned int>>>(
            r_nodes, [&](ModelPart::NodeType& rNode) -> std::tuple<unsigned int, unsigned int> {
                double& r_value = rNode.FastGetSolutionStepValue(rVariable);

                if (r_value < MinimumValue) {
                    r_value = MinimumValue;
                    return std::make_tuple<unsigned int, unsigned int>(1, 0);
                } else if (r_value > MaximumValue) {
                    r_value = MaximumValue;
                    return std::make_tuple<unsigned int, unsigned int>(0, 1);
                }

                return std::make_tuple<unsigned int, unsigned int>(0, 0);
            });

    r_communicator.SynchronizeVariable(rVariable);

    // Stores followings
    // index - 0 : number_of_nodes_below_minimum
    // index - 1 : number_of_nodes_above_maximum
    std::vector<unsigned int> nodes_count = {number_of_nodes_below_minimum,
                                             number_of_nodes_above_maximum};
    const std::vector<unsigned int>& total_nodes_count =
        r_communicator.GetDataCommunicator().SumAll(nodes_count);

    return std::tuple<unsigned int, unsigned int>(total_nodes_count[0], total_nodes_count[1]);

    KRATOS_CATCH("")
}

double GetMinimumScalarValue(
    const ModelPart& rModelPart,
    const Variable<double>& rVariable)
{
    KRATOS_TRY

    const auto& r_communicator = rModelPart.GetCommunicator();
    const auto& r_nodes = r_communicator.LocalMesh().Nodes();

    const int number_of_nodes = r_nodes.size();
    const double min_value =
        IndexPartition<int>(number_of_nodes).for_each<MinReduction<double>>([&](const int i) -> double {
            return (r_nodes.begin() + i)->FastGetSolutionStepValue(rVariable);
        });

    return r_communicator.GetDataCommunicator().MinAll(min_value);

    KRATOS_CATCH("");
}

double GetMaximumScalarValue(
    const ModelPart& rModelPart,
    const Variable<double>& rVariable)
{
    KRATOS_TRY

    const auto& r_communicator = rModelPart.GetCommunicator();
    const auto& r_nodes = r_communicator.LocalMesh().Nodes();

    const int number_of_nodes = r_nodes.size();
    const double max_value =
        IndexPartition<int>(number_of_nodes).for_each<MaxReduction<double>>([&](const int i) -> double {
            return (r_nodes.begin() + i)->FastGetSolutionStepValue(rVariable);
        });

    return r_communicator.GetDataCommunicator().MaxAll(max_value);

    KRATOS_CATCH("");
}

void GetNodalVariablesVector(
    Vector& rValues,
    const ModelPart::NodesContainerType& rNodes,
    const Variable<double>& rVariable)
{
    const int number_of_nodes = rNodes.size();

    if (static_cast<int>(rValues.size()) != number_of_nodes) {
        rValues.resize(number_of_nodes);
    }

    IndexPartition<int>(number_of_nodes).for_each([&](const int i_node) {
        rValues[i_node] = (rNodes.begin() + i_node)->FastGetSolutionStepValue(rVariable);
    });
}

void GetNodalArray(
    Vector& rNodalValues,
    const Element& rElement,
    const Variable<double>& rVariable)
{
    const auto& r_geometry = rElement.GetGeometry();
    std::size_t number_of_nodes = r_geometry.PointsNumber();

    if (rNodalValues.size() != number_of_nodes) {
        rNodalValues.resize(number_of_nodes);
    }

    for (std::size_t i_node = 0; i_node < number_of_nodes; ++i_node) {
        rNodalValues[i_node] = r_geometry[i_node].FastGetSolutionStepValue(rVariable);
    }
}

void SetNodalVariables(
    ModelPart::NodesContainerType& rNodes,
    const Vector& rValues,
    const Variable<double>& rVariable)
{
    const int number_of_nodes = rNodes.size();

    KRATOS_ERROR_IF(static_cast<int>(rValues.size()) != number_of_nodes)
        << "rValues vector size mismatch with rNodes size in "
           "SetNodalVariables. [ rValues.size = "
        << rValues.size() << ", rNodes.size = " << rNodes.size() << " ]\n";

    IndexPartition<int>(number_of_nodes).for_each([&](const int i_node) {
        (rNodes.begin() + i_node)->FastGetSolutionStepValue(rVariable) = rValues[i_node];
    });
}

void CopyNodalSolutionStepVariablesList(
    ModelPart& rOriginModelPart,
    ModelPart& rDestinationModelPart)
{
    KRATOS_TRY

    rDestinationModelPart.GetNodalSolutionStepVariablesList() =
        rOriginModelPart.GetNodalSolutionStepVariablesList();

    KRATOS_CATCH("");
}

template <typename TDataType>
void AssignConditionVariableValuesToNodes(
    ModelPart& rModelPart,
    const Variable<TDataType>& rVariable,
    const Flags& rFlag,
    const bool FlagValue)
{
    auto& r_nodes = rModelPart.Nodes();
    VariableUtils().SetHistoricalVariableToZero(rVariable, r_nodes);

    block_for_each(rModelPart.Conditions(), [&](ModelPart::ConditionType& rCondition) {
        if (rCondition.Is(rFlag) == FlagValue) {
            const int number_of_nodes = rCondition.GetGeometry().PointsNumber();
            const auto& r_value = rCondition.GetValue(rVariable);
            for (int i_node = 0; i_node < number_of_nodes; ++i_node) {
                auto& r_node = rCondition.GetGeometry()[i_node];
                r_node.SetLock();
                r_node.FastGetSolutionStepValue(rVariable) +=
                    r_value * (1.0 / static_cast<double>(number_of_nodes));
                r_node.UnSetLock();
            }
        }
    });

    rModelPart.GetCommunicator().AssembleCurrentData(rVariable);
}

void AddAnalysisStep(
    ModelPart& rModelPart,
    const std::string& rStepName)
{
    auto& r_process_info = rModelPart.GetProcessInfo();
    if (!r_process_info.Has(ANALYSIS_STEPS))
    {
        r_process_info.SetValue(ANALYSIS_STEPS, std::vector<std::string>());
    }
    r_process_info[ANALYSIS_STEPS].push_back(rStepName);
}

bool IsAnalysisStepCompleted(
    const ModelPart& rModelPart,
    const std::string& rStepName)
{
    const auto& r_process_info = rModelPart.GetProcessInfo();
    if (r_process_info.Has(ANALYSIS_STEPS))
    {
        const std::vector<std::string>& r_steps = r_process_info[ANALYSIS_STEPS];
        return (std::find(r_steps.begin(), r_steps.end(), rStepName) != r_steps.end());
    }
    else
    {
        return false;
    }
}

void AssignBoundaryFlagsToGeometries(
    ModelPart& rModelPart)
{
    block_for_each(rModelPart.Conditions(), [&](ModelPart::ConditionType& rCondition) {
        rCondition.SetValue(RANS_IS_INLET, rCondition.Is(INLET));
        rCondition.SetValue(RANS_IS_OUTLET, rCondition.Is(OUTLET));
        rCondition.SetValue(RANS_IS_STRUCTURE, rCondition.Is(STRUCTURE));
    });
}

template <>
KRATOS_API(RANS_APPLICATION)
double GetVariableValueNorm(const double& rValue)
{
    return rValue;
}

template <>
KRATOS_API(RANS_APPLICATION)
double GetVariableValueNorm(const array_1d<double, 3>& rValue)
{
    return norm_2(rValue);
}

template <typename TDataType>
std::tuple<double, double> CalculateTransientVariableConvergence(
    const ModelPart& rModelPart,
    const Variable<TDataType>& rVariable)
{
    KRATOS_TRY

    const auto& r_communicator = rModelPart.GetCommunicator();
    const auto& r_nodes = r_communicator.LocalMesh().Nodes();
    const int number_of_nodes = r_nodes.size();

    KRATOS_ERROR_IF(rModelPart.GetBufferSize() < 2)
        << rModelPart.Name() << " buffer size is "
        << rModelPart.GetBufferSize() << ". Buffer size of 2 or greater is required to calculate transient variable convergence for "
        << rVariable.Name() << ".\n";

    double dx, solution, number_of_dofs;
    std::tie(dx, solution, number_of_dofs) =
        IndexPartition<int>(number_of_nodes)
            .for_each<CombinedReduction<SumReduction<double>, SumReduction<double>, SumReduction<double>>>(
                [&](const int iNode) -> std::tuple<double, double, double> {
                    const auto p_node = r_nodes.begin() + iNode;
                    const auto& r_old_value =
                        p_node->FastGetSolutionStepValue(rVariable, 1);
                    const auto& r_new_value = p_node->FastGetSolutionStepValue(rVariable);
                    return std::make_tuple<double, double, double>(
                        std::pow(GetVariableValueNorm<TDataType>(r_new_value - r_old_value), 2),
                        std::pow(GetVariableValueNorm<TDataType>(r_new_value), 2),
                        ((p_node->HasDofFor(rVariable) && !p_node->IsFixed(rVariable)) ||
                         !p_node->HasDofFor(rVariable)));
                });

    // to improve mpi communication performance
    const std::vector<double> process_values = {dx, solution, number_of_dofs};
    const std::vector<double>& global_values =
        r_communicator.GetDataCommunicator().SumAll(process_values);

    dx = std::sqrt(global_values[0]);
    solution = std::sqrt(global_values[1]);
    number_of_dofs = std::max(1.0, global_values[2]);
    solution = (solution > 0.0 ? solution : 1.0);

    return std::make_tuple<double, double>(std::forward<double>(dx / solution),
                                           std::forward<double>(dx / number_of_dofs));

    KRATOS_CATCH("");
}

void SetElementConstitutiveLaws(ModelPart::ElementsContainerType& rElements)
{
    KRATOS_TRY

    block_for_each(rElements, [](ModelPart::ElementType& rElement){
        if (!rElement.Has(CONSTITUTIVE_LAW)) {
            const auto& r_properties = rElement.GetProperties();
            KRATOS_ERROR_IF_NOT(r_properties.Has(CONSTITUTIVE_LAW))
                << "In initialization of entity " << rElement.Info()
                << ": No CONSTITUTIVE_LAW defined for property "
                << r_properties.Id() << "." << std::endl;

            const auto rans_cl_name = r_properties[CONSTITUTIVE_LAW]->Info();

            KRATOS_ERROR_IF(rans_cl_name.substr(0, 4) != "Rans")
                << "Incompatible constitutive law is used. Please use constitutive "
                "laws which starts with \"Rans*\" [ Constitutive law "
                "name = "
                << rans_cl_name << " ].\n";

            // get the fluid constitutive law here because, turbulence models need the mu of fluid
            auto p_constitutive_law =
                KratosComponents<ConstitutiveLaw>::Get(rans_cl_name.substr(4)).Clone();

            const auto& r_geometry = rElement.GetGeometry();
            const auto& r_shape_functions =
                r_geometry.ShapeFunctionsValues(GeometryData::GI_GAUSS_1);
            p_constitutive_law->InitializeMaterial(r_properties, r_geometry,
                                                row(r_shape_functions, 0));

            rElement.SetValue(CONSTITUTIVE_LAW, p_constitutive_law);
        }
    });

    KRATOS_CATCH("");
}

// template instantiations
template KRATOS_API(RANS_APPLICATION) void AssignConditionVariableValuesToNodes<double>(
    ModelPart&, const Variable<double>&, const Flags&, const bool);

template KRATOS_API(RANS_APPLICATION) void AssignConditionVariableValuesToNodes<array_1d<double, 3>>(
    ModelPart&, const Variable<array_1d<double, 3>>&, const Flags&, const bool);

template KRATOS_API(RANS_APPLICATION) std::tuple<double, double> CalculateTransientVariableConvergence<double>(
    const ModelPart&, const Variable<double>&);

template KRATOS_API(RANS_APPLICATION)
    std::tuple<double, double> CalculateTransientVariableConvergence<array_1d<double, 3>>(
        const ModelPart&, const Variable<array_1d<double, 3>>&);

} // namespace RansVariableUtilities

} /* namespace Kratos.*/
