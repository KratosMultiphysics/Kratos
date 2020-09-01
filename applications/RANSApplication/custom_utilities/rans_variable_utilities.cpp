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

/* System includes */
#include <cmath>
#include <limits>

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "utilities/parallel_utilities.h"

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

    class ClippingReducer{
        public:
            typedef std::tuple<unsigned int,unsigned int> value_type;
            unsigned int number_of_nodes_below_minimum = 0;
            unsigned int number_of_nodes_above_maximum = 0;

            value_type GetValue()
            {
                value_type values;
                std::get<0>(values) = number_of_nodes_below_minimum;
                std::get<1>(values) = number_of_nodes_above_maximum;
                return values;
            }

            void LocalReduce(const std::tuple<unsigned int, unsigned int>& node_values){
                number_of_nodes_below_minimum += std::get<0>(node_values);
                number_of_nodes_above_maximum += std::get<1>(node_values);
            }
            void ThreadSafeReduce(ClippingReducer& rOther){
                #pragma omp critical
                {
                    this->number_of_nodes_below_minimum += rOther.number_of_nodes_below_minimum;
                    this->number_of_nodes_above_maximum += rOther.number_of_nodes_above_maximum;
                }
            }
    };

    unsigned int number_of_nodes_below_minimum, number_of_nodes_above_maximum;
    std::tie(number_of_nodes_below_minimum, number_of_nodes_above_maximum) = BlockPartition<ModelPart::NodesContainerType>(r_nodes).for_each<ClippingReducer>(
        [&](ModelPart::NodeType& rNode) -> std::tuple<unsigned int, unsigned int> {
            double& r_value = rNode.FastGetSolutionStepValue(rVariable);

            if (r_value < MinimumValue) {
                r_value = MinimumValue;
                return std::tuple<unsigned int, unsigned int>(1, 0);
            } else if (r_value > MaximumValue) {
                r_value = MaximumValue;
                return std::tuple<unsigned int, unsigned int>(0, 1);
            }

            return std::tuple<unsigned int, unsigned int>(0, 0);
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

    class MinReducer{
        public:
            typedef double value_type;
            double min_value = std::numeric_limits<double>::max();

            value_type GetValue()
            {
                return min_value;
            }

            void LocalReduce(const double min_value){
                this->min_value = std::min(this->min_value, min_value);
            }
            void ThreadSafeReduce(MinReducer& rOther){
                #pragma omp critical
                {
                    this->min_value = std::min(this->min_value, rOther.min_value);
                }
            }
    };

    const int number_of_nodes = r_nodes.size();
    const double min_value =
        IndexPartition<int>(number_of_nodes).for_each<MinReducer>(
            [&](const int i) -> double {
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

    class MaxReducer{
        public:
            typedef double value_type;
            double max_value = std::numeric_limits<double>::lowest();

            value_type GetValue()
            {
                return max_value;
            }

            void LocalReduce(const double max_value){
                this->max_value = std::max(this->max_value, max_value);
            }
            void ThreadSafeReduce(MaxReducer& rOther){
                #pragma omp critical
                {
                    this->max_value = std::max(this->max_value, rOther.max_value);
                }
            }
    };

    const int number_of_nodes = r_nodes.size();
    const double max_value =
        IndexPartition<int>(number_of_nodes).for_each<MaxReducer>(
            [&](const int i) -> double {
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
    const Geometry<ModelPart::NodeType>& r_geometry = rElement.GetGeometry();
    const std::size_t number_of_nodes = r_geometry.PointsNumber();

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
} // namespace RansVariableUtilities

} /* namespace Kratos.*/
