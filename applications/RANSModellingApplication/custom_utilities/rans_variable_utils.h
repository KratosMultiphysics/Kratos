//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//
//

#if !defined(KRATOS_RANS_VARIABLE_UTILS)
#define KRATOS_RANS_VARIABLE_UTILS

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "includes/model_part.h"
#include "rans_modelling_application_variables.h"
#include "utilities/openmp_utils.h"
#include "utilities/variable_utils.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class RansVariableUtils
 * @ingroup KratosRANSModellingApplication
 * @brief This class extends @KratosCore VariableUtils class to perform some RANS unique
 * tasks related with the variable values.
 * @details The methods are exported to python in order to add this improvements to the python interface
 * @author Suneth Warnakulasuriya
 * @see VariableUtils
 */
class RansVariableUtils : public VariableUtils
{
public:
    ///@name Type Definitions
    ///@{

    /// We create the Pointer related to RansVariableUtils
    KRATOS_CLASS_POINTER_DEFINITION(RansVariableUtils);

    void ClipScalarVariable(unsigned int& rNumberOfNodesBelowMinimum,
                            unsigned int& rNumberOfNodesAboveMaximum,
                            unsigned int& rNumberOfSelectedNodes,
                            const double MinimumValue,
                            const double MaximumValue,
                            const Variable<double>& rVariable,
                            ModelPart& rModelPart) const
    {
        KRATOS_TRY

        Communicator& r_communicator = rModelPart.GetCommunicator();
        ModelPart::NodesContainerType& r_nodes = r_communicator.LocalMesh().Nodes();

        const int number_of_nodes = r_nodes.size();

        unsigned int number_of_nodes_below_minimum = 0;
        unsigned int number_of_nodes_above_maximum = 0;
        unsigned int number_of_nodes_selected = 0;

#pragma omp parallel for reduction( +: number_of_nodes_below_minimum, number_of_nodes_above_maximum, number_of_nodes_selected)
        for (int i = 0; i < number_of_nodes; i++)
        {
            ModelPart::NodeType& r_node = *(r_nodes.begin() + i);
            double& r_value = r_node.FastGetSolutionStepValue(rVariable);

            if (r_value < MinimumValue)
            {
                number_of_nodes_below_minimum++;
                r_value = MinimumValue;
            }
            else if (r_value > MaximumValue)
            {
                number_of_nodes_above_maximum++;
                r_value = MaximumValue;
            }
            number_of_nodes_selected++;
        }

        r_communicator.SynchronizeVariable(rVariable);

        // Stores followings
        // index - 0 : number_of_nodes_below_minimum
        // index - 1 : number_of_nodes_above_maximum
        // index - 2 : number_of_nodes_selected
        std::vector<unsigned int> nodes_count = {number_of_nodes_below_minimum,
                                                 number_of_nodes_above_maximum,
                                                 number_of_nodes_selected};
        const std::vector<unsigned int>& total_nodes_count =
            r_communicator.GetDataCommunicator().SumAll(nodes_count);

        rNumberOfNodesBelowMinimum = total_nodes_count[0];
        rNumberOfNodesAboveMaximum = total_nodes_count[1];
        rNumberOfSelectedNodes = total_nodes_count[2];

        KRATOS_CATCH("")
    }

    double GetMinimumScalarValue(const ModelPart& rModelPart,
                                 const Variable<double>& rVariable) const
    {
        KRATOS_TRY

        double min_value = std::numeric_limits<double>::max();

        const Communicator& r_communicator = rModelPart.GetCommunicator();

        const ModelPart::NodesContainerType& r_nodes =
            r_communicator.LocalMesh().Nodes();

        const int number_of_nodes = r_nodes.size();

        if (number_of_nodes != 0)
        {
            int number_of_threads = OpenMPUtils::GetNumThreads();
            OpenMPUtils::PartitionVector node_partition;
            OpenMPUtils::DivideInPartitions(number_of_nodes, number_of_threads, node_partition);
            Vector min_values(number_of_threads);

#pragma omp parallel
            {
                int k = OpenMPUtils::ThisThread();

                auto nodes_begin = r_nodes.begin() + node_partition[k];
                auto nodes_end = r_nodes.begin() + node_partition[k + 1];
                min_values[k] = nodes_begin->FastGetSolutionStepValue(rVariable);

                for (auto itNode = nodes_begin; itNode != nodes_end; itNode++)
                {
                    const double value = itNode->FastGetSolutionStepValue(rVariable);
                    min_values[k] = std::min(min_values[k], value);
                }

#pragma omp critical
                {
                    for (int i = 0; i < number_of_threads; ++i)
                    {
                        min_value = std::min(min_value, min_values[i]);
                    }
                }
            }
        }

        return r_communicator.GetDataCommunicator().MinAll(min_value);

        KRATOS_CATCH("");
    }

    double GetFlaggedMinimumScalarValue(const ModelPart& rModelPart,
                                        const Variable<double>& rVariable,
                                        const Flags& rCheckFlag,
                                        const bool CheckFlagValue) const
    {
        KRATOS_TRY

        const Communicator& r_communicator = rModelPart.GetCommunicator();
        const ModelPart::NodesContainerType& r_nodes =
            r_communicator.LocalMesh().Nodes();

        double min_value = std::numeric_limits<double>::max();

        const int number_of_nodes = r_nodes.size();

        if (number_of_nodes != 0)
        {
            int number_of_threads = OpenMPUtils::GetNumThreads();
            OpenMPUtils::PartitionVector node_partition;
            OpenMPUtils::DivideInPartitions(number_of_nodes, number_of_threads, node_partition);
            Vector min_values(number_of_threads);

#pragma omp parallel
            {
                int k = OpenMPUtils::ThisThread();

                auto nodes_begin = r_nodes.begin() + node_partition[k];
                auto nodes_end = r_nodes.begin() + node_partition[k + 1];
                min_values[k] = nodes_begin->FastGetSolutionStepValue(rVariable);

                for (auto itNode = nodes_begin; itNode != nodes_end; itNode++)
                {
                    if ((itNode->Is(rCheckFlag) && CheckFlagValue) ||
                        (!itNode->Is(rCheckFlag) && !CheckFlagValue))
                    {
                        const double value = itNode->FastGetSolutionStepValue(rVariable);
                        min_values[k] = std::min(min_values[k], value);
                    }
                }

#pragma omp critical
                {
                    for (int i = 0; i < number_of_threads; ++i)
                    {
                        min_value = std::min(min_value, min_values[i]);
                    }
                }
            }
        }

        return r_communicator.GetDataCommunicator().MinAll(min_value);

        KRATOS_CATCH("");
    }

    double GetMaximumScalarValue(const ModelPart& rModelPart,
                                 const Variable<double>& rVariable) const
    {
        KRATOS_TRY

        double max_value = std::numeric_limits<double>::min();

        const Communicator& r_communicator = rModelPart.GetCommunicator();

        const ModelPart::NodesContainerType& r_nodes =
            r_communicator.LocalMesh().Nodes();

        const int number_of_nodes = r_nodes.size();

        if (number_of_nodes != 0)
        {
            int number_of_threads = OpenMPUtils::GetNumThreads();
            OpenMPUtils::PartitionVector node_partition;
            OpenMPUtils::DivideInPartitions(number_of_nodes, number_of_threads, node_partition);
            Vector max_values(number_of_threads);

#pragma omp parallel
            {
                int k = OpenMPUtils::ThisThread();

                auto nodes_begin = r_nodes.begin() + node_partition[k];
                auto nodes_end = r_nodes.begin() + node_partition[k + 1];
                max_values[k] = nodes_begin->FastGetSolutionStepValue(rVariable);

                for (auto itNode = nodes_begin; itNode != nodes_end; itNode++)
                {
                    const double value = itNode->FastGetSolutionStepValue(rVariable);
                    max_values[k] = std::max(max_values[k], value);
                }

#pragma omp critical
                {
                    for (int i = 0; i < number_of_threads; ++i)
                    {
                        max_value = std::max(max_value, max_values[i]);
                    }
                }
            }
        }

        return r_communicator.GetDataCommunicator().MaxAll(max_value);

        KRATOS_CATCH("");
    }

    double GetFlaggedMaximumScalarValue(const ModelPart& rModelPart,
                                        const Variable<double>& rVariable,
                                        const Flags& rCheckFlag,
                                        const bool CheckFlagValue) const
    {
        KRATOS_TRY

        const Communicator& r_communicator = rModelPart.GetCommunicator();
        const ModelPart::NodesContainerType& r_nodes =
            r_communicator.LocalMesh().Nodes();

        double max_value = std::numeric_limits<double>::min();

        const int number_of_nodes = r_nodes.size();

        if (number_of_nodes != 0)
        {
            int number_of_threads = OpenMPUtils::GetNumThreads();
            OpenMPUtils::PartitionVector node_partition;
            OpenMPUtils::DivideInPartitions(number_of_nodes, number_of_threads, node_partition);
            Vector max_values(number_of_threads);

#pragma omp parallel
            {
                int k = OpenMPUtils::ThisThread();

                auto nodes_begin = r_nodes.begin() + node_partition[k];
                auto nodes_end = r_nodes.begin() + node_partition[k + 1];
                max_values[k] = nodes_begin->FastGetSolutionStepValue(rVariable);

                for (auto itNode = nodes_begin; itNode != nodes_end; itNode++)
                {
                    if ((itNode->Is(rCheckFlag) && CheckFlagValue) ||
                        (!itNode->Is(rCheckFlag) && !CheckFlagValue))
                    {
                        const double value = itNode->FastGetSolutionStepValue(rVariable);
                        max_values[k] = std::max(max_values[k], value);
                    }
                }

#pragma omp critical
                {
                    for (int i = 0; i < number_of_threads; ++i)
                    {
                        max_value = std::max(max_value, max_values[i]);
                    }
                }
            }
        }

        return r_communicator.GetDataCommunicator().MaxAll(max_value);

        KRATOS_CATCH("");
    }

    void GetNodalVariablesVector(Vector& rValues,
                                 const ModelPart::NodesContainerType& rNodes,
                                 const Variable<double>& rVariable) const
    {
        const int number_of_nodes = rNodes.size();

        if (static_cast<int>(rValues.size()) != number_of_nodes)
            rValues.resize(number_of_nodes);

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
            rValues[i_node] =
                (rNodes.begin() + i_node)->FastGetSolutionStepValue(rVariable);
    }
    void GetNodalArray(Vector& rNodalValues, const Element& rElement, const Variable<double>& rVariable)
    {
        const Geometry<ModelPart::NodeType>& r_geometry = rElement.GetGeometry();
        std::size_t number_of_nodes = r_geometry.PointsNumber();

        if (rNodalValues.size() != number_of_nodes)
            rNodalValues.resize(number_of_nodes);

        for (std::size_t i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            rNodalValues[i_node] = r_geometry[i_node].FastGetSolutionStepValue(rVariable);
        }
    }

    void SetNodalVariables(ModelPart::NodesContainerType& rNodes,
                           const Vector& rValues,
                           const Variable<double>& rVariable) const
    {
        const int number_of_nodes = rNodes.size();

        KRATOS_ERROR_IF(static_cast<int>(rValues.size()) != number_of_nodes)
            << "rValues vector size mismatch with rNodes size in "
               "SetNodalVariables. [ rValues.size = "
            << rValues.size() << ", rNodes.size = " << rNodes.size() << " ]\n";

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
            (rNodes.begin() + i_node)->FastGetSolutionStepValue(rVariable) =
                rValues[i_node];
    }

    void CopyNodalSolutionStepVariablesList(ModelPart& rOriginModelPart, ModelPart& rDestinationModelPart)
    {
        KRATOS_TRY

        rDestinationModelPart.GetNodalSolutionStepVariablesList() =
            rOriginModelPart.GetNodalSolutionStepVariablesList();

        KRATOS_CATCH("");
    }
}; /* Class RansVariableUtils */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_RANS_VARIABLE_UTILS  defined */