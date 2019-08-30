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
                            ModelPart::NodesContainerType& rNodes)
    {
        KRATOS_TRY

        CheckVariableExists(rVariable, rNodes);

        const int number_of_nodes = rNodes.size();

        unsigned int number_of_nodes_below_minimum = 0;
        unsigned int number_of_nodes_above_maximum = 0;
        unsigned int number_of_nodes_selected = 0;

#pragma omp parallel for reduction( +: number_of_nodes_below_minimum, number_of_nodes_above_maximum, number_of_nodes_selected)
        for (int i = 0; i < number_of_nodes; i++)
        {
            ModelPart::NodeType& r_node = *(rNodes.begin() + i);
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

        rNumberOfNodesBelowMinimum = number_of_nodes_below_minimum;
        rNumberOfNodesAboveMaximum = number_of_nodes_above_maximum;
        rNumberOfSelectedNodes = number_of_nodes_selected;

        KRATOS_CATCH("")
    }

    unsigned int GetNumberOfNegativeScalarValueNodes(const ModelPart::NodesContainerType& rNodes,
                                                     const Variable<double>& rVariable)
    {
        KRATOS_TRY

        CheckVariableExists(rVariable, rNodes);

        const int number_of_nodes = rNodes.size();
        unsigned int number_of_negative_nodes = 0;

#pragma omp parallel for reduction(+ : number_of_negative_nodes)
        for (int i = 0; i < number_of_nodes; i++)
        {
            const double value = (rNodes.begin() + i)->FastGetSolutionStepValue(rVariable);
            number_of_negative_nodes += (value < 0.0);
        }

        return number_of_negative_nodes;

        KRATOS_CATCH("");
    }

    double GetMinimumScalarValue(const ModelPart::NodesContainerType& rNodes,
                                 const Variable<double>& rVariable)
    {
        KRATOS_TRY

        CheckVariableExists(rVariable, rNodes);

        const int number_of_nodes = rNodes.size();

        if (number_of_nodes == 0)
            return 0.0;

        double min_value = rNodes.begin()->FastGetSolutionStepValue(rVariable);

        int number_of_threads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector node_partition;
        OpenMPUtils::DivideInPartitions(number_of_nodes, number_of_threads, node_partition);
        Vector min_values(number_of_threads);

#pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();

            auto NodesBegin = rNodes.begin() + node_partition[k];
            auto NodesEnd = rNodes.begin() + node_partition[k + 1];
            min_values[k] = NodesBegin->FastGetSolutionStepValue(rVariable);

            for (auto itNode = NodesBegin; itNode != NodesEnd; itNode++)
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

        return min_value;

        KRATOS_CATCH("");
    }

    double GetFlaggedMinimumScalarValue(const ModelPart::NodesContainerType& rNodes,
                                        const Variable<double>& rVariable,
                                        const Flags& rCheckFlag,
                                        const bool CheckFlagValue)
    {
        KRATOS_TRY

        CheckVariableExists(rVariable, rNodes);

        const int number_of_nodes = rNodes.size();

        if (number_of_nodes == 0)
            return 0.0;

        double min_value = rNodes.begin()->FastGetSolutionStepValue(rVariable);

        int number_of_threads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector node_partition;
        OpenMPUtils::DivideInPartitions(number_of_nodes, number_of_threads, node_partition);
        Vector min_values(number_of_threads);

#pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();

            auto NodesBegin = rNodes.begin() + node_partition[k];
            auto NodesEnd = rNodes.begin() + node_partition[k + 1];
            min_values[k] = NodesBegin->FastGetSolutionStepValue(rVariable);

            for (auto itNode = NodesBegin; itNode != NodesEnd; itNode++)
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

        return min_value;

        KRATOS_CATCH("");
    }

    double GetMaximumScalarValue(const ModelPart::NodesContainerType& rNodes,
                                 const Variable<double>& rVariable)
    {
        KRATOS_TRY

        CheckVariableExists(rVariable, rNodes);

        const int number_of_nodes = rNodes.size();

        if (number_of_nodes == 0)
            return 0.0;

        double max_value = rNodes.begin()->FastGetSolutionStepValue(rVariable);

        int number_of_threads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector node_partition;
        OpenMPUtils::DivideInPartitions(number_of_nodes, number_of_threads, node_partition);
        Vector max_values(number_of_threads);

#pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();

            auto NodesBegin = rNodes.begin() + node_partition[k];
            auto NodesEnd = rNodes.begin() + node_partition[k + 1];
            max_values[k] = NodesBegin->FastGetSolutionStepValue(rVariable);

            for (auto itNode = NodesBegin; itNode != NodesEnd; itNode++)
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

        return max_value;

        KRATOS_CATCH("");
    }

    void GetNodalVariablesVector(Vector& rValues,
                                 const ModelPart::NodesContainerType& rNodes,
                                 const Variable<double>& rVariable)
    {
        const int number_of_nodes = rNodes.size();

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

    void SetNodalVariables(const Vector& rValues,
                           ModelPart::NodesContainerType& rNodes,
                           const Variable<double>& rVariable)
    {
        const int number_of_nodes = rNodes.size();

#pragma omp parallel for
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
            (rNodes.begin() + i_node)->FastGetSolutionStepValue(rVariable) =
                rValues[i_node];
    }

    double GetScalarVariableDifferenceNormSquare(const ModelPart::NodesContainerType& rNodes,
                                                 const Variable<double>& rVariableA,
                                                 const Variable<double>& rVariableB)
    {
        KRATOS_TRY

        CheckVariableExists(rVariableA, rNodes);
        CheckVariableExists(rVariableB, rNodes);

        const int number_of_nodes = rNodes.size();
        double increase_norm_square = 0.0;

#pragma omp parallel for reduction(+ : increase_norm_square)
        for (int i = 0; i < number_of_nodes; i++)
        {
            const ModelPart::NodeType& r_node = *(rNodes.begin() + i);
            const double value_a = r_node.FastGetSolutionStepValue(rVariableA);
            const double value_b = r_node.FastGetSolutionStepValue(rVariableB);
            increase_norm_square += std::pow(value_a - value_b, 2);
        }

        return increase_norm_square;

        KRATOS_CATCH("");
    }

    double GetScalarVariableSolutionNormSquare(const ModelPart::NodesContainerType& rNodes,
                                               const Variable<double>& rVariable)
    {
        KRATOS_TRY

        CheckVariableExists(rVariable, rNodes);

        const int number_of_nodes = rNodes.size();
        double solution_norm_square = 0.0;

#pragma omp parallel for reduction(+ : solution_norm_square)
        for (int i = 0; i < number_of_nodes; i++)
        {
            const double solution_value =
                (rNodes.begin() + i)->FastGetSolutionStepValue(rVariable);
            solution_norm_square += std::pow(solution_value, 2);
        }

        return solution_norm_square;

        KRATOS_CATCH("");
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