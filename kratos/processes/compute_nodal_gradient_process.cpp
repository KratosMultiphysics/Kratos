//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Vicente Mataix Ferrandiz
//
//

/* System includes */

/* External includes */

/* Project includes */
#include "utilities/variable_utils.h"
#include "utilities/geometry_utilities.h"
#include "processes/compute_nodal_gradient_process.h"

namespace Kratos
{
template<bool THistorical>
void ComputeNodalGradientProcess<THistorical>::Execute()
{
    KRATOS_TRY;

    // Set to zero
    ClearGradient();

    // Auxiliar containers
    Matrix DN_DX, J0;
    Vector N;

    // First element iterator
    const auto it_element_begin = mrModelPart.ElementsBegin();

    // Current domain size
    const std::size_t dimension = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];

    // Iterate over the elements
    #pragma omp parallel for firstprivate(DN_DX,  N, J0)
    for(int i_elem=0; i_elem<static_cast<int>(mrModelPart.Elements().size()); ++i_elem) {
        auto it_elem = it_element_begin + i_elem;
        auto& r_geometry = it_elem->GetGeometry();

        // Current geometry information
        const std::size_t local_space_dimension = r_geometry.LocalSpaceDimension();
        const std::size_t number_of_nodes = r_geometry.PointsNumber();

        // Resize if needed
        if (DN_DX.size1() != number_of_nodes || DN_DX.size2() != dimension)
            DN_DX.resize(number_of_nodes, dimension);
        if (N.size() != number_of_nodes)
            N.resize(number_of_nodes);
        if (J0.size1() != dimension || J0.size2() != local_space_dimension)
            J0.resize(dimension, local_space_dimension);

        // The integration points
        const auto& integration_method = r_geometry.GetDefaultIntegrationMethod();
        const auto& integration_points = r_geometry.IntegrationPoints(integration_method);
        const std::size_t number_of_integration_points = integration_points.size();

        Vector values(number_of_nodes);
        if (mrOriginVariableDoubleList.size() > 0) {
            for(std::size_t i_node=0; i_node<number_of_nodes; ++i_node)
                values[i_node] = r_geometry[i_node].FastGetSolutionStepValue(*mrOriginVariableDoubleList[0]);
        } else {
            for(std::size_t i_node=0; i_node<number_of_nodes; ++i_node)
                values[i_node] = r_geometry[i_node].FastGetSolutionStepValue(*mrOriginVariableComponentsList[0]);
        }

        // The containers of the shape functions and the local gradients
        const auto& rNcontainer = r_geometry.ShapeFunctionsValues(integration_method);
        const auto& rDN_DeContainer = r_geometry.ShapeFunctionsLocalGradients(integration_method);

        for ( IndexType point_number = 0; point_number < number_of_integration_points; ++point_number ) {
            // Getting the shape functions
            noalias(N) = row(rNcontainer, point_number);

            // Getting the jacobians and local gradients
            GeometryUtils::JacobianOnInitialConfiguration(r_geometry, integration_points[point_number], J0);
            double detJ0;
            Matrix InvJ0;
            MathUtils<double>::InvertMatrix(J0, InvJ0, detJ0);
            const Matrix& rDN_De = rDN_DeContainer[point_number];
            GeometryUtils::ShapeFunctionsGradients(rDN_De, InvJ0, DN_DX);

            const Vector grad = prod(trans(DN_DX), values);
            const double gauss_point_volume = integration_points[point_number].Weight() * detJ0;

            for(std::size_t i_node=0; i_node<number_of_nodes; ++i_node) {
                array_1d<double, 3>& r_gradient = GetGradient(r_geometry, i_node);
                for(std::size_t k=0; k<dimension; ++k) {
                    #pragma omp atomic
                    r_gradient[k] += N[i_node] * gauss_point_volume*grad[k];
                }

                double& vol = r_geometry[i_node].GetValue(mrAreaVariable);

                #pragma omp atomic
                vol += N[i_node] * gauss_point_volume;
            }
        }
    }

    PonderateGradient();

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<>
ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsHistoricalVariable>::ComputeNodalGradientProcess(
    ModelPart& rModelPart,
    Variable<double>& rOriginVariable,
    Variable<array_1d<double,3> >& rGradientVariable,
    Variable<double>& rAreaVariable)
    :mrModelPart(rModelPart), mrGradientVariable(rGradientVariable), mrAreaVariable(rAreaVariable)
{
    KRATOS_TRY

    // We push the list of double variables
    mrOriginVariableDoubleList.push_back(&rOriginVariable);

    VariableUtils().CheckVariableExists(rOriginVariable, mrModelPart.Nodes());
    VariableUtils().CheckVariableExists(rGradientVariable, mrModelPart.Nodes());
    // In case the area or gradient variable is not initialized we initialize it
    auto& r_nodes = rModelPart.Nodes();
    if (!r_nodes.begin()->Has( rAreaVariable )) {
        VariableUtils().SetNonHistoricalVariable(rAreaVariable, 0.0, r_nodes);
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<>
ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsNonHistoricalVariable>::ComputeNodalGradientProcess(
    ModelPart& rModelPart,
    Variable<double>& rOriginVariable,
    Variable<array_1d<double,3> >& rGradientVariable,
    Variable<double>& rAreaVariable)
    :mrModelPart(rModelPart), mrGradientVariable(rGradientVariable), mrAreaVariable(rAreaVariable)
{
    KRATOS_TRY

    // We push the list of double variables
    mrOriginVariableDoubleList.push_back(&rOriginVariable);

    VariableUtils().CheckVariableExists(rOriginVariable, mrModelPart.Nodes());
    // In case the area or gradient variable is not initialized we initialize it
    auto& r_nodes = rModelPart.Nodes();
    if (!r_nodes.begin()->Has( rGradientVariable )) {
        const array_1d<double,3> zero_vector = ZeroVector(3);
        VariableUtils().SetNonHistoricalVariable(rGradientVariable, zero_vector, r_nodes);
    }
    if (!r_nodes.begin()->Has( rAreaVariable )) {
        VariableUtils().SetNonHistoricalVariable(rAreaVariable, 0.0, r_nodes);
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<>
ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsHistoricalVariable>::ComputeNodalGradientProcess(
    ModelPart& rModelPart,
    ComponentType& rOriginVariable,
    Variable<array_1d<double,3> >& rGradientVariable,
    Variable<double>& rAreaVariable)
    :mrModelPart(rModelPart), mrGradientVariable(rGradientVariable), mrAreaVariable(rAreaVariable)
{
    KRATOS_TRY

    // We push the components list
    mrOriginVariableComponentsList.push_back(&rOriginVariable);

    VariableUtils().CheckVariableExists(rOriginVariable, mrModelPart.Nodes());
    VariableUtils().CheckVariableExists(rGradientVariable, mrModelPart.Nodes());
    // In case the area or gradient variable is not initialized we initialize it
    auto& r_nodes = rModelPart.Nodes();
    if (!r_nodes.begin()->Has( rAreaVariable )) {
        VariableUtils().SetNonHistoricalVariable(rAreaVariable, 0.0, r_nodes);
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<>
ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsNonHistoricalVariable>::ComputeNodalGradientProcess(
    ModelPart& rModelPart,
    ComponentType& rOriginVariable,
    Variable<array_1d<double,3> >& rGradientVariable,
    Variable<double>& rAreaVariable)
    :mrModelPart(rModelPart), mrGradientVariable(rGradientVariable), mrAreaVariable(rAreaVariable)
{
    KRATOS_TRY

    // We push the components list
    mrOriginVariableComponentsList.push_back(&rOriginVariable);

    VariableUtils().CheckVariableExists(rOriginVariable, mrModelPart.Nodes());
    // In case the area or gradient variable is not initialized we initialize it
    auto& r_nodes = rModelPart.Nodes();
    if (!r_nodes.begin()->Has( rGradientVariable )) {
        const array_1d<double,3> zero_vector = ZeroVector(3);
        VariableUtils().SetNonHistoricalVariable(rGradientVariable, zero_vector, r_nodes);
    }
    if (!r_nodes.begin()->Has( rAreaVariable )) {
        VariableUtils().SetNonHistoricalVariable(rAreaVariable, 0.0, r_nodes);
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsHistoricalVariable>::ClearGradient()
{
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i) {
        auto it_node=mrModelPart.NodesBegin()+i;
        it_node->SetValue(mrAreaVariable, 0.0);
        it_node->FastGetSolutionStepValue(mrGradientVariable).clear();
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsNonHistoricalVariable>::ClearGradient()
{
    const array_1d<double, 3> aux_zero_vector = ZeroVector(3);

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i) {
        auto it_node=mrModelPart.NodesBegin()+i;
        it_node->SetValue(mrAreaVariable, 0.0);
        it_node->SetValue(mrGradientVariable, aux_zero_vector);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <>
array_1d<double, 3>& ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsHistoricalVariable>::GetGradient(
    Element::GeometryType& rThisGeometry,
    unsigned int i
    )
{
    array_1d<double, 3>& val = rThisGeometry[i].FastGetSolutionStepValue(mrGradientVariable);
    return val;
}

/***********************************************************************************/
/***********************************************************************************/

template <>
array_1d<double, 3>& ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsNonHistoricalVariable>::GetGradient(
    Element::GeometryType& rThisGeometry,
    unsigned int i
    )
{
    array_1d<double, 3>& val = rThisGeometry[i].GetValue(mrGradientVariable);
    return val;
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsHistoricalVariable>::PonderateGradient()
{
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i) {
        auto it_node = mrModelPart.NodesBegin()+i;
        it_node->FastGetSolutionStepValue(mrGradientVariable) /= it_node->GetValue(mrAreaVariable);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <>
void ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsNonHistoricalVariable>::PonderateGradient()
{
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i)
    {
        auto it_node=mrModelPart.NodesBegin()+i;
        it_node->GetValue(mrGradientVariable) /= it_node->GetValue(mrAreaVariable);
    }
}
/***********************************************************************************/
/***********************************************************************************/

template class ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsHistoricalVariable>;
template class ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsNonHistoricalVariable>;

} /* namespace Kratos.*/
