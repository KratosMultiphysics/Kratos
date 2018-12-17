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
//  Collaborators:   Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "utilities/geometry_utilities.h"
#include "processes/calculate_nodal_area_process.h"

namespace Kratos
{

template<bool THistorical>
void CalculateNodalAreaProcess<THistorical>::Execute()
{
    KRATOS_TRY

    // Check if variables are available
    KRATOS_ERROR_IF(mrModelPart.Nodes().size() == 0) << "No nodes in the model part" << std::endl;
    if (THistorical) {
        KRATOS_ERROR_IF_NOT(mrModelPart.NodesBegin()->SolutionStepsDataHas( NODAL_AREA )) << "Variable NODAL_AREA not in the model part!" << std::endl;
    }

    // Set to zero the nodal area
    const auto it_node_begin = mrModelPart.NodesBegin();
    #pragma omp parallel for
    for(int i=0; i<static_cast<int>(mrModelPart.Nodes().size()); ++i) {
        auto it_node = it_node_begin + i;
        SetInitialValue(it_node);
    }

    const auto& it_element_begin = mrModelPart.GetCommunicator().LocalMesh().ElementsBegin();
    const auto& r_first_element_geometry = it_element_begin->GetGeometry();
    const std::size_t local_space_dimension = r_first_element_geometry.LocalSpaceDimension();
    const std::size_t number_of_nodes = r_first_element_geometry.PointsNumber();

    // The integration points
    const auto& integration_method = r_first_element_geometry.GetDefaultIntegrationMethod();
    const auto& integration_points = r_first_element_geometry.IntegrationPoints(integration_method);
    const std::size_t number_of_integration_points = integration_points.size();

    Vector N = ZeroVector(number_of_nodes);
    Matrix J0 = ZeroMatrix(mDomainSize, local_space_dimension);

    #pragma omp parallel for firstprivate(N, J0)
    for(int i=0; i<static_cast<int>(mrModelPart.Elements().size()); ++i) {
        auto it_elem = it_element_begin + i;
        auto& r_geometry = it_elem->GetGeometry();

        // The containers of the shape functions
        const auto& rNcontainer = r_geometry.ShapeFunctionsValues(integration_method);

        for ( IndexType point_number = 0; point_number < number_of_integration_points; ++point_number ) {
            // Getting the shape functions
            noalias(N) = row(rNcontainer, point_number);

            // Getting the jacobians and local gradients
            GeometryUtils::JacobianOnInitialConfiguration(r_geometry, integration_points[point_number], J0);
            const double detJ0 = MathUtils<double>::GeneralizedDet(J0);
            const double gauss_point_volume = integration_points[point_number].Weight() * detJ0;

            for(std::size_t i_node =0; i_node < number_of_nodes; ++i_node) {
                double& nodal_area = GetAreaValue(r_geometry[i_node]);
                #pragma omp atomic
                nodal_area += N[i_node] * gauss_point_volume;
            }
        }
    }

    // Synchronize data
    if (THistorical) {
        mrModelPart.GetCommunicator().AssembleCurrentData(NODAL_AREA);
    } else {
        mrModelPart.GetCommunicator().AssembleNonHistoricalData(NODAL_AREA);
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<>
double& CalculateNodalAreaProcess<true>::GetAreaValue(NodeType& rNode)
{
    return rNode.FastGetSolutionStepValue(NODAL_AREA);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
double& CalculateNodalAreaProcess<false>::GetAreaValue(NodeType& rNode)
{
    return rNode.GetValue(NODAL_AREA);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void CalculateNodalAreaProcess<true>::SetAreaValue(
    NodeType& rNode,
    const double Value
    )
{
    rNode.FastGetSolutionStepValue(NODAL_AREA) = Value;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void CalculateNodalAreaProcess<false>::SetAreaValue(
    NodeType& rNode,
    const double Value
    )
{
    rNode.SetValue(NODAL_AREA, Value);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void CalculateNodalAreaProcess<true>::SetInitialValue(NodeIterator itNode)
{
    itNode->FastGetSolutionStepValue(NODAL_AREA) = 0.0;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void CalculateNodalAreaProcess<false>::SetInitialValue(NodeIterator itNode)
{
    itNode->SetValue(NODAL_AREA, 0.0);
}

/***********************************************************************************/
/***********************************************************************************/

template class CalculateNodalAreaProcess<true>;
template class CalculateNodalAreaProcess<false>;

}  // namespace Kratos.


