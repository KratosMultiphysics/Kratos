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
#include "utilities/variable_utils.h"
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
    if (THistorical) {
        VariableUtils().SetVariable(NODAL_AREA, 0.0, mrModelPart.Nodes());
    } else {
        VariableUtils().SetNonHistoricalVariable(NODAL_AREA, 0.0, mrModelPart.Nodes());
    }

    // We create the N and J0
    struct tls_type
    {
        Vector N;
        Matrix J0;
    };

    //Using Parallel Utilities
    block_for_each(mrModelPart.GetCommunicator().LocalMesh().Elements(), tls_type(),[&](Element& rElem, tls_type& rTLS){
        auto& r_geometry = rElem.GetGeometry();

        const std::size_t local_space_dimension = r_geometry.LocalSpaceDimension();
        const std::size_t number_of_nodes = r_geometry.PointsNumber();

        // The integration points
        const auto& integration_method = r_geometry.GetDefaultIntegrationMethod();
        const auto& integration_points = r_geometry.IntegrationPoints(integration_method);
        const std::size_t number_of_integration_points = integration_points.size();

        // Resize the N and J0
        if ((rTLS.N).size() != number_of_nodes)
            (rTLS.N).resize(number_of_nodes);
        if ((rTLS.J0).size1() != mDomainSize || (rTLS.J0).size2() != local_space_dimension)
            (rTLS.J0).resize(mDomainSize, local_space_dimension);

        // The containers of the shape functions
        const auto& rNcontainer = r_geometry.ShapeFunctionsValues(integration_method);

        for ( IndexType point_number = 0; point_number < number_of_integration_points; ++point_number ) {
            // Getting the shape functions
            noalias((rTLS.N)) = row(rNcontainer, point_number);

            // Getting the jacobians and local gradients
            GeometryUtils::JacobianOnInitialConfiguration(r_geometry, integration_points[point_number], (rTLS.J0));
            const double detJ0 = MathUtils<double>::GeneralizedDet((rTLS.J0));
            const double gauss_point_volume = integration_points[point_number].Weight() * detJ0;

            for(std::size_t i_node =0; i_node < number_of_nodes; ++i_node) {
                double& nodal_area = GetAreaValue(r_geometry[i_node]);
                AtomicAdd(nodal_area, ((rTLS.N)[i_node] * gauss_point_volume));
            }
        }
    });

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

template class CalculateNodalAreaProcess<true>;
template class CalculateNodalAreaProcess<false>;

}  // namespace Kratos.


