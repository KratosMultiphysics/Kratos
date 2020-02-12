//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    ME
//
//

// System includes

// External includes

// Project includes
// See header file

// Application includes
#include "lumped_interface_positive_negative_pressure_gradient.h"


namespace Kratos
{

/* Public functions *******************************************************/

/// constructor
LumpedInterfacePositiveNegativePressureGradient::LumpedInterfacePositiveNegativePressureGradient(
        ModelPart& rModelPart)
    : Process(), mrModelPart(rModelPart) {}

/***********************************************************************************/
/***********************************************************************************/

void LumpedInterfacePositiveNegativePressureGradient::Execute(){

    KRATOS_TRY;

    // First node iterator
    const auto it_node_begin = mrModelPart.NodesBegin();

    // Set to zero
    ClearVariables();

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

        Vector distances( number_of_nodes );
        for(std::size_t i_node=0; i_node<number_of_nodes; ++i_node){
            distances(i_node) = r_geometry[i_node].FastGetSolutionStepValue(DISTANCE);
        }

        unsigned int nneg=0, npos=0;
        for(unsigned int i = 0; i < number_of_nodes; ++i)
        {
            if(distances(i) > 0) npos += 1;
            else /* if(distances(i) < 0) */ nneg += 1;
        }
        
        if(nneg == 0 || npos == 0)
        {
            // Resize if needed
            if (DN_DX.size1() != number_of_nodes || DN_DX.size2() != dimension)
                DN_DX.resize(number_of_nodes, dimension);
            if (N.size() != number_of_nodes)
                N.resize(number_of_nodes);
            if (J0.size1() != dimension || J0.size2() != local_space_dimension)
                J0.resize(dimension, local_space_dimension);

            // The integration points
            const auto& r_integration_method = r_geometry.GetDefaultIntegrationMethod();
            const auto& r_integration_points = r_geometry.IntegrationPoints(r_integration_method);
            const std::size_t number_of_integration_points = r_integration_points.size();

            Vector values(number_of_nodes);
            for(std::size_t i_node=0; i_node<number_of_nodes; ++i_node)
                values[i_node] = r_geometry[i_node].FastGetSolutionStepValue(PRESSURE);

            // The containers of the shape functions and the local gradients
            const Matrix& rNcontainer = r_geometry.ShapeFunctionsValues(r_integration_method);
            const auto& rDN_DeContainer = r_geometry.ShapeFunctionsLocalGradients(r_integration_method);

            for ( IndexType point_number = 0; point_number < number_of_integration_points; ++point_number ) {
                // Getting the shape functions
                noalias(N) = row(rNcontainer, point_number);

                // Getting the jacobians and local gradients
                GeometryUtils::JacobianOnInitialConfiguration(r_geometry, r_integration_points[point_number], J0);
                double detJ0;
                Matrix InvJ0;
                MathUtils<double>::GeneralizedInvertMatrix(J0, InvJ0, detJ0);
                const Matrix& rDN_De = rDN_DeContainer[point_number];
                GeometryUtils::ShapeFunctionsGradients(rDN_De, InvJ0, DN_DX);

                const Vector grad = prod(trans(DN_DX), values);
                const double gauss_point_volume = r_integration_points[point_number].Weight() * detJ0;

                for(std::size_t i_node=0; i_node<number_of_nodes; ++i_node) {
                    array_1d<double, 3>& r_gradient = r_geometry[i_node].FastGetSolutionStepValue(PRESSURE_GRADIENT_AUX);
                    for(std::size_t k=0; k<dimension; ++k) {
                        #pragma omp atomic
                        r_gradient[k] += N[i_node] * gauss_point_volume*grad[k];
                    }

                    double& vol = r_geometry[i_node].GetValue(AREA_VARIABLE_AUX);

                    #pragma omp atomic
                    vol += N[i_node] * gauss_point_volume;
                }
            }   
        }        
    }

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i) {
        auto it_node = it_node_begin + i;
        if (it_node->GetValue(AREA_VARIABLE_AUX) > 1.0e-15)
            it_node->FastGetSolutionStepValue(PRESSURE_GRADIENT_AUX) /= it_node->GetValue(AREA_VARIABLE_AUX);
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void LumpedInterfacePositiveNegativePressureGradient::ClearVariables()
{
    const auto it_node_begin = mrModelPart.NodesBegin();

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i) {
        auto it_node=it_node_begin + i;
        it_node->SetValue(AREA_VARIABLE_AUX, 0.0);
        it_node->FastGetSolutionStepValue(PRESSURE_GRADIENT_AUX).clear();
        /* it_node->SetValue(IS_NEAR_CUT, 0.0); */
    }
}

};  // namespace Kratos.
