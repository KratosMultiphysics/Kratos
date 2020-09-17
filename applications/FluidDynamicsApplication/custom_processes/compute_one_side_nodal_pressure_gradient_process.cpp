//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Mohammad R. Hashemi
//

// Project includes
#include "includes/define.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "compute_one_side_nodal_pressure_gradient_process.h"


namespace Kratos
{

/* Public functions *******************************************************/

/// constructor
ComputeOneSideNodalPressureGradientProcess::ComputeOneSideNodalPressureGradientProcess(
        ModelPart& rModelPart)
    : Process(), mrModelPart(rModelPart) {}

/// Constructor with Kratos parameters.
ComputeOneSideNodalPressureGradientProcess::ComputeOneSideNodalPressureGradientProcess(
    ModelPart& rModelPart,
    Parameters Parameters)
    : ComputeOneSideNodalPressureGradientProcess(
        rModelPart
    ){}

/// Constructor with Kratos model
ComputeOneSideNodalPressureGradientProcess::ComputeOneSideNodalPressureGradientProcess(
    Model& rModel,
    Parameters Parameters)
    : ComputeOneSideNodalPressureGradientProcess(
        rModel.GetModelPart(Parameters["model_part_name"].GetString())
    ){}

void ComputeOneSideNodalPressureGradientProcess::Execute(){

    KRATOS_TRY;

    // Set to zero
    block_for_each(mrModelPart.Nodes(), [&](Node<3>& rNode){
        rNode.SetValue(NODAL_AREA, 0.0);
        rNode.FastGetSolutionStepValue(PRESSURE_GRADIENT).clear();
    });

    // Current domain size
    const unsigned int num_dim = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];
    const unsigned int num_nodes = num_dim + 1; //For tetrahedra and triangles

    // Auxiliar containers
    double detJ0 = 0.0;
    Matrix DN_DX, J0, InvJ0;
    Vector N(num_nodes), pressures(num_nodes), distances(num_nodes), grad(num_dim);

    // First element iterator
    const auto it_element_begin = mrModelPart.ElementsBegin();

    // Iterate over the elements
    #pragma omp parallel for firstprivate(DN_DX, N, J0, InvJ0, detJ0, pressures, distances, grad)
    for(int i_elem=0; i_elem<static_cast<int>(mrModelPart.Elements().size()); ++i_elem) {

        auto it_elem = it_element_begin + i_elem;
        auto& r_geometry = it_elem->GetGeometry();

        // Current geometry information
        const unsigned number_of_nodes = r_geometry.PointsNumber();

        for(std::size_t i_node=0; i_node<number_of_nodes; ++i_node){
            distances(i_node) = r_geometry[i_node].FastGetSolutionStepValue(DISTANCE);
        }

        unsigned int nneg=0, npos=0;
        for(unsigned int i = 0; i < number_of_nodes; ++i)
        {
            if(distances(i) > 0) {
                npos += 1;
            } else {
                nneg += 1;
            }
        }

        if(nneg == 0 || npos == 0)
        {
            if (N.size() != number_of_nodes){
                N.resize(number_of_nodes);
                pressures.resize(number_of_nodes);
                distances.resize(number_of_nodes);
            }

            // The integration points
            const auto& r_integration_method = r_geometry.GetDefaultIntegrationMethod();
            const auto& r_integration_points = r_geometry.IntegrationPoints(r_integration_method);
            const unsigned int number_of_integration_points = r_integration_points.size();

            for(std::size_t i_node=0; i_node<number_of_nodes; ++i_node)
                pressures[i_node] = r_geometry[i_node].FastGetSolutionStepValue(PRESSURE);

            // The containers of the shape functions and the local gradients
            const auto& rNcontainer = r_geometry.ShapeFunctionsValues(r_integration_method);
            const auto& rDN_DeContainer = r_geometry.ShapeFunctionsLocalGradients(r_integration_method);

            for ( IndexType point_number = 0; point_number < number_of_integration_points; ++point_number ) {
                // Getting the shape functions
                noalias(N) = row(rNcontainer, point_number);

                // Getting the jacobians and local gradients
                GeometryUtils::JacobianOnInitialConfiguration(r_geometry, r_integration_points[point_number], J0);
                MathUtils<double>::GeneralizedInvertMatrix(J0, InvJ0, detJ0);
                const auto& rDN_De = rDN_DeContainer[point_number];
                GeometryUtils::ShapeFunctionsGradients(rDN_De, InvJ0, DN_DX);

                noalias(grad) = prod(trans(DN_DX), pressures);
                const double gauss_point_volume = r_integration_points[point_number].Weight() * detJ0;

                for(unsigned int i_node=0; i_node<number_of_nodes; ++i_node) {
                    array_1d<double, 3>& r_gradient = r_geometry[i_node].FastGetSolutionStepValue(PRESSURE_GRADIENT);
                    for(unsigned int k=0; k<num_dim; ++k) {
                        #pragma omp atomic
                        r_gradient[k] += N[i_node] * gauss_point_volume*grad[k];
                    }

                    double& r_vol = r_geometry[i_node].GetValue(NODAL_AREA);

                    #pragma omp atomic
                    r_vol += N[i_node] * gauss_point_volume;
                }
            }
        }
    }

    block_for_each(mrModelPart.Nodes(), [&](Node<3>& rNode){
        if (rNode.GetValue(NODAL_AREA) > 1.0e-12){
            rNode.FastGetSolutionStepValue(PRESSURE_GRADIENT) /= rNode.GetValue(NODAL_AREA);}
    });

    KRATOS_CATCH("")
}

};  // namespace Kratos.
