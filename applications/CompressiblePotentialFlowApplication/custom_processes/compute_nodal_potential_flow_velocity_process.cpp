//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Marc Nunez, based on R.Rossi and V.Mataix work
//
//
//

/* System includes */

/* External includes */

/* Project includes */
#include "utilities/variable_utils.h"
#include "compute_nodal_potential_flow_velocity_process.h"
#include "compressible_potential_flow_application_variables.h"
#include "custom_utilities/potential_flow_utilities.h"


namespace Kratos
{

void ComputeNodalPotentialFlowVelocityProcess::Execute()
{
    KRATOS_TRY;

    // Set to zero
    ClearGradient();

    // Auxiliar containers
    Matrix DN_DX;
    Vector N;

    // First element iterator
    const auto it_element_begin = mrModelPart.ElementsBegin();

    // Current domain size
    const std::size_t dimension = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];
    KRATOS_ERROR_IF(dimension != 2 && dimension !=3) << "Dimension has to be either 2 or 3! Current dimension: " << dimension << std::endl;

    // Iterate over the elements
    #pragma omp parallel for firstprivate(N)
    for(int i_elem=0; i_elem<static_cast<int>(mrModelPart.Elements().size()); ++i_elem) {
        auto it_elem = it_element_begin + i_elem;
        auto& r_geometry = it_elem->GetGeometry();

        // Current geometry information
        const std::size_t number_of_nodes = r_geometry.PointsNumber();

        // Resize if needed
        if (N.size() != number_of_nodes)
            N.resize(number_of_nodes);

        // The integration points
        const auto& r_integration_method = r_geometry.GetDefaultIntegrationMethod();
        const auto& r_integration_points = r_geometry.IntegrationPoints(r_integration_method);
        const std::size_t number_of_integration_points = r_integration_points.size();

        Vector values(number_of_nodes);

        // The containers of the shape functions and the local gradients
        const Matrix& rNmatrix = r_geometry.ShapeFunctionsValues(r_integration_method);

        for ( IndexType i_gauss = 0; i_gauss < number_of_integration_points; ++i_gauss ) {
            // Getting the shape functions
            noalias(N) = row(rNmatrix, i_gauss);

            Vector detJ0;
            GeometryData::ShapeFunctionsGradientsType DN_DX;
            r_geometry.ShapeFunctionsIntegrationPointsGradients(DN_DX, detJ0, r_integration_method);

            const auto r_element = *it_elem;
            Vector velocity_vector(3);
            if (dimension == 2) { // 2D
                velocity_vector = PotentialFlowUtilities::ComputeVelocity<2, 3>(r_element);
            } else if (dimension == 3) { // 3D
                velocity_vector = PotentialFlowUtilities::ComputeVelocity<3, 4>(r_element);
            }
            const Vector grad = velocity_vector;


            const double gauss_point_volume = r_integration_points[i_gauss].Weight() * detJ0[i_gauss];

            for(std::size_t i_node=0; i_node<number_of_nodes; ++i_node) {
                array_1d<double, 3>& r_gradient = GetGradient(r_geometry, i_node);
                for(std::size_t k=0; k<dimension; ++k) {
                    #pragma omp atomic
                    r_gradient[k] += N[i_node] * gauss_point_volume*grad[k];
                }

                double& vol = r_geometry[i_node].GetValue(NODAL_AREA);

                #pragma omp atomic
                vol += N[i_node] * gauss_point_volume;
            }
        }
    }

    PonderateGradient();

    KRATOS_CATCH("")
}

ComputeNodalPotentialFlowVelocityProcess::ComputeNodalPotentialFlowVelocityProcess(
    ModelPart& rModelPart)
    :mrModelPart(rModelPart)
{
    KRATOS_TRY

    // In case the area or gradient variable is not initialized we initialize it
    auto& r_nodes = rModelPart.Nodes();
    if (!r_nodes.begin()->Has(VELOCITY)) {
        const array_1d<double,3> zero_vector = ZeroVector(3);
        VariableUtils().SetNonHistoricalVariable(VELOCITY, zero_vector, r_nodes);
    }
    if (!r_nodes.begin()->Has(NODAL_AREA)) {
        VariableUtils().SetNonHistoricalVariable(NODAL_AREA, 0.0, r_nodes);
    }

    KRATOS_CATCH("")
}

void ComputeNodalPotentialFlowVelocityProcess::ClearGradient()
{
    const array_1d<double, 3> aux_zero_vector = ZeroVector(3);

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i) {
        auto it_node=mrModelPart.NodesBegin()+i;
        it_node->SetValue(NODAL_AREA, 0.0);
        it_node->SetValue(VELOCITY, aux_zero_vector);
    }
}

array_1d<double, 3>& ComputeNodalPotentialFlowVelocityProcess::GetGradient(
    Element::GeometryType& rThisGeometry,
    unsigned int i
    )
{
    array_1d<double, 3>& val = rThisGeometry[i].GetValue(VELOCITY);
    return val;
}

void ComputeNodalPotentialFlowVelocityProcess::PonderateGradient()
{
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(mrModelPart.Nodes().size()); ++i)
    {
        auto it_node=mrModelPart.NodesBegin()+i;
        it_node->GetValue(VELOCITY) /= it_node->GetValue(NODAL_AREA);
    }
}

} /* namespace Kratos.*/
