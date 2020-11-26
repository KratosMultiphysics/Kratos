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
#include "utilities/atomic_utilities.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "calulate_levelset_consistent_nodal_gradient_process.h"


namespace Kratos
{

/* Public functions *******************************************************/

/// constructor
CalulateLevelsetConsistentNodalGradientProcess::CalulateLevelsetConsistentNodalGradientProcess(
        ModelPart& rModelPart)
    : Process(), mrModelPart(rModelPart) {}

/// Constructor with Kratos parameters.
CalulateLevelsetConsistentNodalGradientProcess::CalulateLevelsetConsistentNodalGradientProcess(
    ModelPart& rModelPart,
    Parameters Parameters)
    : CalulateLevelsetConsistentNodalGradientProcess(
        rModelPart
    ){}

/// Constructor with Kratos model
CalulateLevelsetConsistentNodalGradientProcess::CalulateLevelsetConsistentNodalGradientProcess(
    Model& rModel,
    Parameters Parameters)
    : CalulateLevelsetConsistentNodalGradientProcess(
        rModel.GetModelPart(Parameters["model_part_name"].GetString())
    ){}

void CalulateLevelsetConsistentNodalGradientProcess::Execute(){

    KRATOS_TRY;

    const auto zero_vector = ZeroVector(3);
    // Set to zero
    block_for_each(mrModelPart.Nodes(), [&](Node<3>& rNode){
        rNode.SetValue(NODAL_AREA, 0.0);
        rNode.SetValue(PRESSURE_GRADIENT, zero_vector);
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

    const int elements_number = mrModelPart.Elements().size();

    // Iterate over the elements
    #pragma omp parallel for firstprivate(DN_DX, N, J0, InvJ0, detJ0, pressures, distances, grad)
    for(int i_elem = 0; i_elem < elements_number; ++i_elem) {

        const auto it_elem = it_element_begin + i_elem;
        auto& r_geometry = it_elem->GetGeometry();

        // Current geometry information
        const unsigned number_of_nodes = r_geometry.PointsNumber();

        for(std::size_t i_node=0; i_node<number_of_nodes; ++i_node){
            distances(i_node) = r_geometry[i_node].FastGetSolutionStepValue(DISTANCE);
        }

        if(!IsSplit(distances))
        {
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
                    auto& r_gradient = r_geometry[i_node].GetValue(PRESSURE_GRADIENT);
                    AtomicAdd(r_gradient, N[i_node]*gauss_point_volume*grad);

                    double& r_vol = r_geometry[i_node].GetValue(NODAL_AREA);
                    AtomicAdd(r_vol, N[i_node] * gauss_point_volume);
                }
            }
        }
    }

    block_for_each(mrModelPart.Nodes(), [&](Node<3>& rNode){
        if (rNode.GetValue(NODAL_AREA) > 1.0e-12){
            rNode.GetValue(PRESSURE_GRADIENT) /= rNode.GetValue(NODAL_AREA);}
    });

    KRATOS_CATCH("")
}

bool CalulateLevelsetConsistentNodalGradientProcess::IsSplit(const Vector& rDistances)
{
    bool is_split = false;

    unsigned int nneg=0, npos=0;
    for(unsigned int i = 0; i < rDistances.size(); ++i)
    {
        if(rDistances[i] > 0) {
            npos += 1;
        } else {
            nneg += 1;
        }
    }

    if(nneg > 0 && npos > 0)
        is_split = true;

    return is_split;
}

const Parameters CalulateLevelsetConsistentNodalGradientProcess::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
    {
        "model_part_name"             : "please_specify_model_part_name"
    })" );

    return default_parameters;
}

};  // namespace Kratos.
