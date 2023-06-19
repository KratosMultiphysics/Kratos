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
    block_for_each(mrModelPart.Nodes(), [&](Node& rNode){
        rNode.SetValue(NODAL_AREA, 0.0);
        rNode.SetValue(PRESSURE_GRADIENT, zero_vector);
    });

    // Current domain size
    const unsigned int num_dim = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];
    const unsigned int num_nodes = num_dim + 1; //For tetrahedra and triangles

    // Calculate elemental gradient contribution
    if (num_nodes == 3 && num_dim == 2) {
        auto tls_container_2d = SetTLSContainer2D();
        auto elemental_function_2d = GetScalarNodalGradientElementFunction2D();
        block_for_each(mrModelPart.Elements(), tls_container_2d, elemental_function_2d);
    } else if (num_nodes == 4 && num_dim == 3) {
        auto tls_container_3d = SetTLSContainer3D();
        auto elemental_function_3d = GetScalarNodalGradientElementFunction3D();
        block_for_each(mrModelPart.Elements(), tls_container_3d, elemental_function_3d);
    } else {
        KRATOS_ERROR << "Asking for a non-implemented geometry type." << std::endl;
    }

    mrModelPart.GetCommunicator().AssembleNonHistoricalData(PRESSURE_GRADIENT);
    mrModelPart.GetCommunicator().AssembleNonHistoricalData(NODAL_AREA);

    block_for_each(mrModelPart.Nodes(), [&](Node& rNode){
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

CalulateLevelsetConsistentNodalGradientProcess::TLSContainerType2D CalulateLevelsetConsistentNodalGradientProcess::SetTLSContainer2D()
{
    BoundedMatrix<double,3,2> aux_mat;
    array_1d<double,3> aux_vect_1;
    array_1d<double,3> aux_vect_2;
    array_1d<double,3> aux_vect_3;
    array_1d<double,3> aux_vect_4;
    array_1d<double,3> aux_vect_5;
    TLSContainerType2D tls_container_2d = std::make_tuple(aux_mat, aux_vect_1, aux_vect_2, aux_vect_3, aux_vect_4, aux_vect_5);

    return tls_container_2d;
}

CalulateLevelsetConsistentNodalGradientProcess::TLSContainerType3D CalulateLevelsetConsistentNodalGradientProcess::SetTLSContainer3D()
{
    BoundedMatrix<double,4,3> aux_mat;
    array_1d<double,4> aux_vect_1;
    array_1d<double,4> aux_vect_2;
    array_1d<double,4> aux_vect_3;
    array_1d<double,3> aux_vect_4;
    array_1d<double,4> aux_vect_5;
    TLSContainerType3D tls_container_3d = std::make_tuple(aux_mat, aux_vect_1, aux_vect_2, aux_vect_3, aux_vect_4, aux_vect_5);

    return tls_container_3d;
}

std::function<void(Element& rElement, CalulateLevelsetConsistentNodalGradientProcess::TLSContainerType2D& rTLSContainer)> CalulateLevelsetConsistentNodalGradientProcess::GetScalarNodalGradientElementFunction2D()
{
    std::function<void(Element& rElement, CalulateLevelsetConsistentNodalGradientProcess::TLSContainerType2D& rTLSContainer)> aux_func = [&, this](Element& rElement, CalulateLevelsetConsistentNodalGradientProcess::TLSContainerType2D& rTLSContainer){
        this->CalculateScalarNodalGradientElementContribution(rElement, rTLSContainer);
    };

    return aux_func;
}

std::function<void(Element& rElement, CalulateLevelsetConsistentNodalGradientProcess::TLSContainerType3D& rTLSContainer)> CalulateLevelsetConsistentNodalGradientProcess::GetScalarNodalGradientElementFunction3D()
{
    std::function<void(Element& rElement, CalulateLevelsetConsistentNodalGradientProcess::TLSContainerType3D& rTLSContainer)> aux_func = [&, this](Element& rElement, CalulateLevelsetConsistentNodalGradientProcess::TLSContainerType3D& rTLSContainer){
        this->CalculateScalarNodalGradientElementContribution(rElement, rTLSContainer);
    };

    return aux_func;
}

template<class TTLSContainer>
void CalulateLevelsetConsistentNodalGradientProcess::CalculateScalarNodalGradientElementContribution(
    Element& rElement,
    TTLSContainer& rTLSContainer)
{
    // Get auxiliary arrays from the TLS container
    auto& r_DN_DX = std::get<0>(rTLSContainer);
    auto& r_N = std::get<1>(rTLSContainer);
    auto& r_pressures = std::get<2>(rTLSContainer);
    auto& r_distances = std::get<3>(rTLSContainer);
    auto& r_grad = std::get<4>(rTLSContainer);
    auto& r_midpoint_N = std::get<5>(rTLSContainer);

    // Current geometry information
    auto& r_geometry = rElement.GetGeometry();
    const unsigned number_of_nodes = r_geometry.PointsNumber();

    for(std::size_t i_node=0; i_node<number_of_nodes; ++i_node){
        r_distances(i_node) = r_geometry[i_node].FastGetSolutionStepValue(DISTANCE);
    }

    if(!IsSplit(r_distances))
    {
        // The integration points
        const auto& r_integration_method = r_geometry.GetDefaultIntegrationMethod();
        const auto& r_integration_points = r_geometry.IntegrationPoints(r_integration_method);
        const unsigned int number_of_integration_points = r_integration_points.size();

        // Set nodal scalar variable values to calculate the gradient
        for(std::size_t i_node=0; i_node<number_of_nodes; ++i_node) {
            r_pressures[i_node] = r_geometry[i_node].FastGetSolutionStepValue(PRESSURE);
        }

        // The containers of the shape functions and the local gradients
        double volume;
        const auto& rNcontainer = r_geometry.ShapeFunctionsValues(r_integration_method);
        GeometryUtils::CalculateGeometryData(r_geometry, r_DN_DX, r_midpoint_N, volume);

        for ( IndexType point_number = 0; point_number < number_of_integration_points; ++point_number ) {
            // Getting the shape functions and calculate values
            noalias(r_N) = row(rNcontainer, point_number);
            noalias(r_grad) = ZeroVector(3);
            for (unsigned int i = 0; i < r_DN_DX.size1(); ++i) {
                for (unsigned int j = 0; j < r_DN_DX.size2(); ++j) {
                    r_grad[j] += r_DN_DX(i,j) * r_pressures(i);
                }
            }
            const double gauss_point_volume = r_integration_points[point_number].Weight() * volume;
            // Atomic addition of the gradient and the weight
            for(unsigned int i_node=0; i_node<number_of_nodes; ++i_node) {
                auto& r_gradient = r_geometry[i_node].GetValue(PRESSURE_GRADIENT);
                AtomicAddVector(r_gradient, r_N[i_node]*gauss_point_volume*r_grad);

                double& r_vol = r_geometry[i_node].GetValue(NODAL_AREA);
                AtomicAdd(r_vol, r_N[i_node] * gauss_point_volume);
            }
        }
    }
}

};  // namespace Kratos.
