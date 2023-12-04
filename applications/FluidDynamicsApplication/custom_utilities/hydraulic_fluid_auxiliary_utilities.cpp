//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Suneth Warnakulasuriya
//

// System includes


// External includes


// Project includes
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"
#include "processes/find_global_nodal_neighbours_process.h"
#include "includes/global_pointer_variables.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "spatial_containers/bins_dynamic.h"
#include "utilities/rbf_shape_functions_utility.h"

// Application includes
#include "fluid_auxiliary_utilities.h"
#include "hydraulic_fluid_auxiliary_utilities.h"

namespace Kratos
{
double HydraulicFluidAuxiliaryUtilities::CalculateWettedPetimeter(
    ModelPart &rModelPart,
    const Flags &rSkinFlag)
{
   
    // Check that there are conditions and distance_inlet variable in the nodal database
    KRATOS_ERROR_IF(rModelPart.GetProcessInfo()[DOMAIN_SIZE] < 3 )<< "Wetted perimeter is only implemented for 3D." << std::endl;
    KRATOS_ERROR_IF(rModelPart.GetCommunicator().GlobalNumberOfConditions() == 0) << "There are no conditions in the provided model part. Wetted perimeter cannot be computed." << std::endl;
    const auto &r_communicator = rModelPart.GetCommunicator();
    double wedded_perimeter = 0.0;

    // if (r_communicator.LocalMesh().NumberOfNodes() != 0)

    // {
    //     KRATOS_ERROR_IF_NOT(r_communicator.LocalMesh().NodesBegin()->SolutionStepsDataHas(AUX_DISTANCE)) << "Nodal solution step data has no \'AUX_DISTANCE\' variable. Wetted perimeter cannot be computed" << std::endl;
    //     // Create a vector where all the face edges id are stored. 
    //     std::vector<std::string> edge_boundary_id;
    //     for (auto& cond : r_communicator.LocalMesh().Conditions())
    //     {
    //         // if (cond.Is(rSkinFlag))
    //         // {
    //         //     const auto& r_geom = cond.GetGeometry();
    //         //     const std::size_t n_nodes = r_geom.PointsNumber();
    //         //     for (std::size_t i_node = 0; i_node < n_nodes; ++i_node){
    //         //         for (auto &r_neigh : r_geom[i_node].GetValue(NEIGHBOUR_NODES)){
    //         //             // if (r_neigh.Is(rSkinFlag)){
    //         //             //     // std::string edge_id;
    //         //             //     // // if (r_geom[i_node].Id()<r_neigh.Id()){
    //         //             //     // //     std::string edge_id = std::to_string(r_geom[i_node].Id()) + std::to_string(r_neigh.Id());
    //         //             //     // // }
    //         //             //     // // else{
    //         //             //     // //     std::string edge_id = std::to_string(r_neigh.Id()) + std::to_string(r_geom[i_node].Id());
    //         //             //     // // }

    //         //             //     // edge_boundary_id.push_back(edge_id);
    //         //             // }
    //         //             r_neigh.Id();
    //         //         }
    //         //     }
    //         }
    //     }

    //     // Create a maps structure with edge_id and the frecuency of repetition 
    //     std::map<std::string, int> edge_frecuency;
    //     for (const auto &edge : edge_boundary_id)
    //     {
    //         edge_frecuency[edge]++;
    //     }

    //     // Delete from the edge_id container all the edges which are repeated since they do not belong to the perimeter. 
    //     auto it = edge_boundary_id.begin();
    //     while (it != edge_boundary_id.end())
    //     {
    //         if (edge_frecuency[*it] > 1)
    //         {
    //             it = edge_boundary_id.erase(it);
    //         }
    //         else
    //         {
    //             ++it; 
    //         }
    //     }

    // }

return wedded_perimeter;
}



double HydraulicFluidAuxiliaryUtilities::CalculateConditionArea(const GeometryType &rGeometry)
{
    // Calculate the condition area normal
    GeometryType::CoordinatesArrayType point_local;
    rGeometry.PointLocalCoordinates(point_local, rGeometry.Center()) ;
    const array_1d<double,3> area_normal = rGeometry.Normal(point_local);
    // Check condition area and calculate the condition average flow rate
    double condition_area = 0.0;
    if (norm_2(area_normal) > std::numeric_limits<double>::epsilon()) {
            condition_area = norm_2(area_normal);
   
    } else {
        KRATOS_WARNING("CalculateConditionArea") << "Condition area is close to zero." << std::endl;
    }
    return condition_area;
}
// bool HydraulicFluidAuxiliaryUtilities::CheckNegativeSubdomain(const GeometryType){

// }
// bool HydraulicFluidAuxiliaryUtilities::CheckSplitSubdomain(const GeometryType){

// }
double HydraulicFluidAuxiliaryUtilities::CalculateWettedArea(ModelPart &rModelPart, const Flags &rSkinFlag)
{
//     // Check that there are conditions and distance variable in the nodal database
//     KRATOS_ERROR_IF(rModelPart.GetCommunicator().GlobalNumberOfConditions() == 0) << "There are no conditions in the provided model part. Flow rate cannot be computed." << std::endl;
//     const auto& r_communicator = rModelPart.GetCommunicator();
//     if (r_communicator.LocalMesh().NumberOfNodes() !=0) {
//         KRATOS_ERROR_IF_NOT(r_communicator.LocalMesh().NodesBegin()->SolutionStepsDataHas(AUX_DISTANCE)) << "Nodal solution step data has no \'AUX_DISTANCE\' variable. Flow rate cannot be computed" << std::endl;
//     }
    double wetted_area = 0.0;
//     if (r_communicator.LocalMesh().NumberOfConditions() != 0) {
//         // Create the modified shape functions factory with the first condition parent as prototype
//         const auto& r_cond_begin = r_communicator.LocalMesh().ConditionsBegin();
//         const auto& r_parent_cond_begin = r_cond_begin->GetValue(NEIGHBOUR_ELEMENTS)[0];
//         auto mod_sh_func_factory = FluidAuxiliaryUtilities::GetStandardModifiedShapeFunctionsFactory(r_parent_cond_begin.GetGeometry());

//         std::size_t n_dim = rModelPart.GetProcessInfo().GetValue(DOMAIN_SIZE);
//         Vector nodal_distances(r_cond_begin->GetGeometry().PointsNumber());
//         wetted_area = block_for_each<SumReduction<double>>(r_communicator.LocalMesh().Conditions(), nodal_distances, [&](Condition& rCondition, Vector& rNodalDistancesTLS){
//                 // Check if the condition is to be added to the flow contribution
//                 double cond_wetted_area = 0.0;
//                 if (rCondition.Is(rSkinFlag)) {
//                     // Get geometry data
//                     const auto& r_geom = rCondition.GetGeometry();
//                     const std::size_t n_nodes = r_geom.PointsNumber();

//                     // Set up distances vector
//                     for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
//                         rNodalDistancesTLS(i_node) = r_geom[i_node].FastGetSolutionStepValue(AUX_DISTANCE);
//                     }
//                     // Check if the condition is in the negative subdomain or intersected
//                     if (CheckNegativeSubdomain(rNodalDistancesTLS))
//                     {
//                         cond_flow_rate = CalculateConditionArea(r_geom);
//                     }
//                     else if (CheckSplitSubdomain(rNodalDistancesTLS))
//                     {
//                         // Get the current condition parent
//                         const auto p_parent_element = rCondition.GetValue(NEIGHBOUR_ELEMENTS)(0);
//                         KRATOS_ERROR_IF_NOT(p_parent_element.get()) << "Condition " << rCondition.Id() << " has no parent element assigned." << std::endl;
//                         const auto& r_parent_geom = p_parent_element->GetGeometry();

//                         // Get the corresponding face id of the current condition
//                         const std::size_t n_parent_faces = r_parent_geom.FacesNumber();
//                         DenseMatrix<unsigned int> nodes_in_faces(n_parent_faces, n_parent_faces);
//                         r_parent_geom.NodesInFaces(nodes_in_faces);
//                         std::size_t face_id = 0;
//                         for (std::size_t i_face = 0; i_face < n_parent_faces; ++i_face) {
//                             std::size_t match_nodes = 0;
//                             for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
//                                 std::size_t parent_local_id = nodes_in_faces(i_node + 1, i_face);
//                                 for (std::size_t j_node = 0; j_node < n_nodes; ++j_node) {
//                                     if (r_geom[j_node].Id() == r_parent_geom[parent_local_id].Id()) {
//                                         match_nodes++;
//                                         break;
//                                     }
//                                 }
//                             }
//                             if (match_nodes == n_nodes) {
//                                 face_id = i_face;
//                                 break;
//                             }
//                         }

//                         // Calculate the modified shape functions in the face of interest
//                         const std::size_t n_nodes_parent = r_parent_geom.PointsNumber();
//                         Vector parent_distances(n_nodes_parent);
//                         for (std::size_t i_node = 0; i_node < n_nodes_parent; ++i_node) {
//                             parent_distances(i_node) = r_parent_geom[i_node].FastGetSolutionStepValue(AUX_DISTANCE);
//                         }
//                         auto p_mod_sh_func = mod_sh_func_factory(p_parent_element->pGetGeometry(), parent_distances);

//                         //                 Vector w_vect;
//                         //                 Matrix N_container;
//                         //                 std::vector<array_1d<double,3>> normals_vect;
//                         //                 CalculateSplitConditionGeometryData<IsPositiveSubdomain>(p_mod_sh_func, face_id, N_container, normals_vect, w_vect);

//                         //                 // Interpolate the flow rate in the positive subdomain
//                         //                 Vector i_normal(n_dim);
//                         //                 Vector i_N(n_nodes_parent);
//                         //                 array_1d<double,3> aux_vel;
//                         //                 const std::size_t n_gauss = w_vect.size();
//                         //                 for (std::size_t i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
//                         //                     aux_vel = ZeroVector(3);
//                         //                     i_N = row(N_container, i_gauss);
//                         //                     i_normal = normals_vect[i_gauss];
//                         //                     const double i_normal_norm = norm_2(i_normal);
//                         //                     KRATOS_WARNING_IF("CalculateFlowRatePositiveSkin", i_normal_norm < 1.0e-12) << "Condition " << rCondition.Id() << " normal close to zero." << std::endl;
//                         //                     i_normal /= i_normal_norm;
//                         //                     for (std::size_t i_parent_node = 0; i_parent_node < n_nodes_parent; ++i_parent_node) {
//                         //                         noalias(aux_vel) += i_N[i_parent_node] * r_parent_geom[i_parent_node].FastGetSolutionStepValue(VELOCITY);
//                         //                     }
//                         //                     cond_flow_rate += w_vect[i_gauss] * inner_prod(aux_vel, i_normal);
//                         //                 }
//                     }
//                 }

//                 return cond_wetted_area;
//             });
//     }

//     // Synchronize among processors
//     wetted_area = r_communicator.GetDataCommunicator().SumAll(wetted_area);

    return wetted_area;
}

double HydraulicFluidAuxiliaryUtilities::CalculateInletDischarge(ModelPart &rModelPart, const Flags &rSkinFlag)
{
    // Assume that the water depth in the inlet condition is going to be the critical depth. Fr=1.0
    double water_area = CalculateWettedArea(rModelPart, rSkinFlag);
    double wetted_perimeter = CalculateWettedPetimeter(rModelPart, rSkinFlag);
    const double gravity = 9.81;
    double inlet_discharge = std::sqrt(gravity * water_area / wetted_perimeter) * water_area;

    return inlet_discharge;
}
} // namespace Kratos
