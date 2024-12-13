//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Uxue Chasco
//

// System includes


// External includes


// Project includes
#include "processes/find_global_nodal_neighbours_process.h"
#include "includes/global_pointer_variables.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "spatial_containers/bins_dynamic.h"
#include "utilities/rbf_shape_functions_utility.h"
#include "utilities/divide_triangle_3d_3.h"
#include "../../FluidDynamicsApplication/custom_utilities/fluid_auxiliary_utilities.h"
#include "../../FluidDynamicsApplication/fluid_dynamics_application_variables.h"

// Application includes
#include "hydraulic_fluid_auxiliary_utilities.h"

namespace Kratos
{
double HydraulicFluidAuxiliaryUtilities::CalculateWettedPetimeter(
    ModelPart &rModelPart,
    const Flags &rSkinFlag,
    const Variable<double>& rDistanceVariable,
    const bool IsHistorical)
{
    // Auxiliary function to have the posssibility of having non historical or hitorical inlet distance function.
    std::function<double(NodeType&, const Variable<double>&)> distance_getter;
    if (IsHistorical)
    {
        distance_getter = [](NodeType& rNode, const Variable<double>& rDistanceVariable) -> double {return rNode.FastGetSolutionStepValue(rDistanceVariable);};
    } else {
        distance_getter = [](NodeType& rNode, const Variable<double>& rDistanceVariable) -> double {return rNode.GetValue(rDistanceVariable);};
    }

    // Check that there are conditions and distance_inlet variable in the nodal database
    KRATOS_ERROR_IF(rModelPart.GetProcessInfo()[DOMAIN_SIZE] < 3 )<< "Wetted perimeter is only implemented for 3D." << std::endl;

    const auto &r_communicator = rModelPart.GetCommunicator();
    double wedded_perimeter = 0.0;
    if (r_communicator.LocalMesh().NumberOfNodes() != 0)
    {
        KRATOS_ERROR_IF(IsHistorical && !r_communicator.LocalMesh().NodesBegin()->SolutionStepsDataHas(AUX_DISTANCE)) << "Nodal solution step data has no \'AUX_DISTANCE\' variable. Wetted perimeter cannot be computed" << std::endl;

        // Create a vector where all the face edges id are stored.
        std::pair<int,int> aux_pair;
        std::unordered_map<std::pair<int, int>, EdgeDataContainer> edges_map;

        for (auto& r_cond : rModelPart.Conditions())
        {
            if (r_cond.Is(rSkinFlag))
            {
                auto& r_geom = r_cond.GetGeometry();
                const SizeType n_nodes = r_geom.PointsNumber();
                KRATOS_ERROR_IF(n_nodes != 3) << "This function only supports Triangle3D3N geometries." << std::endl;

                // Loop through pairs of nodes to identify unique edges
                for (IndexType i = 0; i < n_nodes - 1; ++i)
                {
                    const IndexType i_id = r_geom[i].Id();
                    for (IndexType j = i + 1; j < n_nodes; ++j)
                    {
                        const IndexType j_id = r_geom[j].Id();

                        // Ensure consistent ordering of node pairs for uniqueness
                        if (i_id < j_id) {
                            aux_pair = std::make_pair<int, int>(i_id, j_id);
                        } else{
                            aux_pair = std::make_pair<int, int>(j_id, i_id);
                        }

                        // Check if the edge is already in the map
                        auto found = edges_map.find(aux_pair);
                        // If not, create a new entry in the map
                        if (found == edges_map.end())
                        {
                            EdgeDataContainer edge_data;
                            edge_data.pNodeI = i_id < j_id ? r_geom(i) : r_geom(j);
                            edge_data.pNodeJ = i_id < j_id ? r_geom(j) : r_geom(i);
                            edge_data.NumberOfRepetitions = 1;
                            edges_map.insert(std::make_pair(aux_pair, edge_data));
                        }
                        // If the edge is already in the map, update the repetition count
                        else
                        {
                            auto &r_edge_data = found->second;
                            r_edge_data.NumberOfRepetitions += 1;
                        }
                    }
                }
            }
        }

        for (auto it = edges_map.begin(); it != edges_map.end();)
        {   // Delete all edges that have more than one repetition, since they are not part of the perimeter
            const auto& r_edge_data = it->second;
            if (r_edge_data.NumberOfRepetitions > 1) {
                it = edges_map.erase(it);
            } else {
                ++it;
            }
        }
        // Calculate de distance of each edge belonging to the perimeter.
        array_1d<double, 3> edge_vector;
        for (auto it = edges_map.begin(); it != edges_map.end(); ++it)
        {
            const auto& r_edge_data = it->second;
            const double distance_value_j = distance_getter(*(r_edge_data.pNodeJ), rDistanceVariable);
            const double distance_value_i = distance_getter(*(r_edge_data.pNodeI), rDistanceVariable);
            // Case 1: The edge is completly wetted
            if ((distance_value_i<0.0) && (distance_value_j<0.0) ){
                edge_vector = r_edge_data.pNodeJ->Coordinates() - r_edge_data.pNodeI->Coordinates();
                const double edge_length = norm_2(edge_vector);
                wedded_perimeter += edge_length;
            }
            // Case 2: The edge is cut. Interpolate the wetted distance.
            else if (distance_value_i * distance_value_j<0){

                const double phi_neg = distance_value_i<0? distance_value_i: distance_value_j;
                const double phi_pos = distance_value_i>0? distance_value_i: distance_value_j;
                const double distance_int = std::abs(phi_neg)/ (std::abs(phi_neg)+  std::abs(phi_pos));
                edge_vector = r_edge_data.pNodeJ->Coordinates() - r_edge_data.pNodeI->Coordinates();
                const double edge_length = norm_2(edge_vector);
                wedded_perimeter += distance_int*edge_length;
            }
        }
    }

    return wedded_perimeter;
}

double HydraulicFluidAuxiliaryUtilities::CalculateWettedArea(
    ModelPart &rModelPart,
    const Flags &rSkinFlag,
    const Variable<double> &rDistanceVariable,
    bool IsHistorical)
{
    // Auxiliary function to have the posssibility of having non historical or hitorical inlet distance function.
    std::function<double(NodeType &, const Variable<double> &)> distance_getter;
    if (IsHistorical)
    {
        distance_getter = [](NodeType &rNode, const Variable<double> &rDistanceVariable) -> double
        { return rNode.FastGetSolutionStepValue(rDistanceVariable); };
    }
    else
    {
        distance_getter = [](NodeType &rNode, const Variable<double> &rDistanceVariable) -> double
        { return rNode.GetValue(rDistanceVariable); };
    }

    KRATOS_ERROR_IF(rModelPart.GetProcessInfo()[DOMAIN_SIZE] < 3) << "Wetted perimeter is only implemented for 3D." << std::endl;

    // Auxiliary container for the fake Triangle2D3 geometries' points
    GeometryType::PointsArrayType aux_points;
    double wetted_area = 0.0;

    for (auto &r_cond : rModelPart.Conditions())
    {   // Create an auxiliary element based on the rskinflag condition.
        if (r_cond.Is(rSkinFlag))
        {
            std::vector<ModelPart::IndexType> elem_nodes_id;
            auto &r_geom = r_cond.GetGeometry();


            // Fill the points array and distances vector
            Vector aux_distances(3);
            for (IndexType i_nodes = 0; i_nodes < r_geom.PointsNumber(); i_nodes++)
            {
                aux_points.push_back(r_geom(i_nodes));
                aux_distances[i_nodes] = r_geom[i_nodes].GetValue(rDistanceVariable);
                // TODO: It should be posible to have an historical distance variable.
                // aux_distances[i_nodes] =distance_getter;
            }
            // Calculate the water area (wetted area) of the cut conditions
            if (FluidAuxiliaryUtilities::IsSplit(aux_distances))
            {
                auto p_splitting_util = Kratos::make_unique<DivideTriangle3D3<NodeType>>(r_geom, aux_distances);
                p_splitting_util->GenerateDivision();
                const auto& r_neg_subdivisions = p_splitting_util->GetNegativeSubdivisions();
                for (const auto& rp_neg_subdivision : r_neg_subdivisions) {
                    wetted_area += rp_neg_subdivision->Area();
                }
            }
            // Calculate the water area (wetted area)
            else if (FluidAuxiliaryUtilities::IsNegative(aux_distances))
            {
                wetted_area += r_geom.Area();
            }
            aux_points.clear();
        }
    }
    return wetted_area;
}

double HydraulicFluidAuxiliaryUtilities::InitialWaterDepth(ModelPart &rModelPart)
{

   //Determine the initial estimate for water depth by considering the average of the maximum and minimum coordinates.
    KRATOS_INFO("HydraulicFluidAuxiliaryUtilities") << "The water depth is assumed to be in the Z direction" << std::endl;
    // double max_value = std::numeric_limits<double>::lowest();
    // double min_value = std::numeric_limits<double>::max();
    // for (auto &r_node : rModelPart.Nodes())
    // {
    //     if (r_node.Z() > max_value)
    //     {
    //         max_value = r_node.Z();
    //     }
    //     if (r_node.Z() < min_value)
    //     {
    //         min_value = r_node.Z();
    //     }
    // }

    double min_value, max_value;
    std::tie(min_value, max_value) = block_for_each<CombinedReduction<MinReduction<double>,MaxReduction<double>>>(rModelPart.Nodes(), [](NodeType& rNode){
        return std::make_tuple(rNode.Z(), rNode.Z());
    });

    return 0.5 * std::abs(max_value - min_value);
}

void HydraulicFluidAuxiliaryUtilities::SetInletVelocity(
    ModelPart &rModelPart,
    double InletVelocity,
    const Variable<double> &rDistanceVariable)
{
    struct AuxTLS
    {
        array_1d<double,3> InletNorm;
        array_1d<double,3> InletVelocity;
    };

    block_for_each(rModelPart.Nodes(), AuxTLS(), [&](NodeType &rNode, AuxTLS &rTLS)
    {
        // Get TLS variables
        auto& inlet_norm = rTLS.InletNorm;
        auto& inlet_velocity = rTLS.InletVelocity;

        inlet_norm = rNode.GetValue(INLET_NORMAL);

        const double n_norm = norm_2(inlet_norm);
        // Orient the velocity vector in the opposite direction of the outer normal to ensure it is directed inwards.

        if (n_norm > 1.0e-12)
        {
            inlet_norm /= -n_norm;
        }
        else
        {
            KRATOS_WARNING("SetInletVelocity") << "Node " << rNode.Id() << " INLET_NORMAL is close to zero." << std::endl;
            inlet_norm /= -1.0;
        }
        //  Inlet velocity vector.
        inlet_velocity = inlet_norm * InletVelocity;

        //  Assign the velocity vector to each node representing water in the inlet condition and fix its value
        if (rNode.GetValue(rDistanceVariable) < 0.0)
        {
            rNode.FastGetSolutionStepValue(VELOCITY) =  inlet_velocity;
            rNode.Fix(VELOCITY_X);
            rNode.Fix(VELOCITY_Y);
            rNode.Fix(VELOCITY_Z);
        }
        else{
            // The air velocity in the inlet node is assumed to be null.
            rNode.FastGetSolutionStepValue(VELOCITY_X) = 0.0;
            rNode.FastGetSolutionStepValue(VELOCITY_Y) = 0.0;
            rNode.FastGetSolutionStepValue(VELOCITY_Z) = 0.0;
            rNode.Fix(VELOCITY_X);
            rNode.Fix(VELOCITY_Y);
            rNode.Fix(VELOCITY_Z);
            // rNode.FastGetSolutionStepValue(VELOCITY) = inlet_velocity;
            // rNode.Fix(VELOCITY_X);
            // rNode.Fix(VELOCITY_Y);
            // rNode.Fix(VELOCITY_Z);
        }
    });
}
void HydraulicFluidAuxiliaryUtilities::FreeInlet(ModelPart& rModelPart)
{
    // Free the velocity.
    block_for_each(rModelPart.Nodes(), [](NodeType &rNode){
        rNode.Free(VELOCITY_X);
        rNode.Free(VELOCITY_Y);
        rNode.Free(VELOCITY_Z);
        rNode.Free(DISTANCE);
    });
}
void HydraulicFluidAuxiliaryUtilities::SetInletFreeSurface(ModelPart &rModelPart, const Flags &rSkinFlag,  const Variable<double> &rDistanceVariable)
{
    // Assign the water depth (DISTANCE) to all nodes within the inlet model part to be equal to the water depth corresponding to a Froude number of 1 (&rDistanceVariable).
    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
        if (rNode.Is(rSkinFlag)){
            double inlet_dist = rNode.GetValue(rDistanceVariable);
            rNode.FastGetSolutionStepValue(DISTANCE) = inlet_dist;
            if (inlet_dist<0.0){
                rNode.Fix(DISTANCE);
            }
        }
    });
}
void HydraulicFluidAuxiliaryUtilities::TurnOffGravityOnAirElements(ModelPart &rModelPart){
    block_for_each(rModelPart.Elements(), [&](Element &rElement){
        auto &r_geom = rElement.GetGeometry();
        Vector distances(3);
        Vector gravity = ZeroVector(3);
        for (IndexType i_nodes = 0; i_nodes < r_geom.PointsNumber(); i_nodes++){
            distances[i_nodes] = r_geom[i_nodes].FastGetSolutionStepValue(DISTANCE);
        }
        if (FluidAuxiliaryUtilities::IsPositive(distances)){
            for (IndexType i_nodes = 0; i_nodes < r_geom.PointsNumber(); i_nodes++){
                r_geom[i_nodes].FastGetSolutionStepValue(BODY_FORCE) =gravity;
            }
        }
    });
}
// void HydraulicFluidAuxiliaryUtilities::FixCornerNodeVelocity(ModelPart &rModelPart, double angle_corner)
// {
//     // Obtain for each condition each neighbor condition
//         block_for_each(rModelPart.Nodes(), [&](NodeType &rNode){
//             rNode.SetValue(AUX_INDEX,0);
//         });
//     ModelPart::ConditionsContainerType &r_cond = rModelPart.Conditions();
//     double angle_corner_rad = angle_corner / 180.0 * 3.14159; //
//     double acceptable_cos = cos(angle_corner_rad);

//     for (ModelPart::ConditionsContainerType::iterator cond_it = r_cond.begin(); cond_it != r_cond.end(); cond_it++)
//     {

//         // get geometry data of the face
//         Geometry<Node> &face_geometry = cond_it->GetGeometry();
//         // reference for area normal of the face
//         const array_1d<double, 3> &face_normal = cond_it->GetValue(NORMAL);
//         double An = norm_2(face_normal);
//         const GlobalPointersVector<Condition> &neighb = cond_it->GetValue(NEIGHBOUR_CONDITIONS);
//         KRATOS_WATCH(neighb.size())
//         for (unsigned int c_itr = 0; c_itr < neighb.size(); c_itr++)
//         {
//             Condition rCondition = neighb[c_itr];

//             const array_1d<double, 3> &face_neig_normal = rCondition.GetValue(NORMAL);
//             double cos_angle = 1 / (An * norm_2(face_neig_normal)) * inner_prod(face_normal,face_neig_normal);
//             if (cos_angle <acceptable_cos)
//             {
//                 //it is considered corner angle
//                 for (auto &node:face_geometry){
//                     double repeated = node.GetValue(AUX_INDEX);
//                     repeated+=1;
//                     node.SetValue(AUX_INDEX, repeated);
//                 }
//             }
//         }
//      }
//      block_for_each(rModelPart.Nodes(), [&](NodeType &rNode){
//         if (rNode.GetValue(AUX_INDEX)>=6){

//             // std::cout << rNode.GetValue(AUX_INDEX) << " " << rNode.Id() << std::endl;
//             Vector veloctiy = ZeroVector(3);
//             rNode.FastGetSolutionStepValue(VELOCITY) = veloctiy;
//             rNode.Fix(VELOCITY_X);
//             rNode.Fix(VELOCITY_Y);
//             rNode.Fix(VELOCITY_Z);

//         } });
// }

void HydraulicFluidAuxiliaryUtilities::FixCornerNodeVelocity(ModelPart &rModelPart, double angle_corner)
{
    // Obtain for each condition each neighbor condition
    block_for_each(rModelPart.Nodes(), [&](NodeType &rNode)
                   { rNode.SetValue(AUX_INDEX, 0); });
    ModelPart::ConditionsContainerType &r_cond = rModelPart.Conditions();
    double angle_corner_rad = angle_corner / 180.0 * 3.14159; //
    double acceptable_cos = cos(angle_corner_rad);
    double max_height = block_for_each<MaxReduction<double>>(rModelPart.Nodes(), [&](NodeType &rNode){
       return rNode.Z();
    });
    block_for_each(rModelPart.Nodes(),[&](NodeType &rNode){
       if(rNode.Z() >=max_height){
         rNode.Set(MARKER,true);
       }
    });

    for (ModelPart::ConditionsContainerType::iterator cond_it = r_cond.begin(); cond_it != r_cond.end(); cond_it++)
    {


        // auto r_geom = cond_it->GetGeometry();
        // auto edgelist = r_geom.GenerateEdges();
        // get geometry data of the face
        // Geometry<Node> &face_geometry = cond_it->GetGeometry();
        // reference for area normal of the face
        const array_1d<double, 3> &face_normal = cond_it->GetValue(NORMAL);
        double An = norm_2(face_normal);
        const GlobalPointersVector<Condition> &neighb = cond_it->GetValue(NEIGHBOUR_CONDITIONS);
        const GeometryType &r_geom = cond_it->GetGeometry();
        auto edgelist = r_geom.GenerateEdges();


        // auto r_geom = cond_it->GetGeometry();
        // auto edgelist = r_geom.GenerateEdges();
        // // get geometry data of the face
        // // // Geometry<Node> &face_geometry = cond_it->GetGeometry();

        // // // reference for area normal of the face
        // // const array_1d<double, 3> &face_normal = cond_it->GetValue(NORMAL);
        // // double An = norm_2(face_normal);
        // // const GlobalPointersVector<Condition> &neighb = cond_it->GetValue(NEIGHBOUR_CONDITIONS);

        // KRATOS_WATCH(neighb.size())
        for (unsigned int c_itr = 0; c_itr < neighb.size(); c_itr++)
        {
            auto rCondition = neighb[c_itr];
            const array_1d<double, 3> &face_neig_normal = rCondition.GetValue(NORMAL);
            double cos_angle = 1 / (An * norm_2(face_neig_normal)) * inner_prod(face_normal, face_neig_normal);
            auto &r_geom_neig = rCondition.GetGeometry();
            auto edgelist_neig = r_geom_neig.GenerateEdges();

            array_1d<double, 2> ids_neig;
            array_1d<double, 2> ids_cond;
            if (cos_angle < acceptable_cos)
            {
                for (IndexType edge = 0; edge < edgelist.size(); edge++)
                {
                    for (IndexType i = 0; i < edgelist[edge].size(); i++)
                    {
                        ids_cond[i] = edgelist[edge][i].Id();
                    }
                    std::sort(ids_cond.begin(), ids_cond.end());

                    for (IndexType edge_neg = 0; edge_neg < edgelist_neig.size(); edge_neg++)
                    {
                        for (IndexType i = 0; i < edgelist_neig[edge_neg].size(); i++)
                        {
                            ids_neig[i] = edgelist_neig[edge_neg][i].Id();
                        }
                        std::sort(ids_neig.begin(), ids_neig.end());

                        if (ids_cond == ids_neig)
                        {

                            for (IndexType i = 0; i < edgelist[edge].size(); i++){

                                double repeated = edgelist[edge][i].GetValue(AUX_INDEX);
                                repeated += 1;
                                edgelist[edge][i].SetValue(AUX_INDEX, repeated);
                                if (edgelist[edge][i].Id()==1031)
                                {
                                    auto geom = cond_it->GetGeometry();
                                }

                            }
                        }
                    }
                }
            }
        }
    }
    block_for_each(rModelPart.Nodes(),[&](NodeType &rNode){
        if (rNode.GetValue(AUX_INDEX)*0.5>=3 && rNode.Is(MARKER)){

            // std::cout << rNode.GetValue(AUX_INDEX) << " " << rNode.Id() << std::endl;
            std::cout  << rNode.Id() << std::endl;
            Vector veloctiy = ZeroVector(3);
            rNode.FastGetSolutionStepValue(VELOCITY) = veloctiy;
            rNode.Fix(VELOCITY_X);
            rNode.Fix(VELOCITY_Y);
            rNode.Fix(VELOCITY_Z);
            double value= rNode.GetValue(AUX_INDEX);
            rNode.SetValue(AUX_INDEX,value*0.5);
        }
             else{
            rNode.SetValue(AUX_INDEX,0.0);
        }
    });
}

                    // void HydraulicFluidAuxiliaryUtilities::FixCornerNodeVelocity(ModelPart &rModelPart, double angle_corner)
                    // {
                    //     // Obtain for each condition each neighbor condition
                    //     block_for_each(rModelPart.Nodes(), [&](NodeType &rNode){ rNode.SetValue(AUX_INDEX, 0); });
                    //     ModelPart::ConditionsContainerType &r_cond = rModelPart.Conditions();
                    //     double angle_corner_rad = angle_corner / 180.0 * 3.14159; //
                    //     double acceptable_cos = cos(angle_corner_rad);

                    //     for (ModelPart::ConditionsContainerType::iterator cond_it = r_cond.begin(); cond_it != r_cond.end(); cond_it++){
                    //         auto &r_geom = cond_it->GetGeometry();
                    //         // KRATOS_WATCH("ENTRA LA PRIMERA CONDICION")
                    //         r_geom = cond_it->GetGeometry();
                    //         // KRATOS_WATCH(r_geom)
                    //         auto edgelist = r_geom.GenerateEdges();
                    //         // reference for area normal of the face
                    //         const array_1d<double, 3> &face_normal = cond_it->GetValue(NORMAL);
                    //         double An = norm_2(face_normal);
                    //         const GlobalPointersVector<Condition> &neighb = cond_it->GetValue(NEIGHBOUR_CONDITIONS);
                    //         for (unsigned int c_itr = 0; c_itr < neighb.size(); c_itr++){

                    // auto rCondition = neighb[c_itr];
                    // const array_1d<double, 3> &face_neig_normal = rCondition.GetValue(NORMAL);
                    // double cos_angle = 1 / (An * norm_2(face_neig_normal)) * inner_prod(face_normal, face_neig_normal);
                    // auto &r_geom_neig = rCondition.GetGeometry();
                    // auto edgelist_neig = r_geom_neig.GenerateEdges();
                    // KRATOS_WATCH("VAMOS POR CADA NEIGHBOUR")
                    // KRATOS_WATCH(r_geom_neig)
                    // array_1d<double, 2> ids_neig;
                    // array_1d<double, 2> ids_cond;
                    // if (cos_angle < acceptable_cos)
                    // {
                    //     KRATOS_WATCH("EL ANGULO ES MAYOR DEL QUE HEMOS DICHO")
                    //     for (IndexType edge = 0; edge < edgelist.size(); edge++)
                    //     {
                    //         for (IndexType i = 0; i < edgelist[edge].size(); i++)
                    //         {
                    //             ids_cond[i] = edgelist[edge][i].Id();
                    //         }
                    //         std::sort(ids_cond.begin(), ids_cond.end());
                    //         KRATOS_WATCH("CADA EDGE CONDICION")
                    //         KRATOS_WATCH(ids_cond)
                    //         for (IndexType edge_neg = 0; edge_neg < edgelist_neig.size(); edge_neg++)
                    //         {
                    //             for (IndexType i = 0; i < edgelist_neig[edge_neg].size(); i++)
                    //             {
                    //                 ids_neig[i] = edgelist_neig[edge_neg][i].Id();
                    //             }
                    //             std::sort(ids_neig.begin(), ids_neig.end());
                    //             KRATOS_WATCH("CADA EDGE NEIG")
                    //             KRATOS_WATCH(ids_neig)
                    //             if (ids_cond == ids_neig)
                    //             {
                    //                 KRATOS_WATCH("CUMPLE")
                    //                 KRATOS_WATCH(ids_cond)
                    //                 KRATOS_WATCH(ids_neig)
                    //                 for (IndexType i = 0; i < edgelist[edge].size(); i++)
                    //                 {
                    //                     double repeated = edgelist[edge][i].GetValue(AUX_INDEX);
                    //                     repeated += 1;
                    //                     edgelist[edge][i].SetValue(AUX_INDEX, repeated);
                                        //                             }
                                        //                         }
                                        //                     }
                                        //                 }
                                        //             }
                                        //         }
                                        //     }

                                        //     block_for_each(rModelPart.Nodes(), [&](NodeType &rNode){
                                        //         if (rNode.GetValue(AUX_INDEX)>2){
                                        //             std::cout << rNode.GetValue(AUX_INDEX) << " " << rNode.Id() << std::endl;
                                        //             Vector veloctiy = ZeroVector(3);
                                        //             rNode.FastGetSolutionStepValue(VELOCITY) = veloctiy;
                                        //             rNode.Fix(VELOCITY_X);
                                        //             rNode.Fix(VELOCITY_Y);
                                        //             rNode.Fix(VELOCITY_Z);
                                        //         }
                                        //     });
                                        // }
bool HydraulicFluidAuxiliaryUtilities::MaximumWaterDepthChange(ModelPart &rModelPart)
{
    // make sure the Ids of the Nodes start with one
    auto min_Z = block_for_each<MinReduction<double>>(rModelPart.Nodes(), [&](const auto &rNode)
                                                      { return rNode.Z(); });

    ;
    double maximum_water_depth_change = false;
    block_for_each(rModelPart.Nodes(), [&](Node &rNode)
                   {
        if (rNode.Z() == min_Z) {

            double inlet_phi = rNode.GetValue(AUX_DISTANCE);
            double domain_phi = rNode.FastGetSolutionStepValue(DISTANCE);
            if (std::abs(domain_phi) > std::abs(inlet_phi))
            {
                maximum_water_depth_change=true;
            }
            else
            {
                maximum_water_depth_change = false;
            }


        } });
    return maximum_water_depth_change;
}
                                    

void HydraulicFluidAuxiliaryUtilities::CalculateArtificialViscosity(ModelPart &rModelPart,double WaterDynamicViscosityMax){
    const auto &r_process_info = rModelPart.GetProcessInfo();
    block_for_each(rModelPart.Elements(), [&](Element &rElement)
    {
    double artificial_viscosity;
    rElement.Calculate(ARTIFICIAL_DYNAMIC_VISCOSITY,artificial_viscosity ,r_process_info);
        if (artificial_viscosity > WaterDynamicViscosityMax)
        {
            artificial_viscosity = WaterDynamicViscosityMax;
        }
        double neg_nodes = 0.0;
        double pos_nodes=0.0;
        for (auto &r_node : rElement.GetGeometry())
        {
            double distance = r_node.FastGetSolutionStepValue(DISTANCE);
            
            if (distance > 0){
                pos_nodes += 1;
            }
            else{
                neg_nodes += 1;
            }
        }
        if (neg_nodes > 0 && pos_nodes > 0){
            artificial_viscosity = 0.0;
            }
            rElement.SetValue(ARTIFICIAL_DYNAMIC_VISCOSITY, artificial_viscosity);
    });
}

} // namespace Kratos