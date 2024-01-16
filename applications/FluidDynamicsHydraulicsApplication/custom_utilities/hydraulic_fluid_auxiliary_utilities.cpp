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
            rNode.Fix(DISTANCE);
        }
    });
}

} // namespace Kratos
