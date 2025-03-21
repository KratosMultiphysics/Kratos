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

void HydraulicFluidAuxiliaryUtilities::CalculateNonIntersectedElementsArtificialViscosity(
    ModelPart &rModelPart,
    double DynamicViscosityMax)
    {
    const auto &r_process_info = rModelPart.GetProcessInfo();
    block_for_each(rModelPart.Elements(), [&](Element &rElement)
    {
        double artificial_viscosity;
        // Check if the element is cut
        double neg_nodes = 0.0;
        double pos_nodes=0.0;
        for (auto &r_node : rElement.GetGeometry())
        {
            double distance = r_node.FastGetSolutionStepValue(DISTANCE);
            if (distance > 0)
            {
                pos_nodes += 1;
            }
            else
            {
                neg_nodes += 1;
            }
        }
        // Calculate the artificial viscosity
        if (neg_nodes > 0 && pos_nodes > 0)
        {
            // If the element is cut, the artificial viscosity is set to zero.
            artificial_viscosity = 0.0;
        }
        else
        {
            // If the element is NOT cut, the artificial viscosity is calculated.
            rElement.Calculate(ARTIFICIAL_DYNAMIC_VISCOSITY,artificial_viscosity ,r_process_info);
            if (artificial_viscosity > DynamicViscosityMax)
            {
                artificial_viscosity = DynamicViscosityMax;
            }
        }
        rElement.SetValue(ARTIFICIAL_DYNAMIC_VISCOSITY, artificial_viscosity);
    });
}

void HydraulicFluidAuxiliaryUtilities::ApplyOutletInflowLimiter(
    ModelPart &rModelPart,
    const Variable<array_1d<double, 3>>& rVariable,
    const Variable<array_1d<double, 3>>& rVariableNormal)
{

    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode)
    {
        array_1d<double, 3>& r_velocity = rNode.FastGetSolutionStepValue(rVariable);
        // We use a non-historical variable in case rVariableNormal is not the NORMAL variable and an auxiliary variable is used
        const array_1d<double, 3>& r_normal = rNode.GetValue(rVariableNormal);
        const double norm_n = norm_2(r_normal);
        if (norm_n > 0.0)
        {
            array_1d<double, 3> n_unit = r_normal / norm_n;
            const double aux = inner_prod(r_velocity, n_unit);
            if (aux < 0.0)
            {
                noalias(r_velocity) -= aux * n_unit;
            }
        }
    });
}

void HydraulicFluidAuxiliaryUtilities::FixingInflow(
    ModelPart &rModelPart,
    const Variable<array_1d<double, 3>> &rVariable,
    double DomainSize)
{
    const std::array<const Variable<double>*, 3> *components = nullptr;

    // Asignamos el puntero dependiendo de la variable
    if (rVariable == VELOCITY)
    {
        static const std::array<const Variable<double>*, 3> velocity_components = {&VELOCITY_X, &VELOCITY_Y, &VELOCITY_Z};
        components = &velocity_components;
    }
    else if (rVariable == FRACTIONAL_VELOCITY)
    {
        static const std::array<const Variable<double>*, 3> fractional_velocity_components = {&FRACTIONAL_VELOCITY_X, &FRACTIONAL_VELOCITY_Y, &FRACTIONAL_VELOCITY_Z};
        components = &fractional_velocity_components;
    }
    else
    {
        KRATOS_ERROR << "The variable " << rVariable.Name() << " is not supported. Only VELOCITY and FRACTIONAL_VELOCITY are allowed." << std::endl;
    }
    block_for_each(rModelPart.Nodes(), [&](NodeType &rNode)
        {   rNode.Set(OUTLET);
            // array_1d<double, 3>& r_velocity = rNode.FastGetSolutionStepValue(rVariable);
            // // TODO: SE TENDRIA QUE PILLAR LA NORMAL DE ARRIBA NO??
            // const array_1d<double, 3>& r_normal = rNode.FastGetSolutionStepValue(NORMAL);
            // const double norm_n = norm_2(r_normal);
            // if (norm_n > 0.0) {
            //     array_1d<double, 3> n_unit = r_normal / norm_n;
            //     const double aux = inner_prod(r_velocity, n_unit);

            //     if (aux < 1e-7) {
            //         for (std::size_t i_dim = 0; i_dim < (DomainSize - 1); ++i_dim) {
            //             rNode.Fix(*(*components)[i_dim]);
            //         }
            //     }
            // }
        });
    block_for_each(rModelPart.Conditions(), [&](Condition & rCondition){
        rCondition.Set(OUTLET);
    });
}

void HydraulicFluidAuxiliaryUtilities::ImposeOutletPressure(ModelPart &rModelPart, double WaterDepth, const Variable<double> &rDistanceVariable)
{
    {
        // Obtain minimum Z coordinate. Assume that in 3D domain the height is given by the z coordinate.
        double gravity = 9.81;
        double min_z = std::numeric_limits<double>::max();
        min_z = block_for_each<MinReduction<double>>(rModelPart.Nodes(), [&](Node &rNode)
                                                     { return rNode.Z(); });
        // Obtain z of the free surface
        double z_free_surface = WaterDepth + min_z;

        block_for_each(rModelPart.Nodes(), [&](Node &rNode)
        {
            double aux_distance = rNode.Z() - z_free_surface;
            rNode.FastGetSolutionStepValue(rDistanceVariable) = aux_distance;
            double& distance = rNode.FastGetSolutionStepValue(DISTANCE);

            if (aux_distance < 0.0 && distance < 0.0)
            {
                const double& density = rNode.FastGetSolutionStepValue(DENSITY);
                double pressure_external = gravity * density * -aux_distance;
                rNode.FastGetSolutionStepValue(EXTERNAL_PRESSURE) = pressure_external;
            } });
    }
}

void HydraulicFluidAuxiliaryUtilities::SetBoundaryWaterDepth(ModelPart &rModelPart, double WaterDepth, const Variable<double> &rDistanceVariable){
    block_for_each(rModelPart.Nodes(), [&](Node &rNode)
    {
    double phi = rNode.Z() - WaterDepth;
    rNode.SetValue(rDistanceVariable, phi);
    }); }


void HydraulicFluidAuxiliaryUtilities::AssignInletWaterDepth(ModelPart &rModelPart,double InletVelocity,double DeltaTime){

    //  Projected Velocity
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
        array_1d<double,3> n = inlet_velocity/norm_2(inlet_velocity);
        array_1d<double,3> grad_phi = rNode.FastGetSolutionStepValue(DISTANCE_GRADIENT);
        double aux = inner_prod(n,grad_phi);
        if (aux >0.0){
            double phi = rNode.GetValue(AUX_DISTANCE);
            double grad_phi_norm = norm_2(grad_phi);
            double phi_new = phi + grad_phi_norm*InletVelocity*DeltaTime;

            rNode.SetValue(AUX_DISTANCE,phi_new);
        };


    });
}

void HydraulicFluidAuxiliaryUtilities::SummergedInletCheck(ModelPart &rModelPart, double InletVelocity, double DeltaTime)

{
    // Promedio de aux en nodos y condiciones cortadas
    double total_aux_nodes = 0.0;
    int total_nodes_count = 0;
    double total_aux_conditions = 0.0;
    int cut_conditions_count = 0;

    for (auto &rCondition : rModelPart.Conditions())
    {
        const int number_of_nodes = rCondition.GetGeometry().PointsNumber();
        double neg_nodes = 0.0;
        double pos_nodes = 0.0;
        double condition_aux_sum = 0.0; // Para acumular aux en esta condición

        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            auto &r_node = rCondition.GetGeometry()[i_node];
            if (r_node.GetValue(AUX_DISTANCE) > 0.0)
            {
                pos_nodes += 1.0;
            }
            else
            {
                neg_nodes += 1.0;
            }
        }

        if (neg_nodes != 0.0 && pos_nodes != 0.0)

        {
             // La condición está cortada
            for (int i_node = 0; i_node < number_of_nodes; ++i_node)
            {
                auto &r_node = rCondition.GetGeometry()[i_node];
                array_1d<double, 3> inlet_norm = r_node.GetValue(INLET_NORMAL);
                const double n_norm = norm_2(inlet_norm);

                if (n_norm > 1.0e-12)
                {
                    inlet_norm /= -n_norm;
                }
                else
                {
                    KRATOS_WARNING("SetInletVelocity")
                        << "Node " << r_node.Id() << " INLET_NORMAL is close to zero." << std::endl;
                    inlet_norm = ZeroVector(3);
                }

                // Inlet velocity vector.
                // if (r_node.GetValue(AUX_DISTANCE)>0.0){
                //     InletVelocity =0.0;
                // }
                array_1d<double, 3> inlet_velocity = inlet_norm * InletVelocity;
                const double vel_norm = norm_2(inlet_velocity);

                // Corrección del operador ternario
                array_1d<double, 3> n;

                if (vel_norm > 1.0e-12)
                {
                    n = inlet_velocity / vel_norm;
                }
                else
                {
                    n = ZeroVector(3);
                }

                array_1d<double, 3> grad_phi = r_node.FastGetSolutionStepValue(DISTANCE_GRADIENT);
                double aux = inner_prod(grad_phi, n);
                condition_aux_sum += aux; // Acumular aux para esta condición
                total_nodes_count += 1;   // Contar nodos procesados
            }

            total_aux_conditions += condition_aux_sum / number_of_nodes; // Media de aux en esta condición
            cut_conditions_count += 1;
        }
    }

    // Calcular la media final
    double average_aux_conditions = (cut_conditions_count > 0) ? total_aux_conditions / cut_conditions_count : 0.0;



    // Corrección en block_for_each
    block_for_each(rModelPart.Nodes(), [&](ModelPart::NodeType &rNode)
                   {
        double phi_old = rNode.GetValue(AUX_DISTANCE);
        double new_phi = phi_old +average_aux_conditions * DeltaTime;
        rNode.SetValue(AUX_DISTANCE, new_phi); });
}

void HydraulicFluidAuxiliaryUtilities::InletFreeSurface(ModelPart &rModelPart){
    block_for_each(rModelPart.Nodes(), [&](NodeType &rNode) {
        double phi_c = rNode.GetValue(AUX_DISTANCE);
        double phi_domain = rNode.FastGetSolutionStepValue(DISTANCE);
        if (phi_c >0.0 && phi_domain>0.0){
            if (phi_c-phi_domain>0.0){
                    rNode.SetValue(AUX_DISTANCE,phi_domain);
            }

        }
        else if (phi_c>0.0 && phi_domain<0.0){
            rNode.SetValue(AUX_DISTANCE, phi_domain);
        }
        else if (phi_c <0.0 && phi_domain<0.0){
            if (phi_c-phi_domain>0.0){
                    rNode.SetValue(AUX_DISTANCE,phi_domain);
            }

        }
    });
}
void HydraulicFluidAuxiliaryUtilities::ComputeNodalFroudeNumber(ModelPart &rModelPart)
{
    block_for_each(rModelPart.Nodes(), [&](NodeType &rNode)
                   {

        const double water_depth = rNode.GetValue(WATER_DEPTH);
        const array_1d<double, 3>& v_n_1 = rNode.FastGetSolutionStepValue(VELOCITY);
        const array_1d<double, 3>& g = rNode.FastGetSolutionStepValue(BODY_FORCE);

        array_1d<double, 3> local_froude_number;
        const double sqrt_water_depth = std::sqrt(water_depth);

        for (size_t i = 0; i < 3; ++i) {
            local_froude_number[i] = v_n_1[i] / (std::sqrt(g[i]) * sqrt_water_depth);
        }

        // rNode.SetValue(FROUDE_NUMBER, local_froude_number);
         });
}
} // namespace Kratos
