//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main author:    Uxue Chasco
//                  Ruben Zorrilla
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "includes/gid_io.h"
#include "geometries/geometry_data.h"
#include "includes/global_pointer_variables.h"
#include "processes/find_nodal_h_process.h"
#include "spatial_containers/bins_dynamic.h"
#include "utilities/binbased_fast_point_locator.h"
#include "utilities/mls_shape_functions_utility.h"
#include "utilities/variable_utils.h"

// Application includes
#include "utilities/rbf_shape_functions_utility.h"
#include "two_fluid_history_projection_utility.h"

namespace Kratos
{
    /* Public functions *******************************************************/

    void TwoFluidHistoryProjectionUtility::CalculateHistoryProjection(
        ModelPart& rModelPart,
        const bool ComputeNodalH,
        double ParticleLayerThickness,
        double SearchFactor)
    {
        // 1. Nodes that has different density between two time consecutive time steps are selected and element close to that nodes according to a used defined area (ParticleLayerThickness) are selected.
        FlagElementsAndNodes(rModelPart, ParticleLayerThickness);

        // 2. Calculates the final coordinates of the lagrangian particles and each corresponding velocity. Both values are storaged in a data container.
        auto particle_data = SeedAndConvectParticles(rModelPart,ParticleLayerThickness);

        // 3. Nodal edge is calculated.
        if (ComputeNodalH) {
            FindNodalHProcess<FindNodalHSettings::SaveAsNonHistoricalVariable> nodal_h_calculator(rModelPart);
            nodal_h_calculator.Execute();
        }

        // 4. New previous velocity is predicted after interpolating velocity of the lagrangian particles.
        std::size_t domain_size = rModelPart.GetProcessInfo()[DOMAIN_SIZE];
        if (domain_size == 2) {
            CalculateLagrangianVelocityInterpolation<2>(rModelPart, particle_data, ParticleLayerThickness,SearchFactor);
        } else if (domain_size == 3) {
            CalculateLagrangianVelocityInterpolation<3>(rModelPart, particle_data, ParticleLayerThickness,SearchFactor);
        } else {
            KRATOS_ERROR << "Wrong 'DOMAIN_SIZE' " << domain_size << "." << std::endl;
        }
    }

    /* Private functions *******************************************************/

    void TwoFluidHistoryProjectionUtility::FlagElementsAndNodes(
        ModelPart& rModelPart,
        double ParticleLayerThickness)
    {
        // Reset flags
        // SELECTED: Nodes and elements used for creating the seeds
        VariableUtils().SetFlag(SELECTED, false, rModelPart.Nodes());
        VariableUtils().SetFlag(SELECTED, false, rModelPart.Elements());
        VariableUtils().SetFlag(FREE_SURFACE, false, rModelPart.Nodes());
        VariableUtils().SetFlag(FREE_SURFACE, false, rModelPart.Elements());

        // Select the nodes of the area where the lagrangian particles are going to seed according
        //to user defined equidistance to previous time step free surface.
        for (auto& r_node : rModelPart.Nodes()) {
            const double distance_n = r_node.FastGetSolutionStepValue(DISTANCE,1);
            if (std::abs(distance_n) - ParticleLayerThickness < 0.0) {
                r_node.Set(SELECTED);
            }
        }

        // Select the elements of the area where the lagrangian particles are going to seed according
        // to user defined equidistance to previous time step free surface.
        block_for_each(rModelPart.Elements(), [](Element& rElement){
            const auto& r_geom = rElement.GetGeometry();
            for (const auto& rNode : r_geom) {
                if (rNode.Is(SELECTED)) {
                    rElement.Set(SELECTED, true);
                }
                break;
            }
        });
    }

    TwoFluidHistoryProjectionUtility::ParticleDataContainerType TwoFluidHistoryProjectionUtility::SeedAndConvectParticles(
        ModelPart& rModelPart,
        double ParticleLayerThickness)
    {
        // Initialize required data
        const double dt = rModelPart.GetProcessInfo()[DELTA_TIME];
        KRATOS_ERROR_IF(dt < 1.0e-12) << "Found DELTA_TIME close to zero in ProcessInfo." << std::endl;

        // Allocate auxiliary arrays
        ParticleDataContainerType particle_data;

        // Get and allocate arrays according to 1st element
        // Note that we are assuming a unique element type in the mesh
        const auto integration_method = GeometryData::IntegrationMethod::GI_GAUSS_3;
        const auto& elem_begin_geom = rModelPart.ElementsBegin()->GetGeometry();
        const std::size_t n_nodes = elem_begin_geom.PointsNumber();
        const std::size_t n_gauss = elem_begin_geom.IntegrationPointsNumber(integration_method);

        // Set weights
        // For SELECTED nodes, a linear function is used to calculate the weight
        // For not SELECTED nodes a zero weight (Eulerian) is set
        for (auto& r_node : rModelPart.Nodes()) {
            double distance_ratio = 0.0;
            double old_distance=r_node.FastGetSolutionStepValue(DISTANCE,1);
            if (r_node.Is(SELECTED)) {
                distance_ratio = (1 - std::abs(old_distance) / ParticleLayerThickness);
            } else {
                distance_ratio = 0.0;
            }
            r_node.SetValue(TEMPERATURE,1.0);
        }

        // Loop the SELECTED elements
        // TODO: Make this parallel
        array_1d<double,3> a;
        array_1d<double,3> v_n;
        array_1d<double, 3> v_nn;
        array_1d<double,3> aux_coords;
        array_1d<double,3> initial_coords;
        array_1d<double, 3> acceleration;
        array_1d<double, 3> velocity_field_n;
        array_1d<double, 3> aux_coords_node;
        array_1d<double, 3> initial_coords_node;
        array_1d<double, 3> velocity_node_particle;

        std::size_t id = rModelPart.NumberOfNodes();
        auto& r_model = rModelPart.GetModel();
        auto& r_post_process_particles_mp = r_model.CreateModelPart("PostProcessParticles");

        for (const auto& r_elem : rModelPart.Elements()) {
            if (r_elem.Is(SELECTED)) {
                const auto& r_geom = r_elem.GetGeometry();
                const auto& r_N_container = r_geom.ShapeFunctionsValues(integration_method);

                for (std::size_t i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
                    const auto& r_N = row(r_N_container, i_gauss);
                    ProcessInfo &r_process_info = rModelPart.GetProcessInfo();
                    const array_1d<double, 3> &gravity = r_process_info[GRAVITY];
                    // Interpolate the current particle coordinates and velocity
                    // Note that we use the velocity in the first buffer position as this is the une used in the level set convection
                    double distance=0.0;
                    double v_weight = 0.0;
                    double distance_old=0.0;
                    noalias(a) = ZeroVector(3);
                    v_n = ZeroVector(3);
                    v_nn = ZeroVector(3);
                    aux_coords = ZeroVector(3);
                    initial_coords = ZeroVector(3);
                    for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
                        v_weight += r_N[i_node] * r_geom[i_node].GetValue(TEMPERATURE);
                        distance += r_N[i_node] * r_geom[i_node].FastGetSolutionStepValue(DISTANCE, 1);
                        distance_old += r_N[i_node] * r_geom[i_node].FastGetSolutionStepValue(DISTANCE, 2);
                        noalias(a) += r_N[i_node] * r_geom[i_node].GetValue(ACCELERATION);
                        noalias(v_n) += r_N[i_node] * r_geom[i_node].FastGetSolutionStepValue(VELOCITY,1);
                        noalias(v_nn) += r_N[i_node] * r_geom[i_node].FastGetSolutionStepValue(VELOCITY, 2);
                        noalias(aux_coords) += r_N[i_node] * r_geom[i_node].Coordinates();
                        noalias(initial_coords) += r_N[i_node] * r_geom[i_node].Coordinates();
                    }

                    // Calculate the current particle final position
                    if (distance<0.0) {
                        noalias(aux_coords) +=dt * v_weight * v_n;
                    }
                    else{
                        noalias(aux_coords) +=dt * v_weight * v_n+gravity*0.5*dt*dt;

                    }
                    double levelset_increment=distance-distance_old;
                    // Save data in the particle data container
                    auto p_new_particle = Kratos::make_shared<ParticleDataType>(aux_coords, v_n, v_weight * v_n, distance, initial_coords, a, ++id, levelset_increment);
                    particle_data.push_back(p_new_particle);

                    // Create a node from each particle to postprocess
                    auto p_node = r_post_process_particles_mp.CreateNewNode(id, initial_coords[0], initial_coords[1], initial_coords[2]);
                    p_node->Coordinates() = aux_coords;
                    p_node->SetValue(TEMPERATURE, v_weight);
                    p_node->SetValue(DISTANCE, distance);
                    p_node->SetValue(ACCELERATION, a);
                    p_node->SetValue(VELOCITY, v_n);
                    p_node->SetValue(DISPLACEMENT, dt * v_weight * v_n);
                }




            }
        }
        // velocity_field_n = ZeroVector(3);
        // velocity_node_particle = ZeroVector(3);
        // double distance = 0.0;
        // double weight_node = 0.0;
        // acceleration = ZeroVector(3);
        // initial_coords_node = ZeroVector(3);
        // aux_coords_node = ZeroVector(3);

        // for (auto &r_node : rModelPart.Nodes())
        // {
        //     if( r_node.Is(SELECTED))
        //     {

        //     velocity_field_n = r_node.FastGetSolutionStepValue(VELOCITY, 1);
        //     weight_node = r_node.FastGetSolutionStepValue(TEMPERATURE);
        //     velocity_node_particle = velocity_field_n * weight_node;
        //     distance = r_node.FastGetSolutionStepValue(DISTANCE, 1);
        //     initial_coords_node = r_node.Coordinates();
        //     acceleration = r_node.GetValue(ACCELERATION);
        //     aux_coords_node = initial_coords_node + dt * velocity_node_particle;
        //     auto p_new_particle = Kratos::make_shared<ParticleDataType>(aux_coords_node, velocity_field_n, velocity_node_particle, distance, initial_coords_node, acceleration, ++id);

        //     // Create a node from each particle to postprocess
        //     auto p_node = r_post_process_particles_mp.CreateNewNode(++id, initial_coords_node[0], initial_coords_node[1], initial_coords_node[2]);
        //     p_node->Coordinates() = aux_coords_node;
        //     p_node->SetValue(TEMPERATURE, weight_node);
        //     p_node->SetValue(DISTANCE, distance);
        //     p_node->SetValue(ACCELERATION, acceleration);
        //     p_node->SetValue(VELOCITY, velocity_field_n);
        //     p_node->SetValue(DISPLACEMENT, dt * velocity_field_n * weight_node);
        //     }
        // }
        for (auto& r_node : rModelPart.Nodes()) {
            r_node.SetValue(YOUNG_MODULUS, r_node.FastGetSolutionStepValue(DISTANCE, 1));
        }

        // GidIO<> gid_io_fluid("/home/uchasco/Desktop/ProjectionUtilityCppTests/auxiliary_print_" + std::to_string(rModelPart.GetProcessInfo()[STEP]), GiD_PostAscii, SingleFile, WriteUndeformed, WriteConditions);
        GidIO<> gid_io_fluid("/home/uchasco/Desktop/CasosJunioParticle/sloshing_2d/particles_" + std::to_string(rModelPart.GetProcessInfo()[STEP]), GiD_PostAscii, SingleFile, WriteUndeformed, WriteConditions);
        gid_io_fluid.InitializeMesh(0.0);
        gid_io_fluid.WriteMesh(rModelPart.GetMesh());
        gid_io_fluid.WriteNodeMesh(r_post_process_particles_mp.GetMesh());
        gid_io_fluid.FinalizeMesh();
        gid_io_fluid.InitializeResults(0, r_post_process_particles_mp.GetMesh());
        gid_io_fluid.WriteNodalResultsNonHistorical(TEMPERATURE, rModelPart.Nodes(), 0);
        gid_io_fluid.WriteNodalResults(DISTANCE, rModelPart.Nodes(), 0, 0);
        gid_io_fluid.WriteNodalResultsNonHistorical(YOUNG_MODULUS, rModelPart.Nodes(), 0);
        gid_io_fluid.WriteNodalResultsNonHistorical(ACCELERATION, rModelPart.Nodes(), 0);
        gid_io_fluid.WriteNodalResults(VELOCITY, rModelPart.Nodes(), 0, 0);
        gid_io_fluid.WriteNodalFlags(SELECTED, "SELECTED", rModelPart.Nodes(), 0);
        gid_io_fluid.WriteNodalResultsNonHistorical(TEMPERATURE, r_post_process_particles_mp.Nodes(), 0);
        gid_io_fluid.WriteNodalResultsNonHistorical(DISTANCE, r_post_process_particles_mp.Nodes(), 0);
        gid_io_fluid.WriteNodalResultsNonHistorical(ACCELERATION, r_post_process_particles_mp.Nodes(), 0);
        gid_io_fluid.WriteNodalResultsNonHistorical(VELOCITY, r_post_process_particles_mp.Nodes(), 0);
        gid_io_fluid.WriteNodalResultsNonHistorical(DISPLACEMENT, r_post_process_particles_mp.Nodes(), 0);
        gid_io_fluid.FinalizeResults();

        r_model.DeleteModelPart("PostProcessParticles");

        return particle_data;
    }

    template<std::size_t TDim>
    void TwoFluidHistoryProjectionUtility::CalculateLagrangianVelocityInterpolation(
        ModelPart& rModelPart,
        ParticleDataContainerType& rParticleData,
        double ParticleLayerThickness,
        double SearchFactor)
    {
        constexpr std::size_t min_n_particles = TDim + 1;

        // Reset MESH_VELOCITY
        VariableUtils().SetHistoricalVariableToZero(MESH_VELOCITY, rModelPart.Nodes());

        // Set the bins search
        typedef BinsDynamic<TDim, ParticleDataType, ParticleDataContainerType> BinsDynamicSearchType;
        BinsDynamicSearchType bins_search(rParticleData.begin(), rParticleData.end());

        // Create the bin-based point locator to find the FREE_SURFACE nodes
        // These are the ones in which the moved particles lie, that is, where the MLS is applied
        Vector N;
        Element::Pointer p_elem = nullptr;
        // ParticleDataContainerType aux_particle_data;
        BinBasedFastPointLocator<TDim> bin_based_particle_locator(rModelPart);
        bin_based_particle_locator.UpdateSearchDatabase();
        for (auto p_part : rParticleData) {
            const bool found = bin_based_particle_locator.FindPointOnMeshSimplified(p_part->Coordinates, N, p_elem);
            if (found) {
                p_elem->Set(FREE_SURFACE, true);
                auto& r_geom = p_elem->GetGeometry();
                for (auto& r_node : r_geom) {
                    r_node.Set(FREE_SURFACE, true);
                }

                // double dist_int = 0.0;
                // array_1d<double,3> v_int = ZeroVector(3);
                // array_1d<double,3> acc_int = ZeroVector(3);
                // for (std::size_t i = 0; i < r_geom.PointsNumber(); ++i) {
                //     dist_int += N(i) * r_geom[i].FastGetSolutionStepValue(DISTANCE,1);
                //     v_int += N(i) * r_geom[i].FastGetSolutionStepValue(VELOCITY,1);
                //     acc_int += N(i) * r_geom[i].GetValue(ACCELERATION);
                // }

                // p_part->Distance -= dist_int;
                // p_part->OldVelocity -= v_int;
                // p_part->OldAcceleration -= acc_int;

                // aux_particle_data.push_back(p_part);
            }
        }
        // rParticleData.swap(aux_particle_data);
        // aux_particle_data.clear();

        // Compute the Lagrangian velocity interpolation
        const std::size_t max_results = 50;
        for (auto& r_node : rModelPart.Nodes()) {
            if (r_node.Is(FREE_SURFACE)) {
                const auto& r_node_coords = r_node.Coordinates();

                // Search the closest particles to the current node
                ParticleDataType aux_node(r_node_coords);
                ParticleDataContainerType cloud_particles(max_results);
                const double search_dist = SearchFactor * r_node.GetValue(NODAL_H) / 2.0;
                const std::size_t n_cloud_points = bins_search.SearchInRadius(
                    aux_node,
                    search_dist,
                    cloud_particles.begin(),
                    max_results);

                if(n_cloud_points >= min_n_particles)
                {

                    std::vector<std::size_t> part_ids; //TODO: remove after debug

                    // Set the coordinates matrix for the MLS interpolation
                    Matrix cloud_coords = ZeroMatrix(n_cloud_points,3);
                    for (std::size_t i = 0; i < n_cloud_points; ++i) {
                        const auto& r_particle = *(cloud_particles[i]);
                        for (std::size_t d = 0; d < TDim; ++d) {
                            cloud_coords(i,d) = r_particle[d];
                        }
                        part_ids.push_back(r_particle.ParticleId); //TODO: remove after debug
                    }

                    bool repeated;
                    const double tolerance=1e-7;
                    std::vector<std::size_t> to_keep;
                    for (std::size_t i=0; i<n_cloud_points; i++){
                        repeated = false;
                        const auto& r_i_coord = row(cloud_coords, i);
                        for (std::size_t j=i+1; j < n_cloud_points; j++) {
                            const auto& r_j_coord = row(cloud_coords, j);
                            if (norm_2(r_i_coord - r_j_coord) < tolerance) {
                                repeated = true;
                                break;
                            }
                        }
                        if (!repeated) {
                            to_keep.push_back(i);
                        }
                    }

                    //TODO: Think a better solution than the error
                    KRATOS_ERROR_IF(to_keep.size() < min_n_particles) << "Not enough particles to interpolate the velocity" << std::endl;


                    if (!(n_cloud_points == to_keep.size())) {
                        KRATOS_WARNING("TwoFluidHistoryProjectionUtility") << "There are repeated particles." << std::endl;
                        cloud_coords = ZeroMatrix(to_keep.size(), 3);
                        for (std::size_t i_part = 0; i_part< to_keep.size(); i_part++) {
                            const std::size_t part_id = to_keep[i_part];
                            const auto& r_particle = *(cloud_particles[part_id]);
                            for (std::size_t d = 0; d < TDim; ++d) {
                                cloud_coords(i_part,d) = r_particle[d];
                            }
                        }
                    }

                        // Calculate the MLS interpolation in the FREE_SURFACE node
                        // TODO: Use an auxiliary lambda to avoid this if
                    Vector cloud_shape_functions;
                    std::size_t domain_size = rModelPart.GetProcessInfo()[DOMAIN_SIZE];
                    if (domain_size == 2)
                    {
                        MLSShapeFunctionsUtility::CalculateShapeFunctions<2>(cloud_coords, r_node_coords, search_dist, cloud_shape_functions);
                    }
                    else if (domain_size == 3)
                    {
                        MLSShapeFunctionsUtility::CalculateShapeFunctions<3>(cloud_coords, r_node_coords, search_dist, cloud_shape_functions);
                    }
                    else
                    {
                        KRATOS_ERROR << "Wrong DOMAIN_SIZE." << std::endl;
                    }

                    // RBFShapeFunctionsUtility::CalculateShapeFunctions(cloud_coords, r_node_coords, cloud_shape_functions);

                    // Interpolate the FREE_SURFACE node history and mesh velocity
                    array_1d<double,3> a_n = ZeroVector(3);
                    array_1d<double,3> v_n = ZeroVector(3);
                    array_1d<double,3> v_mesh = ZeroVector(3);
                    array_1d<double, 3> x_f = ZeroVector(3);
                    array_1d<double, 3> x_i = ZeroVector(3);
                    array_1d<double, 3> acceleration = ZeroVector(3);
                    double dt=rModelPart.GetProcessInfo()[DELTA_TIME];
                    double dist_n=0.0;
                    double dist_n_1=0.0;
                    for (std::size_t i = 0; i < cloud_coords.size1(); ++i) {
                        const auto& r_particle = *(cloud_particles[i]);
                        noalias(a_n) += cloud_shape_functions(i) * r_particle.OldAcceleration;
                        noalias(v_n) += cloud_shape_functions(i) * r_particle.OldVelocity;
                        noalias(x_f)+= cloud_shape_functions(i) * r_particle.Coordinates;
                        noalias(x_i) += cloud_shape_functions(i) * r_particle.InitialCoordinates;

                        dist_n+= cloud_shape_functions(i) * r_particle.Distance;
                        dist_n_1 += cloud_shape_functions(i)*r_particle.DistanceOld;
                    }
                    v_mesh = (x_f - x_i) / dt;

                    acceleration=a_n*dt;
                    // KRATOS_WATCH(v_n)
                    // In order to avoid changing abruptly mesh velocity between adjacent nodes, a progessive change has been done weighting old velocity and new predicted velocity along all the area defined by the user.
                    // The weights are linearly dependant with the ParticleLayerThickness.Close to the area where the nodes have changes its density, new predictd velocity will have more weight.



                    // Weighted velocity is calculated.

                    // Set new interpolated velocity such that the interface is Lagrangian in an approximate sense.
                    // It should be considered the fixity. For that nodes that ones of each components is fixed it is not applied false fm-ale to that component.Same procedure has been  for SLIP CONDITION,

                    // Do the levelset particle convection
                    // Note that this is done regardless of the wall boundary condition type
                    r_node.FastGetSolutionStepValue(DISTANCE) += dist_n_1;

                    if (r_node.IsNot(SLIP))
                    {
                        if (!r_node.IsFixed(VELOCITY_X))
                        {
                            r_node.GetValue(ACCELERATION_X) -= a_n[0];
                            r_node.FastGetSolutionStepValue(VELOCITY_X) -= v_n[0];
                            r_node.FastGetSolutionStepValue(VELOCITY_X, 1) -= v_n[0];
                            r_node.FastGetSolutionStepValue(MESH_VELOCITY_X) = v_mesh[0];
                        }

                        if (!r_node.IsFixed(VELOCITY_Y))
                        {
                            r_node.GetValue(ACCELERATION_Y) -= a_n[1];
                            r_node.FastGetSolutionStepValue(VELOCITY_Y) += v_n[1];
                            r_node.FastGetSolutionStepValue(VELOCITY_Y, 1) += v_n[1];
                            r_node.FastGetSolutionStepValue(MESH_VELOCITY_Y) = v_mesh[1];
                        }

                        if (!r_node.IsFixed(VELOCITY_Z))
                        {
                            r_node.GetValue(ACCELERATION_Z) -= a_n[2];
                            r_node.FastGetSolutionStepValue(VELOCITY_Z) += v_n[2];
                            r_node.FastGetSolutionStepValue(VELOCITY_Z,1) += v_n[2];
                            r_node.FastGetSolutionStepValue(MESH_VELOCITY_Z) = v_mesh[2];
                        }
                    }
                    else
                    {
                        array_1d<double, 3> unit_normal = r_node.FastGetSolutionStepValue(NORMAL);
                        const double normal_norm = norm_2(unit_normal);
                        if (normal_norm > 1.0e-12)
                        {
                            unit_normal /= normal_norm;
                        }
                        else
                        {
                            KRATOS_WARNING("TwoFluidsHistoryProjectionUtility") << "Normal in node " << r_node.Id() << " is close to zero." << std::endl;
                        }
                        const array_1d<double, 3> a_n_norm = inner_prod(a_n, unit_normal) * unit_normal;
                        const array_1d<double, 3> a_n_tang = a_n - a_n_norm;
                        const array_1d<double, 3> v_n_norm = inner_prod(v_n, unit_normal) * unit_normal;
                        const array_1d<double, 3> v_n_tang = v_n - v_n_norm;
                        const array_1d<double, 3> v_mesh_norm = inner_prod(v_mesh, unit_normal) * unit_normal;
                        const array_1d<double, 3> v_mesh_tang = v_mesh - v_mesh_norm;

                        r_node.GetValue(ACCELERATION) -= a_n_tang;
                        r_node.FastGetSolutionStepValue(VELOCITY) -= v_n_tang;
                        r_node.FastGetSolutionStepValue(VELOCITY, 1) -= v_n_tang;
                        r_node.FastGetSolutionStepValue(MESH_VELOCITY) = v_mesh_tang;
                    }
                }
                else
                {
                    if (r_node.Is(FREE_SURFACE))
                    {
                        KRATOS_WARNING("TwoFluidHistoryProjectionUtility") << "Not enough particles for node " << r_node.Id() << "." << std::endl;
                    }
                }
            }
        }

        // Resetting the distance values in the region far from the level set
        // Note that these are spurious and might affect the redistance calculation
        const double aux_val = 1.0e10;
        block_for_each(rModelPart.Nodes(), [&aux_val](ModelPart::NodeType& r_node) {
            if (!r_node.Is(FREE_SURFACE)) {
                double& r_dist = r_node.FastGetSolutionStepValue(DISTANCE);
                r_dist = r_dist > 0.0 ? aux_val : -aux_val;
            }
        });
    }

    /* External functions *****************************************************/

    /// output stream function for TwoFluidHistoryProjectionUtility
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const TwoFluidHistoryProjectionUtility& rThis)
    {
        rThis.PrintData(rOStream);
        return rOStream;
    }

    /// output stream function for ParticleData class
    inline std::ostream& operator<<(
        std::ostream& rOStream,
        const TwoFluidHistoryProjectionUtility::ParticleData& rThis)
    {
        rOStream << "ParticleData" << std::endl;
        rOStream << "\tCoordinates: [" << rThis.Coordinates[0] << "," << rThis.Coordinates[1] << "," << rThis.Coordinates[2] << "]" << std::endl;
        rOStream << "\tOldVelocity: [" << rThis.OldVelocity[0] << "," << rThis.OldVelocity[1] << "," << rThis.OldVelocity[2] << "]" << std::endl;
        rOStream << "\tParticleVelocity: [" << rThis.ParticleVelocity[0] << "," << rThis.ParticleVelocity[1] << "," << rThis.ParticleVelocity[2] << "]" << std::endl;
        return rOStream;
    }

    /* Explicit template instantiation *****************************************************/
    template void TwoFluidHistoryProjectionUtility::CalculateLagrangianVelocityInterpolation<2>(
        ModelPart& rModelPart,
        ParticleDataContainerType& rParticleData,
        double ParticleLayerThickness,
        double SearchFactor);

    template void TwoFluidHistoryProjectionUtility::CalculateLagrangianVelocityInterpolation<3>(
        ModelPart& rModelPart,
        ParticleDataContainerType& rParticleData,
        double ParticleLayerThickness,
        double SearchFactor);

}

