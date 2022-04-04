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
#include "geometries/geometry_data.h"
#include "includes/global_pointer_variables.h"
#include "processes/find_nodal_h_process.h"
#include "spatial_containers/bins_dynamic.h"
#include "utilities/mls_shape_functions_utility.h"
#include "utilities/variable_utils.h"

// Application includes
#include "radial_basis_functions_utility.h"
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
        auto particle_data = SeedAndConvectParticles(rModelPart);

        // 3. Nodal edge is calculated.
        if (ComputeNodalH) {
            FindNodalHProcess<FindNodalHSettings::SaveAsNonHistoricalVariable> nodal_h_calculator(rModelPart);
            nodal_h_calculator.Execute();
        }

        // 3. New previous velocity is predicted after interpolating velocity of the lagrangian particles.
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
        // FREE_SURFACE: Nodes and elements where the history is to be reset
        VariableUtils().SetFlag(MARKER, false, rModelPart.Nodes());
        VariableUtils().SetFlag(SELECTED, false, rModelPart.Nodes());
        VariableUtils().SetFlag(SELECTED, false, rModelPart.Elements());
        VariableUtils().SetFlag(FREE_SURFACE, false, rModelPart.Nodes());

        // Flag the nodes that have changed its phase
        block_for_each(rModelPart.Nodes(), [](Node<3>& rNode){
            const double new_d = rNode.FastGetSolutionStepValue(DISTANCE);
            const double old_d = rNode.FastGetSolutionStepValue(DISTANCE,1);
            if (new_d * old_d < 0.0) {
                rNode.Set(SELECTED, true);
                rNode.Set(FREE_SURFACE, true);
            }
        });
        // Select the nodes of the area where the lagrangian particles are going to seed according
        //to user defined equidistance to previous time step free surface.
        for (auto& r_node : rModelPart.Nodes()) {
            const double distance_n = r_node.FastGetSolutionStepValue(DISTANCE,1);
            if (distance_n>0.0){
                if(distance_n-ParticleLayerThickness<0.0){
                    r_node.Set(SELECTED,true);
                }
            }
            else{
                if(distance_n+ParticleLayerThickness>0.0){
                    r_node.Set(SELECTED,true);
                }
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
            }
        });
    }

    TwoFluidHistoryProjectionUtility::ParticleDataContainerType TwoFluidHistoryProjectionUtility::SeedAndConvectParticles(ModelPart& rModelPart)
    {
        // Initialize required data
        const double dt = rModelPart.GetProcessInfo()[DELTA_TIME];
        KRATOS_ERROR_IF(dt < 1.0e-12) << "Found DELTA_TIME close to zero in ProcessInfo." << std::endl;

        // Allocate auxiliary arrays
        ParticleDataContainerType particle_data;

        // Get and allocate arrays according to 1st element
        // Note that we are assuming a unique element type in the mesh
        const auto integration_method = GeometryData::IntegrationMethod::GI_GAUSS_1;
        const auto& elem_begin_geom = rModelPart.ElementsBegin()->GetGeometry();
        const std::size_t n_nodes = elem_begin_geom.PointsNumber();
        const std::size_t n_gauss = elem_begin_geom.IntegrationPointsNumber(integration_method);

        // Loop the SELECTED elements
        // TODO: Make this parallel
        array_1d<double,3> aux_v;
        array_1d<double,3> aux_coords;
        for (const auto& r_elem : rModelPart.Elements()) {
            if (r_elem.Is(SELECTED)) {
                const auto& r_geom = r_elem.GetGeometry();
                const auto& r_N_container = r_geom.ShapeFunctionsValues(integration_method);

                for (std::size_t i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
                    const auto& r_N = row(r_N_container, i_gauss);

                    // Interpolate the current particle coordinates and velocity
                    // Note that we use the velocity in the first buffer position as this is the une used in the level set convection
                    aux_v = ZeroVector(3);
                    aux_coords = ZeroVector(3);
                    for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
                        noalias(aux_coords) += r_N[i_node] * r_geom[i_node].Coordinates();
                        noalias(aux_v) += r_N[i_node] * r_geom[i_node].FastGetSolutionStepValue(VELOCITY);
                    }

                    // Calculate the current particle final position
                    noalias(aux_coords) += dt * aux_v;

                    // Save data in the particle data container
                    auto p_new_particle = Kratos::make_shared<ParticleDataType>(aux_coords, aux_v);
                    particle_data.push_back(p_new_particle);
                }
            }
        }

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

        // Compute the Lagrangian velocity interpolation
        const std::size_t max_results = 50;
        for (auto& r_node : rModelPart.Nodes()) {
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
                // Set the coordinates matrix for the MLS interpolation
                Matrix cloud_coords(n_cloud_points, TDim);
                for (std::size_t i = 0; i < n_cloud_points; ++i) {
                    const auto& r_particle = *(cloud_particles[i]);
                    for (std::size_t d = 0; d < TDim; ++d) {
                        cloud_coords(i,d) = r_particle[d];
                    }
                }

                // Calculate the MLS interpolation in the FREE_SURFACE node
                // TODO: Use an auxiliary lambda to avoid this if
                Vector cloud_shape_functions;
                RadialBasisFunctionsUtility::CalculateShapeFunctionsAndShapeParameter(
                    cloud_coords, r_node_coords, cloud_shape_functions);

                // Interpolate the FREE_SURFACE node history and mesh velocity
                array_1d<double,3> v_n = ZeroVector(3);
                for (std::size_t i = 0; i < n_cloud_points; ++i) {
                    const auto& r_particle = *(cloud_particles[i]);
                    noalias(v_n) += cloud_shape_functions(i) * r_particle.OldVelocity;
                }

                // In order to avoid changing abruptly mesh velocity between adjacent nodes, a progessive change has been done weighting old velocity and new predicted velocity along all the area defined by the user.
                // The weights are linearly dependant with the ParticleLayerThickness.Close to the area where the nodes have changes its density, new predictd velocity will have more weight.
                double distance_ratio = 0.0;
                double distance=r_node.FastGetSolutionStepValue(DISTANCE);
                double old_distance=r_node.FastGetSolutionStepValue(DISTANCE,1);

                if (!r_node.Is(FREE_SURFACE)){

                    if (std::abs(distance)<ParticleLayerThickness){

                        if ((distance-old_distance)<0.0){
                            // the level set has positive velocity (up) / weighted air area /(phi_n+1)<-----(phi_n)/weighted water area/
                            if (distance < 0.0 && std::abs(old_distance)<ParticleLayerThickness ){
                                // Slope is calculated : Water in this case is limited by the previous free surface as an upper band
                                distance_ratio=std::abs(old_distance)/ParticleLayerThickness;
                            }
                            else{
                                //Slope is calculated : Air in this case is limited by the actual free surface as an lower band
                                 distance_ratio=std::abs(distance)/ParticleLayerThickness;
                            }
                        }
                        else{
                            // the level set has negative velocity (down)
                           //  / weighted air area /(phi_n)----->(phi_n+1)/weighted water area/
                            if (distance < 0.0 ){
                                // Slope is calculated : Water in this case is limited by the actual free surface as an upper band
                                distance_ratio=std::abs(distance)/ParticleLayerThickness;
                            }
                            else if (old_distance<ParticleLayerThickness){
                                //Slope is calculated : Air in this case is limited by the previous free surface  as an lower band
                                 distance_ratio=std::abs(old_distance)/ParticleLayerThickness;
                            }
                        }
                    }
                }
                else {
                    // For free surface nodes the old step velocity will exactly the velocity predicted using the particle based fm-ale
                    distance_ratio = 1.0;
                }

                const array_1d<double, 3> original_velocity=r_node.FastGetSolutionStepValue(VELOCITY);
                r_node.SetValue(TEMPERATURE,distance_ratio);
                // Weighted velocity is calculated.
                v_n= ((1.0 - distance_ratio)*v_n + distance_ratio*original_velocity);

                // Set new interpolated velocity such that the interface is Lagrangian in an approximate sense.
                // It should be considered the fixity. For that nodes that ones of each components is fixed it is not applied false fm-ale to that component.Same procedure has been  for SLIP CONDITION,

                // double simulation_time= rModelPart.GetProcessInfo()[TIME];

                if (r_node.IsNot(SLIP)) {
                    if (!r_node.IsFixed(VELOCITY_X)){
                        r_node.FastGetSolutionStepValue(VELOCITY_X) = v_n[0];
                        r_node.FastGetSolutionStepValue(VELOCITY_X, 1) = v_n[0];
                        r_node.FastGetSolutionStepValue(MESH_VELOCITY_X) = v_n[0];
                    }

                    if (!r_node.IsFixed(VELOCITY_Y)){
                        r_node.FastGetSolutionStepValue(VELOCITY_Y) = v_n[1];
                        r_node.FastGetSolutionStepValue(VELOCITY_Y, 1) = v_n[1];
                        r_node.FastGetSolutionStepValue(MESH_VELOCITY_Y) = v_n[1];
                    }

                    if (!r_node.IsFixed(VELOCITY_Z)){
                        r_node.FastGetSolutionStepValue(VELOCITY_Z) = v_n[2];
                        r_node.FastGetSolutionStepValue(VELOCITY_Z, 1) = v_n[2];
                        r_node.FastGetSolutionStepValue(MESH_VELOCITY_Z) = v_n[2];
                    }
                } else {
                    const auto& r_normal = r_node.FastGetSolutionStepValue(NORMAL);
                    const array_1d<double,3> v_n_norm = inner_prod(v_n, r_normal) * r_normal;
                    const array_1d<double,3> v_n_tang = v_n - v_n_norm;
                    r_node.FastGetSolutionStepValue(VELOCITY) = v_n_tang;
                    r_node.FastGetSolutionStepValue(VELOCITY, 1) = v_n_tang;
                    r_node.FastGetSolutionStepValue(MESH_VELOCITY) = v_n_tang;
                }

            }
        }
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

