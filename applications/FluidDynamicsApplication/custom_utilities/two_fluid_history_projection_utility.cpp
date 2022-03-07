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
#include "processes/find_nodal_h_process.h"
#include "utilities/mls_shape_functions_utility.h"
#include "utilities/variable_utils.h"

// Application includes
#include "two_fluid_history_projection_utility.h"

namespace Kratos
{
    /* Public functions *******************************************************/

    static void TwoFluidHistoryProjectionUtility::CalculateHistoryProjection(
        ModelPart& rModelPart,
        const bool ComputeNodalH)
    {
        const std::size_t n_extra_layers = 5;
        FlagElementsAndNodes(rModelPart, n_extra_layers);

        const auto particle_data = SeedAndConvectParticles(rModelPart);

        if (ComputeNodalH) {
            FindNodalHProcess<FindNodalHSettings::SaveAsNonHistoricalVariable> nodal_h_calculator(rModelPart);
            nodal_h_calculator.Execute();
        }

        CalculateLagrangianVelocityInterpolation(rModelPart, particle_data);
    }

    /* Private functions *******************************************************/

    static void TwoFluidHistoryProjectionUtility::FlagElementsAndNodes(
        ModelPart& rModelPart,
        const std::size_t NumberOfExtraLayers)
    {
        // Reset flags
        // SELECTED: Nodes and elements used for creating the seeds
        // FREE_SURFACE: Nodes and elements where the history is to be reset
        VariableUtils().SetFlag(SELECTED, false, rModelPart.Nodes());
        VariableUtils().SetFlag(SELECTED, false, rModelPart.Elements());
        VariableUtils().SetFlag(FREE_SURFACE, false, rModelPart.Nodes());
        VariableUtils().SetFlag(FREE_SURFACE, false, rModelPart.Elements());

        // Flag the nodes that have changed its phase
        block_for_each(rModelPart.Nodes(), [](Node<3>& rNode){
            const double new_d = rNode.FastGetSolutionStepValue(DISTANCE);
            const double old_d = rNode.FastGetSolutionStepValue(DISTANCE,1);
            if (new_d * old_d < 0.0) {
                rNode.Set(SELECTED, true);
                rNode.Set(FREE_SURFACE, true);
            }
        });

        // Extra layers for seeding
        // Note that in here we are assuming that the neighbours are already computed
        for (std::size_t i_layer = 0; i_layer < NumberOfExtraLayers; ++i_layer) {
            for (const auto& r_node : rModelPart.Nodes()) {
                auto& r_neighbours = r_node.GetValue(NEIGHBOUR_NODES);
                for (auto& r_neigh : r_neighbours) {
                    if (!r_neigh.Is(SELECTED)) {
                        r_neight.Set(SELECTED, true);
                    }
                }
            }
        }

        // Flag all the elements that own a flagged node
        block_for_each(rModelPart.Elements(), [](Element& rElement){
            const auto& r_geom = rElement.GetGeometry();
            for (const auto& rNode : r_geom) {
                if (rNode.Is(FREE_SURFACE)) {
                    rElement.Set(FREE_SURFACE, true);
                }
                if (rNode.Is(SELECTED)) {
                    rElement.Set(SELECTED, true);
                }
            }
        });
    }

    static ParticleDataContainerType TwoFluidHistoryProjectionUtility::SeedAndConvectParticles(ModelPart& rModelPart)
    {
        // Initialize required data
        const double dt = rModelPart.GetProcessInfo()[DELTA_TIME]
        KRATOS_ERROR_IF(dt < 1.0e-12) << "Found DELTA_TIME close to zero in ProcessInfo." << std::endl;

        // Allocate auxiliary arrays
        ParticleDataContainerType particle_data;

        // Get and allocate arrays according to 1st element
        // Note that we are assuming a unique element type in the mesh
        const auto integration_method = GeometryData::IntegrationMethod::GI_GAUSS_2;
        const auto& elem_begin_geom = rModelPart.BeginElements()->GetGeometry();
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
                    particle_data.push_back(std::make_pair(aux_coords, aux_v));
                }
            }
        }

        return particle_data;
    }

    static void TwoFluidHistoryProjectionUtility::CalculateLagrangianVelocityInterpolation(
        ModelPart& rModelPart,
        const ParticleDataContainerType& rParticleData)
    {
        std::size_t domain_size = rModelPart.GetProcessInfo()[DOMAIN_SIZE];

        const double search_factor = 2.0;
        for (auto& r_node : rModelPart.Nodes()) {
            if (r_node.Is(FREE_SURFACE)) {
                const auto& r_node_coords = r_node.Coordinates();

                // Set the cloud of points inside the search bounding box
                std::vector<array_1d<double,3>> vel_cloud;
                std::vector<array_1d<double,3>> coord_cloud;
                const double search_dist = search_factor * r_node.GetValue(NODAL_H) / 2.0;
                for (const auto& r_data : rParticleData) {
                    const auto& r_part_coords = std::get<0>(r_data);
                    if (std::abs(r_part_coords[0] - r_node_coords[0]) < search_dist) {
                        if (std::abs(r_part_coords[1] - r_node_coords[1]) < search_dist) {
                            if (std::abs(r_part_coords[2] - r_node_coords[2]) < search_dist) {
                                coord_cloud.push_back(r_part_coords);
                                vel_cloud.push_back(std::get<1>(r_data));
                            }
                        }
                    }
                }

                // Set the coordinates matrix for the MLS interpolation
                std::size_t n_cloud_points = coord_cloud.size();
                Matrix cloud_coords(n_cloud_points, domain_size);
                for (std::size_t i = 0; i < n_cloud_points; ++i) {
                    const auto& r_aux_coords = coord_cloud[i];
                    for (std::size_t d = 0; d < domain_size; ++d) {
                        cloud_coords(i,d) = r_aux_coords[d];
                    }
                }

                // Calculate the MLS interpolation in the FREE_SURFACE node
                // TODO: Use an auxiliary lambda to avoid this if
                Vector cloud_shape_functions;
                if (domain_size == 2) {
                    MLSShapeFunctionsUtility::CalculateShapeFunctions<2>(cloud_coords, r_node_coords, search_dist, cloud_shape_functions);
                } else if (domain_size == 3) {
                    MLSShapeFunctionsUtility::CalculateShapeFunctions<3>(cloud_coords, r_node_coords, search_dist, cloud_shape_functions);
                } else {
                    KRATOS_ERROR << "Wrong DOMAIN_SIZE." << std::endl;
                }

                // Interpolate the FREE_SURFACE node history and mesh velocity
                array_1d<double,3> v_n = ZeroVector(3);
                for (std::size_t i_node = 0; i_node < n_cloud_points; ++i_node) {
                    noalias(v_n) += cloud_shape_functions(i_node) * vel_cloud[i_node];
                }

                // Set new interpolated velocity such that the interface is Lagrangian in an approximate sense
                r_node.FastGetSolutionStepValue(VELOCITY) = v_n;
                r_node.FastGetSolutionStepValue(VELOCITY, 1) = v_n;
                r_node.FastGetSolutionStepValue(MESH_VELOCITY) = v_n;
            }
        }
    }

    /* External functions *****************************************************/

    /// output stream function
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const TwoFluidHistoryProjectionUtility& rThis)
    {
        rThis.PrintData(rOStream);
        return rOStream;
    }

}
