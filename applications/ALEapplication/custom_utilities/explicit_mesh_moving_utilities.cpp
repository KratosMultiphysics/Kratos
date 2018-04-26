//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// System includes

// External includes

// Application includes
#include "ale_application.h"
#include "move_mesh_utilities.h"
#include "explicit_mesh_moving_utilities.h"

// Project includes
#include "spatial_containers/spatial_containers.h"
#include "utilities/spatial_containers_configure.h"

namespace Kratos
{
    /* Public functions *******************************************************/

    ExplicitMeshMovingUtilities::ExplicitMeshMovingUtilities(
        ModelPart &rModelPart,
        ModelPart &rStructureModelPart,
        const double SearchRadius) : 
        mrModelPart(rModelPart),
        mrStructureModelPart(rStructureModelPart),
        mSearchRadius(SearchRadius) {
    }

    void ExplicitMeshMovingUtilities::SearchStructureNodes(
        VectorResultNodesContainerType &rSearchResults,
        VectorDistanceType &rSearchDistanceResults) {

        auto &r_nodes = mrModelPart.Nodes();
        auto &r_structure_nodes = mrStructureModelPart.NodesArray();
        const unsigned int max_number_of_nodes = r_structure_nodes.size(); 

        NodeBinsType bins(r_structure_nodes.begin(), r_structure_nodes.end()); 

        #pragma omp parallel
        {
            std::size_t n_results = 0;
            DistanceType local_results_distances(max_number_of_nodes);
            ResultNodesContainerType local_results(max_number_of_nodes);
            
            #pragma omp for
            for(int i = 0; i < static_cast<int>(r_nodes.size()); i++)
            {
                ResultNodesContainerType::iterator it_results = local_results.begin();
                DistanceType::iterator it_results_distances = local_results_distances.begin();
                
                n_results = bins.SearchObjectsInRadius(r_nodes(i), mSearchRadius, it_results, it_results_distances, max_number_of_nodes);

                rSearchResults[i].insert(rSearchResults[i].begin(), local_results.begin(), local_results.begin() + n_results);
                rSearchDistanceResults[i].insert(rSearchDistanceResults[i].begin(), local_results_distances.begin(), local_results_distances.begin() + n_results);
            }
        }
    }

    void ExplicitMeshMovingUtilities::ComputeMeshDisplacement(
        const VectorResultNodesContainerType &rSearchResults,
        const VectorDistanceType &rSearchDistanceResults){

        auto &r_nodes = mrModelPart.Nodes();
        auto &r_structure_nodes = mrStructureModelPart.Nodes();
        
        #pragma omp for
        for(int i_fl = 0; i_fl < static_cast<int>(r_nodes.size()); ++i_fl){
            // Get auxiliar current fluid node info.
            auto &r_node = r_nodes[i_fl];
            const auto i_str_nodes = rSearchResults[i_fl];
            const auto i_str_dists = rSearchDistanceResults[i_fl];
            const std::size_t n_str_nodes = i_str_nodes.size();

            // Compute the average MESH_DISPLACEMENT
            double total_weight = 0.0;
            auto &r_mesh_disp = r_node.FastGetSolutionStepValue(MESH_DISPLACEMENT);
            r_mesh_disp = ZeroVector(3);

            for(unsigned int i_str = 0; i_str < n_str_nodes; ++i_str){
                // Compute the structure point weight according to the kernel function
                const double normalised_distance = i_str_dists[i_str] / mSearchRadius;
                const double weight = this->ComputeKernelValue(normalised_distance);

                // Accumulate the current structure pt. DISPLACEMENT values
                const auto &r_str_disp = r_structure_nodes[i_str].FastGetSolutionStepValue(DISPLACEMENT);
                if (!r_node.IsFixed(MESH_DISPLACEMENT_X))
                    r_mesh_disp[0] += weight * r_str_disp[0];
                if (!r_node.IsFixed(MESH_DISPLACEMENT_Y))
                    r_mesh_disp[1] += weight * r_str_disp[1];
                if (!r_node.IsFixed(MESH_DISPLACEMENT_Z))
                    r_mesh_disp[2] += weight * r_str_disp[2];

                // Accumulate the current structure pt. weight to compute the average
                total_weight += weight;
            }

            r_mesh_disp /= total_weight;
        }
    }

    inline double ExplicitMeshMovingUtilities::ComputeKernelValue(const double NormalisedDistance){
        // Epachenikov (parabolic) kernel function
        return (3.0/4.0)*(1.0-std::pow(NormalisedDistance,2));
    }

    void ExplicitMeshMovingUtilities::ComputeExplicitMeshMovement(const double DeltaTime){

        const int time_order = 1;
        auto &r_nodes = mrModelPart.Nodes();
        ModelPart::Pointer p_model_part = Kratos::make_shared<ModelPart>(mrModelPart);
        VectorResultNodesContainerType search_results;
        VectorDistanceType search_distance_results;

        SearchStructureNodes(search_results, search_distance_results);
        ComputeMeshDisplacement(search_results, search_distance_results);
        MoveMeshUtilities::CalculateMeshVelocities(p_model_part,time_order,DeltaTime);
        MoveMeshUtilities::MoveMesh(r_nodes);
    }

    /* External functions *****************************************************/

    /// output stream function
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const ExplicitMeshMovingUtilities& rThis) {

        rThis.PrintData(rOStream);
        return rOStream;
    }

}
