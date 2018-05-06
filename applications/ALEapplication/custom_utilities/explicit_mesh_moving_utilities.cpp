//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//
//  Main authors:    Ruben Zorrilla
//
//					 Kratos default license: kratos/license.txt
//

// System includes

// External includes

// Application includes
#include "ale_application.h"
#include "move_mesh_utilities.h"
#include "explicit_mesh_moving_utilities.h"

// Project includes
#include "spatial_containers/spatial_containers.h"
#include "utilities/binbased_fast_point_locator.h"
#include "utilities/spatial_containers_configure.h"

namespace Kratos
{
    /* Public functions *******************************************************/

    ExplicitMeshMovingUtilities::ExplicitMeshMovingUtilities(
        ModelPart &rModelPart,
        ModelPart &rStructureModelPart,
        const double SearchRadius) : 
        mSearchRadius(SearchRadius),
        mrModelPart(rModelPart),
        mrStructureModelPart(rStructureModelPart){
    }

    void ExplicitMeshMovingUtilities::FillVirtualModelPart(ModelPart& rOriginModelPart){

        // Check that the origin model part has nodes and elements to be copied
        KRATOS_ERROR_IF(rOriginModelPart.NumberOfNodes() == 0) << "Origin model part has no nodes.";
        KRATOS_ERROR_IF(rOriginModelPart.NumberOfElements() == 0) << "Origin model part has no elements.";

        // Set the buffer size in the virtual model part
        mrModelPart.SetBufferSize(rOriginModelPart.GetBufferSize());

        // Copy the origin model part nodes
        auto &r_nodes_array = rOriginModelPart.NodesArray(); 
        for(auto node : r_nodes_array){
            auto p_node = mrModelPart.CreateNewNode(node->Id(),*node, 0);
        }

        // Copy the origin model part elements
        auto &r_elems = rOriginModelPart.Elements();
        for(auto elem : r_elems){
            // Set the array of virtual nodes to create the element from the original ids.
            PointsArrayType nodes_array;
            auto &r_orig_geom = elem.GetGeometry();
            for (unsigned int i = 0; i < r_orig_geom.PointsNumber(); ++i){
                nodes_array.push_back(mrModelPart.pGetNode(r_orig_geom[i].Id()));
            }

            // Create the same element but using the virtual model part nodes
            auto p_elem = elem.Create(elem.Id(), nodes_array, elem.pGetProperties());
            mrModelPart.AddElement(p_elem);
        }

        // Check that the nodes and elements have been correctly copied
        KRATOS_ERROR_IF(rOriginModelPart.NumberOfNodes() != mrModelPart.NumberOfNodes()) 
            << "Origin and virtual model part have different number of nodes.";
        KRATOS_ERROR_IF(rOriginModelPart.NumberOfElements() != mrModelPart.NumberOfElements()) 
            << "Origin and virtual model part have different number of elements.";
    }

    void ExplicitMeshMovingUtilities::ComputeExplicitMeshMovement(const double DeltaTime){

        const int time_order = 1;
        auto &r_nodes = mrModelPart.Nodes();
        VectorResultNodesContainerType search_results;
        DistanceVectorContainerType search_distance_results;

        SearchStructureNodes(search_results, search_distance_results);
        ComputeMeshDisplacement(search_results, search_distance_results);
        MoveMeshUtilities::CalculateMeshVelocities(mrModelPart, time_order, DeltaTime);
        MoveMeshUtilities::MoveMesh(r_nodes);
    }

    void ExplicitMeshMovingUtilities::UndoMeshMovement(){

        auto &r_nodes = mrModelPart.Nodes();
        MoveMeshUtilities::SetMeshToInitialConfiguration(r_nodes);
    }

    template <unsigned int TDim>
    void ExplicitMeshMovingUtilities::ProjectVirtualValues(
        ModelPart& rOriginModelPart,
        unsigned int BufferSize){

        // Check that the virtual model part has elements
        KRATOS_ERROR_IF(mrModelPart.NumberOfElements() == 0) << "Virtual model part has no elements."; 

        // Set the binbased fast point locator utility
        BinBasedFastPointLocator<TDim> bin_based_point_locator(mrModelPart);
        bin_based_point_locator.UpdateSearchDatabase();

        // Search the origin model part nodes in the virtual mesh elements and 
        // interpolate the values in the virtual element to the origin model part node
        auto &r_nodes_array = rOriginModelPart.NodesArray(); 
        for(auto it_node : r_nodes_array){
            // Find the origin model part node in the virtual mesh
            Vector aux_N;
            Element::Pointer p_elem;
            const bool is_found = bin_based_point_locator.FindPointOnMeshSimplified(it_node->Coordinates(), aux_N, p_elem);

            // Check if the node is found
            if (is_found){
                // Initialize historical data
                // The current step values are also set as a prediction
                auto &r_geom = p_elem->GetGeometry();
                it_node->GetSolutionStepValue(MESH_VELOCITY) = ZeroVector(3);
                for (unsigned int i_step = 0; i_step < BufferSize; ++i_step){
                    it_node->GetSolutionStepValue(PRESSURE, i_step) = 0.0;
                    it_node->GetSolutionStepValue(VELOCITY, i_step) = ZeroVector(3);
                }

                // Interpolate the origin model part node value
                for (std::size_t i_node = 0; i_node < r_geom.PointsNumber(); ++i_node){
                    it_node->GetSolutionStepValue(MESH_VELOCITY) += aux_N(i_node) * r_geom[i_node].GetSolutionStepValue(MESH_VELOCITY);
                    for (unsigned int i_step = 0; i_step < BufferSize; ++i_step){
                        it_node->GetSolutionStepValue(PRESSURE, i_step) += aux_N(i_node) * r_geom[i_node].GetSolutionStepValue(PRESSURE, i_step);
                        it_node->GetSolutionStepValue(VELOCITY, i_step) += aux_N(i_node) * r_geom[i_node].GetSolutionStepValue(VELOCITY, i_step);
                    }
                }
            } else {
                KRATOS_WARNING("ExplicitMeshMovingUtility") << "Origin model part node " << it_node->Id() << " has not been found in any virtual model part element.";
            }
        }
    }

    /* Private functions *******************************************************/

    void ExplicitMeshMovingUtilities::SearchStructureNodes(
        VectorResultNodesContainerType &rSearchResults,
        DistanceVectorContainerType &rSearchDistanceResults) {

        auto &r_nodes = mrModelPart.NodesArray();
        const unsigned int n_nodes = r_nodes.size();
        auto &r_structure_nodes = mrStructureModelPart.NodesArray();
        const unsigned int max_number_of_nodes = r_structure_nodes.size(); 

        rSearchResults.resize(n_nodes);
        rSearchDistanceResults.resize(n_nodes);

        NodeBinsType bins(r_structure_nodes.begin(), r_structure_nodes.end()); 

        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(n_nodes); i++){

            std::size_t n_results = 0;
            ResultNodesContainerType i_node_results(max_number_of_nodes);
            DistanceVectorType i_node_distance_results(max_number_of_nodes);

            auto it_results = i_node_results.begin();
            auto it_results_distances = i_node_distance_results.begin();
            
            auto it_node = r_nodes.begin() + i;
            n_results = bins.SearchObjectsInRadius(*it_node, mSearchRadius, it_results, it_results_distances, max_number_of_nodes);

            rSearchResults[i] = i_node_results;
            rSearchDistanceResults[i] = i_node_distance_results;
        }
    }

    void ExplicitMeshMovingUtilities::ComputeMeshDisplacement(
        const VectorResultNodesContainerType &rSearchResults,
        const DistanceVectorContainerType &rSearchDistanceResults){

        #pragma omp parallel for
        for(int i_fl = 0; i_fl < static_cast<int>(mrModelPart.NumberOfNodes()); ++i_fl){
            // Get auxiliar current fluid node info.
            auto it_node = mrModelPart.NodesBegin() + i_fl;
            const auto i_fl_str_nodes = rSearchResults[i_fl];
            const auto i_fl_str_dists = rSearchDistanceResults[i_fl];
            const std::size_t n_str_nodes = i_fl_str_nodes.size();

            // Compute the average MESH_DISPLACEMENT
            double total_weight = 0.0;
            auto &r_mesh_disp = it_node->FastGetSolutionStepValue(MESH_DISPLACEMENT);
            r_mesh_disp = ZeroVector(3);

            for(unsigned int i_str = 0; i_str < n_str_nodes; ++i_str){
                // Compute the structure point weight according to the kernel function
                const double normalised_distance = i_fl_str_dists[i_str] / mSearchRadius;
                const double weight = this->ComputeKernelValue(normalised_distance);

                // Accumulate the current structure pt. DISPLACEMENT values
                const auto &r_str_disp = (mrStructureModelPart.NodesBegin() + i_str)->FastGetSolutionStepValue(DISPLACEMENT);
                if (!it_node->IsFixed(MESH_DISPLACEMENT_X))
                    r_mesh_disp[0] += weight * r_str_disp[0];
                if (!it_node->IsFixed(MESH_DISPLACEMENT_Y))
                    r_mesh_disp[1] += weight * r_str_disp[1];
                if (!it_node->IsFixed(MESH_DISPLACEMENT_Z))
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

    /* External functions *****************************************************/

    /// output stream function
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const ExplicitMeshMovingUtilities& rThis) {

        rThis.PrintData(rOStream);
        return rOStream;
    }

    template void ExplicitMeshMovingUtilities::ProjectVirtualValues<2>(ModelPart &rOriginModelPart, unsigned int BufferSize);
    template void ExplicitMeshMovingUtilities::ProjectVirtualValues<3>(ModelPart &rOriginModelPart, unsigned int BufferSize);
}
