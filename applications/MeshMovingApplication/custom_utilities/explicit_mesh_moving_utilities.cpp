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
        mrVirtualModelPart(rModelPart),
        mrStructureModelPart(rStructureModelPart)
    {
        if (mrStructureModelPart.GetBufferSize() < 2) {
            (mrStructureModelPart.GetRootModelPart()).SetBufferSize(2);
            KRATOS_WARNING("ExplicitMeshMovingUtilities") << "Structure model part buffer size is 1. Setting buffer size to 2." << std::endl;
        }
    }

    void ExplicitMeshMovingUtilities::FillVirtualModelPart(ModelPart& rOriginModelPart){

        // Check that the origin model part has nodes and elements to be copied
        KRATOS_ERROR_IF(rOriginModelPart.NumberOfNodes() == 0) << "Origin model part has no nodes.";
        KRATOS_ERROR_IF(rOriginModelPart.NumberOfElements() == 0) << "Origin model part has no elements.";

        // Set the buffer size in the virtual model part
        mrVirtualModelPart.SetBufferSize(rOriginModelPart.GetBufferSize());

        // Copy the origin model part nodes and fixity
        auto &r_nodes_array = rOriginModelPart.NodesArray(); 
        for(auto it_node : r_nodes_array){
            // Create a copy of the origin model part node
            auto p_node = mrVirtualModelPart.CreateNewNode(it_node->Id(),*it_node, 0);
            // Check fixity
            if (it_node->IsFixed(MESH_DISPLACEMENT_X))
                p_node->Fix(MESH_DISPLACEMENT_X);
            if (it_node->IsFixed(MESH_DISPLACEMENT_Y))
                p_node->Fix(MESH_DISPLACEMENT_Y);
            if (it_node->IsFixed(MESH_DISPLACEMENT_Z))
                p_node->Fix(MESH_DISPLACEMENT_Z);
        }

        // Copy the origin model part elements
        auto &r_elems = rOriginModelPart.Elements();
        for(auto &elem : r_elems){
            // Set the array of virtual nodes to create the element from the original ids.
            PointsArrayType nodes_array;
            auto &r_orig_geom = elem.GetGeometry();
            for (unsigned int i = 0; i < r_orig_geom.PointsNumber(); ++i){
                nodes_array.push_back(mrVirtualModelPart.pGetNode(r_orig_geom[i].Id()));
            }

            // Create the same element but using the virtual model part nodes
            auto p_elem = elem.Create(elem.Id(), nodes_array, elem.pGetProperties());
            mrVirtualModelPart.AddElement(p_elem);
        }

        // Check that the nodes and elements have been correctly copied
        KRATOS_ERROR_IF(rOriginModelPart.NumberOfNodes() != mrVirtualModelPart.NumberOfNodes()) 
            << "Origin and virtual model part have different number of nodes.";
        KRATOS_ERROR_IF(rOriginModelPart.NumberOfElements() != mrVirtualModelPart.NumberOfElements()) 
            << "Origin and virtual model part have different number of elements.";
    }

    void ExplicitMeshMovingUtilities::ComputeExplicitMeshMovement(const double DeltaTime){

        const int time_order = 1;
        VectorResultNodesContainerType search_results;
        DistanceVectorContainerType search_distance_results;

        SearchStructureNodes(search_results, search_distance_results);
        ComputeMeshDisplacement(search_results, search_distance_results);
        MoveMeshUtilities::CalculateMeshVelocities(mrVirtualModelPart, time_order, DeltaTime);
        MoveMeshUtilities::MoveMesh(mrVirtualModelPart.Nodes());

        // Check that the moved virtual mesh has no negative Jacobian elements
        for (auto it_elem : mrVirtualModelPart.ElementsArray())
            KRATOS_ERROR_IF((it_elem->GetGeometry()).Area() < 0.0) << "Element " << it_elem->Id() << " in virtual model part has negative jacobian." << std::endl;
    }

    void ExplicitMeshMovingUtilities::UndoMeshMovement(){

        auto &r_nodes = mrVirtualModelPart.Nodes();
        MoveMeshUtilities::SetMeshToInitialConfiguration(r_nodes);
    }

    template <unsigned int TDim>
    void ExplicitMeshMovingUtilities::ProjectVirtualValues(
        ModelPart& rOriginModelPart,
        unsigned int BufferSize){

        // Check that the virtual model part has elements
        KRATOS_ERROR_IF(mrVirtualModelPart.NumberOfElements() == 0) << "Virtual model part has no elements."; 

        // Set the binbased fast point locator utility
        BinBasedFastPointLocator<TDim> bin_based_point_locator(mrVirtualModelPart);
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
                it_node->GetSolutionStepValue(MESH_VELOCITY) = ZeroVector(3);
                for (unsigned int i_step = 0; i_step < BufferSize; ++i_step){
                    it_node->GetSolutionStepValue(PRESSURE, i_step) = 0.0;
                    it_node->GetSolutionStepValue(VELOCITY, i_step) = ZeroVector(3);
                }

                // Interpolate the origin model part node value
                auto &r_geom = p_elem->GetGeometry();
                for (std::size_t i_virt_node = 0; i_virt_node < r_geom.PointsNumber(); ++i_virt_node){
                    it_node->GetSolutionStepValue(MESH_VELOCITY) += aux_N(i_virt_node) * r_geom[i_virt_node].GetSolutionStepValue(MESH_VELOCITY);
                    for (unsigned int i_step = 0; i_step < BufferSize; ++i_step){
                        it_node->GetSolutionStepValue(PRESSURE, i_step) += aux_N(i_virt_node) * r_geom[i_virt_node].GetSolutionStepValue(PRESSURE, i_step);
                        it_node->GetSolutionStepValue(VELOCITY, i_step) += aux_N(i_virt_node) * r_geom[i_virt_node].GetSolutionStepValue(VELOCITY, i_step);
                    }
                }
            } else {
                KRATOS_WARNING("ExplicitMeshMovingUtility") << "Origin model part node " << it_node->Id() << " has not been found in any virtual model part element." << std::endl;
            }
        }
    }

    /* Private functions *******************************************************/

    void ExplicitMeshMovingUtilities::SearchStructureNodes(
        VectorResultNodesContainerType &rSearchResults,
        DistanceVectorContainerType &rSearchDistanceResults) {

        auto &r_virt_nodes = mrVirtualModelPart.NodesArray();
        const unsigned int n_virt_nodes = r_virt_nodes.size();
        auto &r_structure_nodes = mrStructureModelPart.NodesArray();
        const unsigned int max_number_of_nodes = r_structure_nodes.size(); 

        rSearchResults.resize(n_virt_nodes);
        rSearchDistanceResults.resize(n_virt_nodes);

        NodeBinsType bins(r_structure_nodes.begin(), r_structure_nodes.end()); 

        #pragma omp parallel for
        for(int i_fl = 0; i_fl < static_cast<int>(n_virt_nodes); ++i_fl){
            // Initialize the search results variables
            std::size_t n_results = 0;
            ResultNodesContainerType i_node_results(max_number_of_nodes);
            DistanceVectorType i_node_distance_results(max_number_of_nodes);
            auto it_results = i_node_results.begin();
            auto it_results_distances = i_node_distance_results.begin();
            
            // Perform the structure nodes search
            auto it_virt_node = r_virt_nodes.begin() + i_fl;
            n_results = bins.SearchObjectsInRadiusExclusive(*it_virt_node, mSearchRadius, it_results, it_results_distances, max_number_of_nodes);

            // Resize and save the current virtual node results
            i_node_results.resize(n_results);
            i_node_distance_results.resize(n_results);

            rSearchResults[i_fl] = i_node_results;
            rSearchDistanceResults[i_fl] = i_node_distance_results;
        }
    }

    void ExplicitMeshMovingUtilities::ComputeMeshDisplacement(
        const VectorResultNodesContainerType &rSearchResults,
        const DistanceVectorContainerType &rSearchDistanceResults){

        #pragma omp parallel for
        for(int i_fl = 0; i_fl < static_cast<int>(mrVirtualModelPart.NumberOfNodes()); ++i_fl){
            // Get auxiliar current fluid node info.
            auto it_node = mrVirtualModelPart.NodesBegin() + i_fl;
            const auto i_fl_str_nodes = rSearchResults[i_fl];
            const auto i_fl_str_dists = rSearchDistanceResults[i_fl];
            const std::size_t n_str_nodes = i_fl_str_nodes.size();

            // Initialize the current virtal model part node MESH_DISPLACEMENT
            auto &r_mesh_disp = it_node->FastGetSolutionStepValue(MESH_DISPLACEMENT);
            r_mesh_disp = ZeroVector(3);

            // Check if any structure node is found
            if (n_str_nodes != 0){
                // Compute the average MESH_DISPLACEMENT
                for(unsigned int i_str = 0; i_str < n_str_nodes; ++i_str){
                    // Compute the structure point weight according to the kernel function
                    const double normalised_distance = i_fl_str_dists[i_str] / mSearchRadius;
                    const double weight = this->ComputeKernelValue(normalised_distance);

                    // Accumulate the current step structure pt. DISPLACEMENT values
                    const auto &r_str_disp_0 = (*(i_fl_str_nodes[i_str])).FastGetSolutionStepValue(DISPLACEMENT,0);
                    const auto &r_str_disp_1 = (*(i_fl_str_nodes[i_str])).FastGetSolutionStepValue(DISPLACEMENT,1);
                    const auto str_disp = r_str_disp_0 - r_str_disp_1;
                    if (!it_node->IsFixed(MESH_DISPLACEMENT_X))
                        r_mesh_disp[0] += weight * str_disp[0];
                    if (!it_node->IsFixed(MESH_DISPLACEMENT_Y))
                        r_mesh_disp[1] += weight * str_disp[1];
                    if (!it_node->IsFixed(MESH_DISPLACEMENT_Z))
                        r_mesh_disp[2] += weight * str_disp[2];
                }

                r_mesh_disp /= n_str_nodes;
            }
        }
    }

    inline double ExplicitMeshMovingUtilities::ComputeKernelValue(const double NormalisedDistance){
        // Epanechnikov (parabolic) kernel function
        return (std::abs(NormalisedDistance) <= 1.0) ? std::abs((3.0/4.0)*(1.0-std::pow(NormalisedDistance,2))) : 0.0;
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
