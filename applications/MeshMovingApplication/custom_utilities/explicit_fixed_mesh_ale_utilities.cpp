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
//
//

// System includes

// External includes

// Application includes
#include "custom_utilities/explicit_fixed_mesh_ale_utilities.h"
#include "custom_utilities/mesh_velocity_calculation.h"
#include "custom_utilities/move_mesh_utilities.h"

// Project includes
#include "includes/mesh_moving_variables.h"
#include "spatial_containers/spatial_containers.h"
#include "utilities/binbased_fast_point_locator.h"
#include "utilities/spatial_containers_configure.h"
#include "utilities/variable_utils.h"
#include "utilities/parallel_utilities.h"


namespace Kratos
{
/* Public functions *******************************************************/

ExplicitFixedMeshALEUtilities::ExplicitFixedMeshALEUtilities(
    ModelPart &rVirtualModelPart,
    ModelPart &rStructureModelPart,
    const double SearchRadius) :
    FixedMeshALEUtilities(
        rVirtualModelPart,
        rStructureModelPart),
    mSearchRadius(SearchRadius)
{
    // Check the structure model part
    if (mrStructureModelPart.GetBufferSize() < 2) {
        (mrStructureModelPart.GetRootModelPart()).SetBufferSize(2);
        KRATOS_WARNING("ExplicitFixedMeshALEUtilities") << "Structure model part buffer size is 1. Setting buffer size to 2." << std::endl;
    }
}

ExplicitFixedMeshALEUtilities::ExplicitFixedMeshALEUtilities(
    Model &rModel,
    Parameters &rParameters) :
    FixedMeshALEUtilities(
        rModel.GetModelPart(rParameters["virtual_model_part_name"].GetString()),
        rModel.GetModelPart(rParameters["structure_model_part_name"].GetString())),
    mSearchRadius(rParameters["search_radius"].GetDouble())
{
    // Validate with default parameters
    Parameters default_parameters(R"(
    {
        "virtual_model_part_name": "",
        "structure_model_part_name": "",
        "search_radius": 0.0
    }  )");
    rParameters.ValidateAndAssignDefaults(default_parameters);

    // Check the input level set type
    if (mSearchRadius <= 0.0) {
        KRATOS_ERROR << "Search radius  is: " << mSearchRadius << ". A value larger than 0 is expected.";
    }

    // Check the structure model part
    if (mrStructureModelPart.GetBufferSize() < 2) {
        (mrStructureModelPart.GetRootModelPart()).SetBufferSize(2);
        KRATOS_WARNING("FixedMeshALEUtilities") << "Structure model part buffer size is 1. Setting buffer size to 2." << std::endl;
    }
}

void ExplicitFixedMeshALEUtilities::Initialize(ModelPart &rOriginModelPart)
{
    // Fill the virtual model part as a copy of the origin mesh
    this->FillVirtualModelPart(rOriginModelPart);
}

void
ExplicitFixedMeshALEUtilities::ComputeMeshMovement(const double DeltaTime)
{
    VectorResultNodesContainerType search_results;
    DistanceVectorContainerType search_distance_results;

    this->SearchStructureNodes(search_results, search_distance_results);
    this->ComputeExplicitMeshDisplacement(search_results, search_distance_results);
    TimeDiscretization::BDF1 time_disc_BDF1;
    mrVirtualModelPart.GetProcessInfo()[DELTA_TIME] = DeltaTime;
    MeshVelocityCalculation::CalculateMeshVelocities(mrVirtualModelPart, time_disc_BDF1);
    MoveMeshUtilities::MoveMesh(mrVirtualModelPart.Nodes());

    // Check that the moved virtual mesh has no negative Jacobian elements
#ifdef KRATOS_DEBUG
    for (auto it_elem : mrVirtualModelPart.ElementsArray()) {
        KRATOS_ERROR_IF((it_elem->GetGeometry()).Area() < 0.0) << "Element " << it_elem->Id() << " in virtual model part has negative jacobian." << std::endl;
    }
#endif
}

/* Private functions *******************************************************/

void ExplicitFixedMeshALEUtilities::SearchStructureNodes(
    VectorResultNodesContainerType &rSearchResults,
    DistanceVectorContainerType &rSearchDistanceResults) {

    auto &r_virt_nodes = mrVirtualModelPart.NodesArray();
    const unsigned int n_virt_nodes = r_virt_nodes.size();
    auto &r_structure_nodes = mrStructureModelPart.NodesArray();
    const unsigned int max_number_of_nodes = r_structure_nodes.size();

    rSearchResults.resize(n_virt_nodes);
    rSearchDistanceResults.resize(n_virt_nodes);

    NodeBinsType bins(r_structure_nodes.begin(), r_structure_nodes.end());

    IndexPartition<size_t>( static_cast<size_t>(n_virt_nodes) ).for_each(
    [&]( size_t index )
    {
        // Initialize the search results variables
        std::size_t n_results = 0;
        ResultNodesContainerType i_node_results(max_number_of_nodes);
        DistanceVectorType i_node_distance_results(max_number_of_nodes);
        auto it_results = i_node_results.begin();
        auto it_results_distances = i_node_distance_results.begin();

        // Perform the structure nodes search
        auto it_virt_node = r_virt_nodes.begin() + index;
        n_results = bins.SearchObjectsInRadiusExclusive(*it_virt_node, mSearchRadius, it_results, it_results_distances, max_number_of_nodes);

        // Resize and save the current virtual node results
        i_node_results.resize(n_results);
        i_node_distance_results.resize(n_results);

        rSearchResults[index] = i_node_results;
        rSearchDistanceResults[index] = i_node_distance_results;
    } );
}

void ExplicitFixedMeshALEUtilities::ComputeExplicitMeshDisplacement(
    const VectorResultNodesContainerType &rSearchResults,
    const DistanceVectorContainerType &rSearchDistanceResults)
{
    IndexPartition<std::size_t>( static_cast<std::size_t>(mrVirtualModelPart.NumberOfNodes()) ).for_each(
    [&]( std::size_t i_fl )
    {
        // Get auxiliar current fluid node info.
        auto it_node = mrVirtualModelPart.NodesBegin() + i_fl;
        const auto i_fl_str_nodes = rSearchResults[i_fl];
        const auto i_fl_str_dists = rSearchDistanceResults[i_fl];
        const std::size_t n_str_nodes = i_fl_str_nodes.size();

        // Check origin model part mesh displacementfixity
        const auto it_orig_node = mpOriginModelPart->NodesBegin() + i_fl;
        if (it_orig_node->IsFixed(MESH_DISPLACEMENT_X)) {
            it_node->Fix(MESH_DISPLACEMENT_X);
        }
        if (it_orig_node->IsFixed(MESH_DISPLACEMENT_Y)) {
            it_node->Fix(MESH_DISPLACEMENT_Y);
        }
        if (it_orig_node->IsFixed(MESH_DISPLACEMENT_Z)) {
            it_node->Fix(MESH_DISPLACEMENT_Z);
        }

        // Initialize the current virtual model part node MESH_DISPLACEMENT
        auto &r_mesh_disp = it_node->FastGetSolutionStepValue(MESH_DISPLACEMENT);
        r_mesh_disp = ZeroVector(3);

        // Check if any structure node is found
        if (n_str_nodes != 0){
            // Compute the closest structure node MESH_DISPLACEMENT
            double min_distance = 1.01;
            Node<3>::Pointer str_node_ptr = nullptr;
            for(unsigned int i_str = 0; i_str < n_str_nodes; ++i_str) {
                // Compute the structure point weight according to the kernel function
                const double normalised_distance = std::sqrt(i_fl_str_dists[i_str]) / mSearchRadius;

                // Check if the point is closer than the current one
                if (normalised_distance < min_distance) {
                    min_distance = normalised_distance;
                    str_node_ptr = i_fl_str_nodes[i_str];
                }
            }

            // Current step structure pt. DISPLACEMENT values
            if (str_node_ptr) {
                const double weight = this->ComputeKernelValue(min_distance);
                const auto &r_str_disp_0 = str_node_ptr->FastGetSolutionStepValue(DISPLACEMENT,0);
                const auto &r_str_disp_1 = str_node_ptr->FastGetSolutionStepValue(DISPLACEMENT,1);
                const auto str_disp = weight * (r_str_disp_0 - r_str_disp_1);
                if (!it_node->IsFixed(MESH_DISPLACEMENT_X)) {
                    r_mesh_disp[0] = str_disp[0];
                }
                if (!it_node->IsFixed(MESH_DISPLACEMENT_Y)) {
                    r_mesh_disp[1] = str_disp[1];
                }
                if (!it_node->IsFixed(MESH_DISPLACEMENT_Z)) {
                    r_mesh_disp[2] = str_disp[2];
                }
            } // if ( str_node_ptr )
        } // if ( n_str_nodes != 0 )
    } ); // IndexPartition.for_each
}

inline double ExplicitFixedMeshALEUtilities::ComputeKernelValue(const double NormalisedDistance){
    // Epanechnikov (parabolic) kernel function
    // return (std::abs(NormalisedDistance) < 1.0) ? std::abs((3.0/4.0)*(1.0-std::pow(NormalisedDistance,2))) : 0.0;
    // Triangle kernel function
    return (std::abs(NormalisedDistance) < 1.0) ? 1.0 - std::abs(NormalisedDistance) : 0.0;
}

void ExplicitFixedMeshALEUtilities::CreateVirtualModelPartElements(const ModelPart &rOriginModelPart)
{
    // Copy the origin model part elements
    auto &r_elems = rOriginModelPart.Elements();
    for(auto &elem : r_elems) {
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
}

/* External functions *****************************************************/

/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const ExplicitFixedMeshALEUtilities& rThis) {

    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos
