//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes
#ifdef INCLUDE_DELAUNATORCPP
    #include "delaunator.hpp"
#endif

// Project includes
#include "utilities/delaunator_utilities.h"

namespace Kratos
{
namespace DelaunatorUtilities
{
void CreateTriangleMeshFromNodes(ModelPart& rModelPart)
{
    KRATOS_TRY

    // Ensure node order
    auto& r_nodes_root_array = rModelPart.GetRootModelPart().Nodes();
    const auto it_node_root_begin = r_nodes_root_array.begin();
    #pragma omp parallel for
    for(int i=0; i<static_cast<int>(r_nodes_root_array.size()); ++i) {
        auto it_node = it_node_root_begin + i;
        it_node->SetId(i + 1);
    }

    // Getting nodes array
    const auto& r_nodes_array = rModelPart.Nodes();
    const auto it_node_begin = r_nodes_array.begin();

    // The vector to fill
    std::vector<double> coordinates;

    // NOTE: 2D asssumed

    #pragma omp parallel
    {
        std::vector<double> coordinates_buffer;

        #pragma omp for
        for(int i=0; i<static_cast<int>(r_nodes_array.size()); ++i) {
            auto it_node = it_node_begin + i;

            // Filling coordinates buffer
            coordinates_buffer.push_back(it_node->X());
            coordinates_buffer.push_back(it_node->Y());
        }

        // Combine buffers together
        #pragma omp critical
        {
            std::move(coordinates_buffer.begin(),coordinates_buffer.end(),back_inserter(coordinates));
        }
    }

#ifdef INCLUDE_MMG
    // Calling the library
    delaunator::Delaunator delaunator(coordinates);

    // Creating the triangles
    Properties::Pointer p_elem_prop = rModelPart.CreateNewProperties(0);
    const auto& r_triangles = delaunator.triangles;
    std::size_t counter = rModelPart.GetRootModelPart().Elements().size() + 1;
    for (std::size_t i = 0; i < r_triangles.size(); i += 3) {
        rModelPart.CreateNewElement("Element2D3N", counter, {{r_triangles[i] + 1,r_triangles[i + 1] + 1, r_triangles[i + 2] + 1}}, p_elem_prop);
        ++counter;
    }
#else

#endif

    KRATOS_CATCH("")
}

} // namespace DelaunatorUtilities
} // namespace Kratos
