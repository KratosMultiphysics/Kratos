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
#include "explicit_mesh_moving_utilities.h"

// Project includes
#include "spatial_containers/spatial_containers.h"
#include "utilities/spatial_containers_configure.h"

namespace Kratos
{
    /* Public functions *******************************************************/

void ExplicitMeshMovingUtilities::SearchStructureNodes(
    ModelPart &rModelPart,
    ModelPart &rStructureModelPart,
    const double SearchRadius,
    VectorResultNodesContainerType &rSearchResults,
    VectorDistanceType &rSearchDistanceResults) {

    auto &r_nodes = rModelPart.Nodes();
    auto &r_structure_nodes = rStructureModelPart.Nodes();
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
            
            n_results = bins.SearchObjectsInRadius(r_nodes[i], SearchRadius, it_results, it_results_distances, max_number_of_nodes);

            rSearchResults[i].insert(rSearchResults[i].begin(), local_results.begin(), local_results.begin() + n_results);
            rSearchDistanceResults[i].insert(rSearchDistanceResults[i].begin(), local_results_distances.begin(), local_results_distances.begin() + n_results);
        }
    }
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
