//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Simon Wenczowski
//
//

// System includes

// External includes

// Project includes
#include "utilities/openmp_utils.h"

// Application includes
#include "cleaning_utilities.h"

namespace Kratos
{
    /* Public functions *******************************************************/

    void CleaningUtilities::CleanIsolatedNodes(){

        const int initial_num = mrModelPart.Nodes().size();
        auto& r_nodes_array = mrModelPart.Nodes();
        const auto& r_elem_array = mrModelPart.Elements();

        // marking all nodes as "superfluous"
        #pragma omp parallel for
        for( int i_node = 0; i_node < static_cast<int>(r_nodes_array.size()); ++i_node ){
            auto node = r_nodes_array.begin() + i_node;
            node->Set(TO_ERASE, true);
        }

        // saving the nodes that belong to an element
        #pragma omp parallel for
        for( int i_elem = 0; i_elem < static_cast<int>(r_elem_array.size()); ++i_elem ){
            const auto elem = r_elem_array.begin() + i_elem;
            auto& r_geom = elem->GetGeometry();

            for (unsigned int i = 0; i < r_geom.size(); ++i){
                r_geom[i].Set(TO_ERASE, false);
            }
        }

        mrModelPart.RemoveNodesFromAllLevels(TO_ERASE);
        const int final_num = mrModelPart.Nodes().size();
        KRATOS_INFO("CleaningUtilities") << "In total " << (initial_num - final_num) <<" superfluous nodes were cleared" << std::endl;

    }


    void CleaningUtilities::ReBuildModelPart(){

        KRATOS_INFO("CleaningUtilities") << "waiting for implementation" << std::endl;

    }

    /* External functions *****************************************************/

    /// output stream function
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const CleaningUtilities& rThis) {

        rThis.PrintData(rOStream);
        return rOStream;
    }

}
