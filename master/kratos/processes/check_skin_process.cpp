//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes
#include <unordered_map>

// External includes

// Project includes
#include "processes/check_skin_process.h"
#include "includes/model_part.h"
#include "includes/key_hash.h"

namespace Kratos
{

void CheckSkinProcess::Execute()
{
    KRATOS_TRY;

    KRATOS_ERROR_IF(mrModelPart.Conditions().size() == 0 && mrModelPart.Elements().size() != 0) << "The number of conditions is zero and the number of elements is not, hence the skin can not envelope the domain" << std::endl;

    using hashmap = std::unordered_map<DenseVector<IndexType>, IndexType, KeyHasherRange<DenseVector<IndexType>>, KeyComparorRange<DenseVector<IndexType>>>;
    hashmap edge_map;

    DenseVector<IndexType> ids(2);

    // Add 1 to the counter for every edge find in the model part
    for (auto& r_cond : mrModelPart.Conditions()) {
        const auto edges = r_cond.GetGeometry().GenerateEdges();

        for(IndexType edge=0; edge<edges.size(); edge++) {
            for(IndexType i=0; i<edges[edge].size(); i++) {
                ids[i] = edges[edge][i].Id();
            }

            //*** THE ARRAY OF IDS MUST BE ORDERED!!! ***
            std::sort(ids.begin(), ids.end());

            edge_map[ids] += 1;
        }
    }

    // Now loop over the entire edge map.
    // All values shall have a value of 2
    // If that is not the case throw an error
    for(auto it=edge_map.begin(); it!=edge_map.end(); ++it) {
        if(it->second > 2) {
            KRATOS_ERROR << "ERROR OVERLAPPING CONDITIONS IN SKIN FOUND : " << std::endl << "The edge between nodes " << it->first[0] << " and " << it->first[1] << std::endl << " belongs to an overlapping condition " << std::endl;
        } else if(it->second < 2) {
            KRATOS_ERROR << "ERROR NON CLOSED SKIN " << std::endl << "The edge between nodes " << it->first[0] << " and " << it->first[1] << std::endl << " only appears once, hence it belongs to a non watertight boundary " << std::endl;
        }
    }

    KRATOS_INFO("CheckSkinProcess") << "Checked " << edge_map.size() << " edges in the skin. No gap or overlap found " << std::endl;

    KRATOS_CATCH("");
}

} // namespace Kratos
