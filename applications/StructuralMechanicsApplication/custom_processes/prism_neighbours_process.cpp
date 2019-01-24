// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "custom_processes/prism_neighbours_process.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
void PrismNeighboursProcess::Execute()
{
    KRATOS_TRY;

    for(IndexType i = 0; i < mrModelPart.Elements().size(); ++i) {
        auto it_elem = mrModelPart.Elements().begin() + i;
        NodePointerVector& r_neighbour_nodes = it_elem->GetValue(NEIGHBOUR_NODES);
        r_neighbour_nodes.reserve(6); // Just in-plane neighbours
        r_neighbour_nodes.erase(r_neighbour_nodes.begin(),r_neighbour_nodes.end() );

        ElementPointerVector& r_neighbour_elements = it_elem->GetValue(NEIGHBOUR_ELEMENTS);
        r_neighbour_elements.reserve(3); // Just in-plane neighbours
        r_neighbour_elements.erase(r_neighbour_elements.begin(),r_neighbour_elements.end() );
    }

    if (mComputeOnNodes) {
        for(IndexType i = 0; i < mrModelPart.Nodes().size(); ++i) {
            auto it_node = mrModelPart.Nodes().begin() + i;
            NodePointerVector& r_neighbour_nodes = it_node->GetValue(NEIGHBOUR_NODES);
            r_neighbour_nodes.reserve(6); // Just it_node-plane neighbours
            r_neighbour_nodes.erase(r_neighbour_nodes.begin(),r_neighbour_nodes.end() );

            ElementPointerVector& r_neighbour_elements = it_node->GetValue(NEIGHBOUR_ELEMENTS);
            r_neighbour_elements.reserve(3); // Just it_node-plane neighbours
            r_neighbour_elements.erase(r_neighbour_elements.begin(),r_neighbour_elements.end() );
        }
    }

    /* NEIGHBOUR ELEMENTS */
    // Create the HashMapVectorIntElementPointerType
    HashMapVectorIntElementPointerType face_map;

    for(IndexType i = 0; i < mrModelPart.Elements().size(); ++i) {
        auto it_elem = mrModelPart.Elements().begin() + i;

        GeometryType& geom = it_elem->GetGeometry();

        // Insert a pointer to the element identified by the hash value ids if it doesn't exist
        ElementPointerVector& r_neighbour_elements = it_elem->GetValue(NEIGHBOUR_ELEMENTS);

        /* IN-PLANE FACES */
        VectorIndexType ids_4(4);

        /* FACE 1 */
        ids_4[0] = geom[0].Id();
        ids_4[1] = geom[1].Id();
        ids_4[2] = geom[3].Id();
        ids_4[3] = geom[4].Id();

        /*** THE ARRAY OF IDS MUST BE ORDERED!!! ***/
        std::sort(ids_4.begin(), ids_4.end());
        // Check if the elements already exist in the HashMapVectorIntElementPointerType
        HashMapVectorIntElementPointerIteratorType it_check = face_map.find(ids_4);

        if(it_check != face_map.end() ) {
            // If it exists the node is added as a neighbour, reciprocally
            r_neighbour_elements.push_back(it_check->second);
            ElementPointerVector& aux_3 = (it_check->second)->GetValue(NEIGHBOUR_ELEMENTS);
            aux_3.push_back(*it_elem.base());
        } else {
            // If it doesn't exist it is added to the database
            face_map.insert( HashMapVectorIntElementPointerType::value_type(ids_4, *it_elem.base()) );
        }

        /* FACE 2 */
        ids_4[0] = geom[1].Id();
        ids_4[1] = geom[2].Id();
        ids_4[2] = geom[4].Id();
        ids_4[3] = geom[5].Id();

        /*** THE ARRAY OF IDS MUST BE ORDERED!!! ***/
        std::sort(ids_4.begin(), ids_4.end());
        // Check if the elements already exist in the HashMapVectorIntElementPointerType
        it_check = face_map.find(ids_4);

        if(it_check != face_map.end() ) {
            // If it exists the node is added as a neighbour, reciprocally
            r_neighbour_elements.push_back(it_check->second);
            ElementPointerVector& aux_el_2 = (it_check->second)->GetValue(NEIGHBOUR_ELEMENTS);
            aux_el_2.push_back(*it_elem.base());
        } else {
            // If it doesn't exist it is added to the database
            face_map.insert( HashMapVectorIntElementPointerType::value_type(ids_4, *it_elem.base()) );
        }

        /* FACE 3 */
        ids_4[0] = geom[0].Id();
        ids_4[1] = geom[2].Id();
        ids_4[2] = geom[3].Id();
        ids_4[3] = geom[5].Id();

        /*** THE ARRAY OF IDS MUST BE ORDERED!!! ***/
        std::sort(ids_4.begin(), ids_4.end());
        // Check if the elements already exist in the HashMapVectorIntElementPointerType
        it_check = face_map.find(ids_4);

        if(it_check != face_map.end() ) {
            // If it exists the node is added as a neighbour, reciprocally
            r_neighbour_elements.push_back(it_check->second);
            ElementPointerVector& aux_el_3 = (it_check->second)->GetValue(NEIGHBOUR_ELEMENTS);
            aux_el_3.push_back(*it_elem.base());
        } else {
            // If it doesn't exist it is added to the database
            face_map.insert( HashMapVectorIntElementPointerType::value_type(ids_4, *it_elem.base()) );
        }
    }

    /* NEIGHBOURS NODES */

    // Create the ids and aux vectors
    VectorIndexType ids(2),  aux_1(2), aux_2(2), aux_3(2);

    // Search the neighbour nodes(for elements)
    for(IndexType i = 0; i < mrModelPart.Elements().size(); ++i) {
        auto it_elem = mrModelPart.Elements().begin() + i;

        GeometryType& geom = it_elem->GetGeometry();
        NodePointerVector& neighb_nodes = it_elem->GetValue(NEIGHBOUR_NODES);
        neighb_nodes.resize(6);

        for (IndexType fill = 0; fill < 6; ++fill) {
            neighb_nodes(fill) = NodeType::WeakPointer(geom(fill));
        }

        // Just upper nodes, the others are +3 IDs
        ElementPointerVector& r_neighbour_elements = it_elem->GetValue(NEIGHBOUR_ELEMENTS);
        for (IndexType j = 0; j < r_neighbour_elements.size(); ++j) {

            auto it_neigh_elem = r_neighbour_elements.begin() + j;
            HashMapVectorIntIntType node_map;
            GeometryType& geom_neig = it_neigh_elem->GetGeometry();

            // Edge 1
            ids[0] = geom_neig[0].Id();
            ids[1] = geom_neig[1].Id();
            std::sort(ids.begin(), ids.end());
            node_map.insert( HashMapVectorIntIntType::value_type(ids, 2) );

            // Edge 2
            ids[0] = geom_neig[1].Id();
            ids[1] = geom_neig[2].Id();
            std::sort(ids.begin(), ids.end());
            node_map.insert( HashMapVectorIntIntType::value_type(ids, 0) );

            // Edge 3
            ids[0] = geom_neig[2].Id();
            ids[1] = geom_neig[0].Id();
            std::sort(ids.begin(), ids.end());
            node_map.insert( HashMapVectorIntIntType::value_type(ids, 1) );

            aux_1[0] = geom[1].Id();
            aux_1[1] = geom[2].Id();

            aux_2[0] = geom[2].Id();
            aux_2[1] = geom[0].Id();

            aux_3[0] = geom[0].Id();
            aux_3[1] = geom[1].Id();

            std::sort(aux_1.begin(), aux_1.end());
            HashMapVectorIntIntIteratorType it_1 = node_map.find(aux_1);
            std::sort(aux_2.begin(), aux_2.end());
            HashMapVectorIntIntIteratorType it_2 = node_map.find(aux_2);
            std::sort(aux_3.begin(), aux_3.end());
            HashMapVectorIntIntIteratorType it_3 = node_map.find(aux_3);

            if(it_1 != node_map.end() ) {
                neighb_nodes(0) = NodeType::WeakPointer(geom_neig(it_1->second));
                neighb_nodes(3) = NodeType::WeakPointer(geom_neig(it_1->second + 3));
            } else if(it_2 != node_map.end() ) {
                neighb_nodes(1) = NodeType::WeakPointer(geom_neig(it_2->second));
                neighb_nodes(4) = NodeType::WeakPointer(geom_neig(it_2->second + 3));
            } else if(it_3 != node_map.end() ) {
                neighb_nodes(2) = NodeType::WeakPointer(geom_neig(it_3->second));
                neighb_nodes(5) = NodeType::WeakPointer(geom_neig(it_3->second + 3));
            }
        }
    }

    // We add the neighbour to the nodes
    if (mComputeOnNodes) {
        // Add the neighbour elements to all the nodes in the mesh
        for(IndexType i = 0; i < mrModelPart.Elements().size(); ++i) {
            auto it_elem = mrModelPart.Elements().begin() + i;

            GeometryType& geom = it_elem->GetGeometry();
            for(IndexType j = 0; j < geom.size(); ++j) {
                (geom[j].GetValue(NEIGHBOUR_ELEMENTS)).push_back( Element::WeakPointer( *(it_elem.base()) ) );
            }
        }

        // Adding the neighbouring nodes (in the same face)
        for(IndexType i = 0; i < mrModelPart.Elements().size(); ++i) {
            auto it_elem = mrModelPart.Elements().begin() + i;

            GeometryType& geom = it_elem->GetGeometry();

            NodePointerVector& neighb_nodes = it_elem->GetValue(NEIGHBOUR_NODES);

            for (IndexType j = 0; j < 3; ++j) {
                NodeType::WeakPointer temp;

                // Adding nodes from the element
                IndexType aux_index1, aux_index2;

                if (j == 0) {
                    aux_index1 = 2;
                    aux_index2 = 1;
                } else if (j == 1) {
                    aux_index1 = 0;
                    aux_index2 = 2;
                } else {
                    aux_index1 = 1;
                    aux_index2 = 0;
                }

                // Lower face
                temp = geom(aux_index1);
                AddUniqueWeakPointer<NodeType>(geom[j].GetValue(NEIGHBOUR_NODES), temp);
                temp = geom(aux_index2);
                AddUniqueWeakPointer<NodeType>(geom[j].GetValue(NEIGHBOUR_NODES), temp);

                // Upper face
                temp = geom(aux_index1 + 3);
                AddUniqueWeakPointer<NodeType>(geom[j + 3].GetValue(NEIGHBOUR_NODES), temp);
                temp = geom(aux_index2 + 3);
                AddUniqueWeakPointer<NodeType>(geom[j + 3].GetValue(NEIGHBOUR_NODES), temp);

                // Adding neighbour elements
                if (neighb_nodes[aux_index1].Id() != geom[j].Id()) {
                    // Lower face
                    temp = neighb_nodes(aux_index1);
                    AddUniqueWeakPointer<NodeType>(geom[j].GetValue(NEIGHBOUR_NODES), temp);
                    // Upper face
                    temp = neighb_nodes(aux_index1 + 3);
                    AddUniqueWeakPointer<NodeType>(geom[j + 3].GetValue(NEIGHBOUR_NODES), temp);
                }
                if (neighb_nodes[aux_index2].Id() != geom[j].Id()) {
                    // Lower face
                    temp = neighb_nodes(aux_index2);
                    AddUniqueWeakPointer<NodeType>(geom[j].GetValue(NEIGHBOUR_NODES), temp);
                    // Upper face
                    temp = neighb_nodes(aux_index2 + 3);
                    AddUniqueWeakPointer<NodeType>(geom[j + 3].GetValue(NEIGHBOUR_NODES), temp);
                }
            }
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void PrismNeighboursProcess::ExecuteFinalize()
{
    ClearNeighbours();
}

/***********************************************************************************/
/***********************************************************************************/

void PrismNeighboursProcess::ClearNeighbours()
{
    for(IndexType i = 0; i < mrModelPart.Elements().size(); ++i) {
        auto it_elem = mrModelPart.Elements().begin() + i;
        ElementPointerVector& r_neighbour_elements = it_elem->GetValue(NEIGHBOUR_ELEMENTS);
        r_neighbour_elements.erase(r_neighbour_elements.begin(),r_neighbour_elements.end());
    }

    if (mComputeOnNodes) {
        for(IndexType i = 0; i < mrModelPart.Nodes().size(); ++i) {
            auto it_node = mrModelPart.Nodes().begin() + i;
            NodePointerVector& r_neighbour_nodes = it_node->GetValue(NEIGHBOUR_NODES);
            r_neighbour_nodes.erase(r_neighbour_nodes.begin(),r_neighbour_nodes.end() );

            ElementPointerVector& r_neighbour_elements = it_node->GetValue(NEIGHBOUR_ELEMENTS);
            r_neighbour_elements.erase(r_neighbour_elements.begin(),r_neighbour_elements.end());
        }
    }
}

// class PrismNeighboursProcess

} // namespace Kratos
