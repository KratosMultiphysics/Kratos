//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Nicola Germano
//
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "utilities/openmp_utils.h"

// Application includes
#include "building_utilities.h"

namespace Kratos
{
    /* Public functions *******************************************************/

    void BuildingUtilities::CheckOverlapElement() {
        // const auto& r_nodes_array = mrModelPart.Nodes();
        auto& r_nodes_array = mrModelPart.Nodes();

        // we set all nodes as not visited
        for (int i_node = 0; i_node < static_cast<int>(r_nodes_array.size()); ++i_node) {
            auto p_node = r_nodes_array.begin() + i_node;
            p_node->Set(VISITED, false);
        }

        /* we check if a node belongs to at least 2 sub model part
         * (if it is, we create a new node and store the old id and the new one in a map) */
        std::unordered_map<IndexType, std::vector<ModelPart::IndexType>> new_nodes_map;      // map{old_node, [new_node_1, ..., new_node_n]}
        IndexType new_node_id = 0;
        FindMaxNodeId(new_node_id);     // we update the new_node_id


        for (auto& r_sub_model_part : mrModelPart.SubModelParts()) {
            const auto& r_nodes_smp_array = r_sub_model_part.Nodes();        // nodes in current sub model part
            for (int i_node = 0; i_node < static_cast<int>(r_nodes_smp_array.size()); ++i_node) {
                auto p_node = r_nodes_smp_array.begin() + i_node;
                if (p_node->Is(VISITED)) {
                    mrModelPart.CreateNewNode(  new_node_id,
                                                p_node->Coordinates()[0],
                                                p_node->Coordinates()[1],
                                                p_node->Coordinates()[2]);
                    
                    std::unordered_map<IndexType, std::vector<ModelPart::IndexType>>::iterator it_map = new_nodes_map.find(p_node->Id());
                    if (it_map != new_nodes_map.end()) {
                        // the value "p_node->Id()" is found on the map
                        it_map->second.push_back(new_node_id);
                    } else {
                        // the value "p_node->Id()" is not found on the map
                        new_nodes_map.insert({p_node->Id(), {new_node_id}});
                    }
                    new_node_id++;

                } else
                    p_node->Set(VISITED, true);
            }
        }
        // *** HERE WE ARE THE new_nodes_map FILLED ***




        // [TODO]: CHECK IF THERE ARE ALL NODES
        // we set all nodes as not visited
        for (int i_node = 0; i_node < static_cast<int>(r_nodes_array.size()); ++i_node) {
            auto p_node = r_nodes_array.begin() + i_node;
            p_node->Set(VISITED, false);
        }



        /*
         * [ITA]
         * loop nei model part
            * loop negli elementi del model part
                * loop nei nodi degli elementi
                * controlliamo se il nodo i-esimo appartiene alla mappa new_nodes_map
                    // se appartiene, modifichiamo il nodo all'inerno dell'elemento prendendo il valore corrispondente al nuovo nodo (dobbiamo quindi ricreare l'elemento)
                    * se appartiene, salvo in una mappa il vecchio nodo e il corrispondente nuovo nodo in modo da poterlo usare successivamente per creare il nuovo elemento
                * prima di uscire dal sub model part i-esimo dobbiamo eliminare il nodo tra i valori della mappa
                  (non possiamo eliminarlo subito perché un nodo può appartenere a più elementi dello stesso sub model part)
                * se il vettore dei valori in new_nodes_map è vuota possiamo eliminare la voce relativa al nodo (per esempio in una situazione del genere {2:[]})
                * creiamo gli elementi con i nuovi nodi appena letti dalla mappa (dobbiamo eliminare anche gli elementi con i vecchi nodi)
            * possiamo quindi passare al model part successivo
         */
        for (auto& r_sub_model_part : mrModelPart.SubModelParts()) {
            const auto& r_nodes_smp_array = r_sub_model_part.Nodes();
            auto& r_elems_smp_array = r_sub_model_part.Elements();

            std::unordered_map<IndexType, std::vector<ModelPart::IndexType>> new_elements_map;    // first: element id; second: node ids vector
            std::unordered_map<IndexType, IndexType> temp_map;  // in this map the nodes that will be replaced will be saved
            std::vector<ModelPart::IndexType> delete_in_map;    // the first node of the vector that points to these values will be deleted [TODO]: improve the name!
            bool create_elems = false;      // it is true if the element already has a visited node (from another sub model part)

            for (int i_elem = 0; i_elem < static_cast<int>(r_elems_smp_array.size()); ++i_elem) {
                auto p_elem = r_elems_smp_array.begin() + i_elem;
                auto& r_geom = p_elem->GetGeometry();     // nodes
                
                for (unsigned int i_node = 0; i_node < r_geom.size(); ++i_node) {
                    auto p_node = r_nodes_smp_array.begin() + i_node;
                    // check if the p_node is in new_nodes_map

                    // std::map<int, std::vector<int>>::iterator it=new_nodes_map.find(1);      // EXAMPLE
                    std::unordered_map<IndexType, std::vector<ModelPart::IndexType>>::iterator it_map = new_nodes_map.find(p_node->Id());
                    if (it_map != new_nodes_map.end()) {
                        // the value "p_node->Id()" is found on the map
                        create_elems = true;                                // at least one element will be recreate

                        r_sub_model_part.AddNodes({static_cast<long unsigned int>(it_map->second.front())});   // we add this new node in the sub model part (is the first id into the vector)
                        delete_in_map.push_back(p_node->Id());              // the first element of the array with key = p_node->Id() will be deleted from the map
                        
                        p_node->SetId(it_map->second.front());              // we update the node id with the id from the map
                        // p_elem->Set(TO_ERASE, true);                        // this element will be deleted and recreate with the new nodes (from the map)
                        p_elem->Set(VISITED, true);                         // this element will be deleted
                    }

                    /*
                     * OLD CODE
                        if (CheckNodeInMap(new_nodes_map, p_node->Id())) {
                            new_nodes_map.second[0]
                            // auto new_node_in_element = new_nodes_map.second()[0]    // the first available node into the value
                            // TODO: we can track this node (because we must delete it after)
                            // create_elems = true;        // now we know that at least one element will be created
                            
                            // temp_map.insert({p_node->Id(), new_nodes_map.second()[0[]]});   // we store the old node id with the first available node into the value in new_nodes_map
                            delete_in_map.push_back(p_node->Id());
                            p_elem->Set(TO_ERASE, true);    // if at least one node is in the map, we must delete and recreate the element
                        }
                    */
                }
            }
            /*
             *** HERE WE MUST DELETE THE NODES (INTO THE MAP) USED RIGHT NOW ***
                * we check if there are some elements into the map with empty vector: in this case we must delete this element from the map
             */
            // now we delete the nodes that are already used from the map
            for (auto& key : delete_in_map) {
                std::vector<ModelPart::IndexType> vect = new_nodes_map[key];    // the vector in the the new_nodes_map with the n-th key
                std::vector<ModelPart::IndexType>::iterator it_vect;
                it_vect = vect.begin();
                vect.erase(it_vect);        // we delete the first value into the vector

                // we check if the vector is empty; if it is, we delete the item from the map
                if (new_nodes_map[key].empty()) {
                    new_nodes_map.erase(key);
                }
            }

            // // we delete all the elements that will be replaced
            // mrModelPart.RemoveElementsFromAllLevels(TO_ERASE);

            // if create_elems is true, we must create new elements
            if (create_elems) {
                KRATOS_INFO("BuildingUtilities") << "Elements will be created!" << std::endl;
                Properties::Pointer p_prop = mrModelPart.pGetProperties(0);

                /*
                 * for each element we will save the information of the nodes
                 * later, we delete the element and recreate it
                 */
                auto& r_elem_array = r_sub_model_part.Elements();
                for( int i_elem = 0; i_elem < static_cast<int>(r_elem_array.size()); ++i_elem) {
                    auto p_elem = r_elem_array.begin() + i_elem;
                    auto& r_geom = p_elem->GetGeometry();     // nodes

                    if (p_elem->Is(VISITED)) {
                        // we recreate only the element with the flag VISITED = true
                        std::vector<ModelPart::IndexType> vect_nodes;       // the nodes we will use to create the new element
                        for (unsigned int i = 0; i < r_geom.size(); ++i) {  // the node ids in n-th element
                            vect_nodes.push_back(r_geom[i]);
                        }
                        IndexType elem_id = p_elem->Id();
                        // we delete the current element...
                        p_elem->Set(TO_ERASE, true);
                        mrModelPart.RemoveElementFromAllLevels(TO_ERASE);
                        // ...and we recreate it
                        Element::Pointer p_new_element = mrModelPart.CreateNewElement(  "Element2D3N",
                                                                                        elem_id,
                                                                                        vect_nodes,
                                                                                        p_prop);
                        // we must add the element into r_sub_model_part
                        r_sub_model_part.AddElement(p_new_element);
                    }
                }
            }
        } // end r_sub_model_part iteration
    }


    /*
     *** [TODO] CHECK IF IT IS NECESSARY ***
    void BuildingUtilities::SplitBuilding() {

    }
    */

    
    /*
     *** [TODO] CHECK IF IT IS NECESSARY ***
    bool BuildingUtilities::CheckNodeInMap(std::unordered_map<IndexType, std::vector<ModelPart::IndexType>> new_nodes_map, IndexType node_id) {
        // for (auto node_it = new_nodes_map.find(node_id); node_it != new_nodes_map.end(); ++node_it)
        //     return true;
        // return false;
        
        if (new_nodes_map.find(node_id) == new_nodes_map.end())
            return false;   // node_id not in new_nodes_map
        else
            return true;    // node_id in new_nodes_map
    }
    */

    /*
     *** [TODO] CHECK IF IT IS NECESSARY ***
    void BuildingUtilities::FillNodesMap(std::unordered_map<IndexType, IndexVector> new_nodes_map, IndexType node_id) {
        if (new_nodes_map.find(node_id) == new_nodes_map.end()) {
            // node_id not in new_nodes_map; so we must create a new component in the map
        } else {
            // node_id in new_nodes_map; so we must append the new node in the map
        }
    }
    */


    void BuildingUtilities::FindMaxNodeId(IndexType new_node_id) {
        for (auto& r_node : mrModelPart.Nodes()) {
            if (r_node.Id() > new_node_id) {
                new_node_id = r_node.Id();
            }
        }
        new_node_id++;
        KRATOS_INFO("BuildingUtilities") << "The new node ID is: " << new_node_id << std::endl;
    }


    /**********************************************************************************************
     * CODE EXAMPLE
    // int findMaxNode(std::vector<element>& elements) {
    //     int max = -1;
    //     for (auto& element : elements)
    //         if (element.nodes[3] > max) max = element.nodes[3];
    //     return max;
    // }

    // void ComputeNormalizedFreeEnergyOnNodesProcess::ObtainMaximumNodeId(std::size_t& rMaxId)
    // {
    //     std::size_t aux = 0;
    //     std::size_t id;

    //     for (auto& r_node : mrModelPart.Nodes()) {
    //         id = r_node.Id();
    //         if (id > aux)
    //             aux = id;
    //     }
    //     rMaxId = aux;
    // }
    **********************************************************************************************/

    /* External functions *****************************************************/

    /// output stream function
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const BuildingUtilities& rThis) {

        rThis.PrintData(rOStream);
        return rOStream;
    }

}
