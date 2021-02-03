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

    /*** THIS FUNCTION WILL BE SPLITTED ***/
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

        /*** new_nodes_map creation ***/
        for (auto& r_sub_model_part : mrModelPart.SubModelParts()) {
            // KRATOS_INFO("\n\n\nBuildingUtilities") << "r_sub_model_part " << r_sub_model_part.Name() << std::endl;
            const auto& r_nodes_smp_array = r_sub_model_part.Nodes();        // nodes in current sub model part
            for (int i_node = 0; i_node < static_cast<int>(r_nodes_smp_array.size()); ++i_node) {
                auto p_node = r_nodes_smp_array.begin() + i_node;
                if (p_node->Is(VISITED)) {
                    mrModelPart.CreateNewNode(  new_node_id,
                                                p_node->Coordinates()[0],
                                                p_node->Coordinates()[1],
                                                p_node->Coordinates()[2]);
                    
                    // we check if "p_node->Id()" is already in the map
                    std::unordered_map<IndexType, std::vector<ModelPart::IndexType>>::iterator it_map = new_nodes_map.find(p_node->Id());
                    if (it_map != new_nodes_map.end()) {
                        // the value "p_node->Id()" is found on the map. We add the "new_node_id" in the vector
                        it_map->second.push_back(new_node_id);
                    } else {
                        // the value "p_node->Id()" is not found on the map. We create a new item in the map
                        new_nodes_map.insert({p_node->Id(), {new_node_id}});
                    }
                    new_node_id++;

                } else
                    p_node->Set(VISITED, true);
            }
        }   // the "new_nodes_map" is filled

        // // we restore all flag of the nodes as not visited
        // std::cout << "\nmrModelPart.Nodes():" << std::endl;
        // r_nodes_array = mrModelPart.Nodes();
        // for (int i_node = 0; i_node < static_cast<int>(r_nodes_array.size()); ++i_node) {
        //     auto p_node = r_nodes_array.begin() + i_node;
        //     p_node->Set(VISITED, false);
        // }



        for (auto& r_sub_model_part : mrModelPart.SubModelParts()) {
            auto& r_elems_smp_array = r_sub_model_part.Elements();

            std::unordered_map<IndexType, std::vector<ModelPart::IndexType>> new_elements_map;    // first: element id; second: node ids vector
            std::vector<ModelPart::IndexType> delete_in_map;

            // std::cout << "\n\n*** r_sub_model_part " << r_sub_model_part.Name() << std::endl;      // [DEBUG]

            for (int i_elem = 0; i_elem < static_cast<int>(r_elems_smp_array.size()); ++i_elem) {
                auto p_elem = r_elems_smp_array.begin() + i_elem;
                auto& r_geom = p_elem->GetGeometry();     // nodes
                std::vector<ModelPart::IndexType> elem_nodes;       // nodes of the n-th element

                bool erase_element = false;     // we set this variable as true if the current element will be deleted and recreated
                
                // std::cout << "---> p_elem->Id() " << p_elem->Id() << std::endl;     // [DEBUG]
                
                for (unsigned int i_node = 0; i_node < r_geom.size(); ++i_node) {
                    // std::cout << "\t---> r_geom[" << i_node << "]Id() " << r_geom[i_node].Id() << std::endl;        
                    std::unordered_map<IndexType, std::vector<ModelPart::IndexType>>::iterator it_map = new_nodes_map.find(r_geom[i_node].Id());
                    if (it_map != new_nodes_map.end()) {
                        // the value "r_geom[i_node].Id()" is found on the map
                        delete_in_map.push_back(r_geom[i_node].Id());     // we add "r_geom[i_node].Id()" in the vector because we will delete it from the map

                        erase_element = true;
                        elem_nodes.push_back(it_map->second.front());     // we get the node from the map

                        p_elem->Set(TO_ERASE, true);          // this element will be deleted and recreated
                        r_geom[i_node].Set(TO_ERASE, true);   // the node will be deleted from current sub model part
                    } else
                        elem_nodes.push_back(r_geom[i_node].Id());
                }
                if (erase_element) {
                    // we fill the "new_elements_map" with all the information that we will need to recreate the element
                    new_elements_map.insert({p_elem->Id(), elem_nodes});
                }
                // std::cout << "\n";      // [DEBUG]
            }

            // // [DEBUG]
            // std::cout << "\nDELETE_IN_MAP: " << std::endl;
            // for (auto &A : delete_in_map) {
            //     std::cout << A << std::endl;
            // }

            // std::cout << "\n\n";    // [DEBUG]
            // WE MUST DELETE THE NODE USED RIGHT NOW FROM THE MAP
            if (!delete_in_map.empty()) {
                // we delete the duplicate in "delete_in_map"
                std::set<ModelPart::IndexType> set_delete_node(delete_in_map.begin(), delete_in_map.end());
                delete_in_map.assign(set_delete_node.begin(), set_delete_node.end());

                // now we delete the nodes that are already used from the map
                for (auto& key : delete_in_map) {
                    // std::cout << "The node " << key << " will be deleted from the map!" << std::endl;       // [DEBUG]
                    std::vector<ModelPart::IndexType>::iterator it_key = new_nodes_map[key].begin();
                    new_nodes_map[key].erase(it_key);        // we delete the first value into the vector

                    // we check if the vector is empty; if it is, we delete the item from the map
                    // std::cout << "new_nodes_map[" << key << "] " << new_nodes_map[key] << std::endl;        // [DEBUG]
                    if (new_nodes_map[key].empty()) {
                        new_nodes_map.erase(key);
                        
                        // restore the flag "VISITED" for this node
                        mrModelPart.GetNode(key).Set(VISITED, false);
                        
                        
                        // // [ONLY FOR DEBUG PURPOSE]
                        // if (!mrModelPart.GetNode(key).Is(VISITED))
                        //     std::cout << "\t*** The node " << mrModelPart.GetNode(key).Id() << " now is VISITED=false" << std::endl;
                    }
                }
            }
            // int num_nodes_before = r_sub_model_part.NumberOfNodes();    // [DEBUG]

            // we delete all the elements that will be replaced
            mrModelPart.RemoveElementsFromAllLevels(TO_ERASE);
            r_sub_model_part.RemoveNodes(TO_ERASE);     // we delete the nodes only from the current sub model part
            
            // int num_nodes_after = r_sub_model_part.NumberOfNodes();     // [DEBUG]
            // std::cout << "\n\n***** In total " << (num_nodes_before-num_nodes_after) << " nodes are deleted! *****\n" << std::endl;     // [DEBUG]
            // std::cout << "Now there are " << num_nodes_after << " in the sub model part " << r_sub_model_part.Name() << std::endl;      // [DEBUG]

            // if "new_elements_map" is empty, we skip these steps; otherwise we create the element with the information in the "new_elements_map"
            if (!new_elements_map.empty()) {
                Properties::Pointer p_prop = mrModelPart.pGetProperties(0);

                // loop over the "new_elements_map"
                for (auto new_elem = new_elements_map.begin(); new_elem != new_elements_map.end(); ++new_elem) { 
                    r_sub_model_part.AddNodes(new_elem->second);
                    Element::Pointer p_new_element = mrModelPart.CreateNewElement(  "Element2D3N",
                                                                                    new_elem->first,
                                                                                    new_elem->second,
                                                                                    p_prop);
                    r_sub_model_part.AddElement(p_new_element);
                    // std::cout << "The element " << new_elem->first << " is added in the sub model part " << r_sub_model_part.Name() << std::endl;      // [DEBUG]
                    // std::cout << "The node " << new_elem->second << " is added in the sub model part " << r_sub_model_part.Name() << std::endl;      // [DEBUG]
                }
                auto& r_nodes_array = mrModelPart.Nodes();
                for (int i_node = 0; i_node < static_cast<int>(r_nodes_array.size()); ++i_node) {
                    auto p_node = r_nodes_array.begin() + i_node;
                    p_node->Set(TO_ERASE, false);
                }
            }
        }
    }


    void BuildingUtilities::DeleteNotValidElements() {
        auto& r_elements_array = mrModelPart.Elements();

        // we set all elements as TO_ERASE=false
        for (int i_elem = 0; i_elem < static_cast<int>(r_elements_array.size()); ++i_elem) {
            auto p_elem = r_elements_array.begin() + i_elem;
            p_elem->Set(TO_ERASE, false);
        }

        // we delete all elements with at last two equal nodes
        for (int i_elem = 0; i_elem < static_cast<int>(r_elements_array.size()); ++i_elem) {
            auto p_elem = r_elements_array.begin() + i_elem;
            auto& r_geom = p_elem->GetGeometry();     // nodes

            std::vector<ModelPart::IndexType> elem_nodes;

            for (unsigned int i_node = 0; i_node < r_geom.size(); ++i_node) {
                std::vector<ModelPart::IndexType>::iterator it_vect = std::find (elem_nodes.begin(), elem_nodes.end(), r_geom[i_node].Id()); 
                if (it_vect != elem_nodes.end()) {
                    p_elem->Set(TO_ERASE, true);
                } else
                    elem_nodes.push_back(r_geom[i_node].Id());
            }
        }
        mrModelPart.RemoveElementsFromAllLevels(TO_ERASE);
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


    void BuildingUtilities::FindMaxNodeId(IndexType &new_node_id) {
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
