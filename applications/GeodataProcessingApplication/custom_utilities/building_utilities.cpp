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

    void BuildingUtilities::CheckOverlap() {
        // const auto& r_nodes_array = mrModelPart.Nodes();
        auto& r_nodes_array = mrModelPart.Nodes();  // probably "const" is not necessary

        // we set all nodes as not visited
        for (int i_node = 0; i_node < static_cast<int>(r_nodes_array.size()); ++i_node) {
            auto p_node = r_nodes_array.begin() + i_node;
            p_node->Set(VISITED, false);
        }

        // we check if a node belongs to at least 2 sub model part
        // (in this case we create a new node and store the old id and the new one in a map)
        std::unordered_map<IndexType, IndexType> new_nodes_map;
        IndexType new_node_id = 0;
        FindMaxNodeId(new_node_id);     // we update the new_node_id

        for (auto& r_sub_model_part : mrModelPart.SubModelParts()) {
            const auto& r_nodes_smp_array = r_sub_model_part.Nodes();        // nodes in sub model part
            for (int i_node = 0; i_node < static_cast<int>(r_nodes_smp_array.size()); ++i_node) {
                auto p_node = r_nodes_smp_array.begin() + i_node;
                if (p_node->Is(VISITED)) {
                    if (!CheckNodeInMap(new_nodes_map, p_node->Id())) {
                        // create new node if the p_node->Id() is already visited and if it is not in the new_nodes_map
                        // we create a new node with different id but same coordinates
                        mrModelPart.CreateNewNode(  new_node_id,
                                                    p_node->Coordinates()[0],
                                                    p_node->Coordinates()[1],
                                                    p_node->Coordinates()[2]);
                        
                        // we insert the new nodes in the map (map{old_node, new_node})
                        new_nodes_map.insert({p_node->Id(), new_node_id});
                        new_node_id++;
                    }
                } else
                    p_node->Set(VISITED, true);
            }
        }

        // TODO: CHECK IF THERE ARE ALL NODES
        // we set all nodes as not visited
        for (int i_node = 0; i_node < static_cast<int>(r_nodes_array.size()); ++i_node) {
            auto p_node = r_nodes_array.begin() + i_node;
            p_node->Set(VISITED, false);
        }

        for (auto& r_sub_model_part : mrModelPart.SubModelParts()) {
            const auto& r_nodes_smp_array = r_sub_model_part.Nodes();
            auto& r_elems_smp_array = r_sub_model_part.Elements();

            for (int i_node = 0; i_node < static_cast<int>(r_nodes_smp_array.size()); ++i_node) {
                auto p_node = r_nodes_smp_array.begin() + i_node;
                if (p_node->Is(VISITED)) {
                    // HERE WE MUST CREATE A NEW ELEMENT WITH THE NEW NODE REFERENCES
                    
                } else
                    p_node->Set(VISITED, true);
            }
        }

    }


    void BuildingUtilities::SplitBuilding() {

    }


    bool BuildingUtilities::CheckNodeInMap(std::unordered_map<IndexType, IndexType> new_nodes_map, IndexType node_id) {
        for (auto node_it = new_nodes_map.find(node_id); node_it != new_nodes_map.end(); ++node_it)
            return true;
    }


    void BuildingUtilities::FindMaxNodeId(IndexType new_node_id) {
        for (auto& r_node : mrModelPart.Nodes()) {
            if (r_node.Id() > new_node_id) {
                new_node_id = r_node.Id();
            }
        }
        new_node_id++;
        KRATOS_INFO("BuildingUtilities") << "the new node ID is: " << new_node_id << std::endl;
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
