//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes


// External includes

// Project includes
#include "processes/find_nodal_neighbours_process.h"

namespace Kratos
{

void FindNodalNeighboursProcess::Execute()
{
    // The array of nodes
    const NodesContainerType& r_nodes_array = mrModelPart.Nodes();

    // First of all the neighbour nodes and elements array are initialized to the guessed size
    // and empties the old entries
    for(auto& r_node : r_nodes_array) {
        (r_node.GetValue(NEIGHBOUR_NODES)).reserve(mAverageNodes);
        auto& r_neighbour_nodes = r_node.GetValue(NEIGHBOUR_NODES);
        r_neighbour_nodes.erase(r_neighbour_nodes.begin(),r_neighbour_nodes.end());

        (r_node.GetValue(NEIGHBOUR_ELEMENTS)).reserve(mAverageElements);
        auto& r_neighbour_elements = r_node.GetValue(NEIGHBOUR_ELEMENTS);
        r_neighbour_elements.erase(r_neighbour_elements.begin(),r_neighbour_elements.end());
    }

    // Add the neighbour elements to all the nodes in the mesh
    ElementsContainerType& r_elements_array = mrModelPart.Elements();
    for(auto& r_elem : r_elements_array) {
        GeometryType& r_geometry = r_elem.GetGeometry();
        for(IndexType i = 0; i < r_geometry.size(); i++) {
            (r_geometry[i].GetValue(NEIGHBOUR_ELEMENTS)).push_back(Element::WeakPointer(&r_elem));
        }
    }

    // Adding the neighbouring nodes
    for(auto& r_node : r_nodes_array) {
        auto& r_neighbour_elements = r_node.GetValue(NEIGHBOUR_ELEMENTS);

        for(IndexType i_elem = 0; i_elem < r_neighbour_elements.size(); i_elem++) {
            GeometryType& r_geometry = r_neighbour_elements[i_elem].GetGeometry();
            for(IndexType i = 0; i < r_geometry.size(); ++i) {
                if(r_geometry[i].Id() != r_node.Id() ) {
                    GlobalPointer<NodeType> temp(r_geometry(i));
                    AddUniqueWeakPointer<NodeType>(r_node.GetValue(NEIGHBOUR_NODES), temp);
                }
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void FindNodalNeighboursProcess::ClearNeighbours()
{
    NodesContainerType& r_nodes_array = mrModelPart.Nodes();
    for(auto& r_node : r_nodes_array) {
        auto& r_neighbour_elements = r_node.GetValue(NEIGHBOUR_ELEMENTS);
        r_neighbour_elements.erase(r_neighbour_elements.begin(),r_neighbour_elements.end());

        auto& r_neighbour_nodes = r_node.GetValue(NEIGHBOUR_NODES);
        r_neighbour_nodes.erase(r_neighbour_nodes.begin(),r_neighbour_nodes.end() );
    }
}
  
}  // namespace Kratos.


