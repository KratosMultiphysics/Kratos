//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//


// System includes


// External includes 


// Project includes
#include "includes/define.h"
#include "nested_refinement_utility.h"


namespace Kratos
{
/// Default constructor
template< unsigned int TDim>
NestedRefinementUtility<TDim>::NestedRefinementUtility(ModelPart& rModelPart) :
    mrModelPart(rModelPart) 
{
    // Initialize the member variables storing the Id
    // Get the last node id
    const int nnodes = mrModelPart.Nodes().size();
    for (int i = 0; i < nnodes; i++)
    {
        ModelPart::NodesContainerType::iterator inode = mrModelPart.NodesBegin() + i;
        if (inode->Id() > mLastNodeId)
            mLastNodeId = inode->Id();
    }

    // Get the elements id
    const int n_elements = mrModelPart.Elements().size();
    for (int i = 0; i < n_elements; i++)
    {
        ModelPart::ElementsContainerType::iterator ielement = mrModelPart.ElementsBegin() + i;
        if (ielement->Id() > mLastElemId)
            mLastElemId = ielement->Id();
    }

    // Get the conditions id
    const int n_conditions = mrModelPart.Conditions().size();
    for (int i = 0; i < n_conditions; i++)
    {
        ModelPart::ConditionsContainerType::iterator icondition = mrModelPart.ConditionsBegin() + i;
        if (icondition->Id() > mLastCondId)
            mLastCondId = icondition->Id();
    }
}

/// Destructor
template< unsigned int TDim>
NestedRefinementUtility<TDim>::~NestedRefinementUtility() {}

/// Turn back information as a string.
template< unsigned int TDim>
std::string NestedRefinementUtility<TDim>::Info() const {
    return "Nested refinement utility.";
}

/// Print information about this object.
template< unsigned int TDim>
void NestedRefinementUtility<TDim>::PrintInfo(std::ostream& rOStream) const {
    rOStream << "Nested refinement utility.";
}

/// Print object's data.
template< unsigned int TDim>
void NestedRefinementUtility<TDim>::PrintData(std::ostream& rOStream) const {
    rOStream << "Nested refinement utility constructed with:\n";
    rOStream << "   Model part: " << mrModelPart.Info() << "\n";
}

/// Execute the refinement
template< unsigned int TDim>
void NestedRefinementUtility<TDim>::Refine() 
{
    // Initialize the entities Id lists
    std::vector<int> elements_id;
    std::vector<int> conditions_id;


    // Get the elements id
    const int n_elements = mrModelPart.Elements().size();
    for (int i = 0; i < n_elements; i++)
    {
        ModelPart::ElementsContainerType::iterator ielement = mrModelPart.ElementsBegin() + i;
        elements_id.push_back(ielement->Id());
    }

    // Get the conditions id
    const int n_conditions = mrModelPart.Conditions().size();
    for (int i = 0; i < n_conditions; i++)
    {
        ModelPart::ConditionsContainerType::iterator icondition = mrModelPart.ConditionsBegin() + i;
        conditions_id.push_back(icondition->Id());
    }

    // Loop the origin elements. Get the middle node on each edge and create the nodes
    for (const int id : elements_id)
    {
        Element::Pointer sub_element;
        // Get the element and its geometry
        Element::Pointer p_element = mrModelPart.Elements()(id);
        Geometry<Node<3>>& geom = p_element->GetGeometry();

        // Get the nodes on the edges of the father element
        array_1d<Node<3>::Pointer, 3> p_middle_nodes;
        p_middle_nodes[0] = GetNodeBetween(geom.pGetPoint(0), geom.pGetPoint(1));
        p_middle_nodes[1] = GetNodeBetween(geom.pGetPoint(1), geom.pGetPoint(2));
        p_middle_nodes[2] = GetNodeBetween(geom.pGetPoint(0), geom.pGetPoint(2));
                
        // And now, create the sub triangles
        std::vector<Node<3>::Pointer> sub_element_nodes (3);
        // First sub element
        sub_element_nodes[0] = geom.pGetPoint(0);
        sub_element_nodes[1] = p_middle_nodes[0];
        sub_element_nodes[2] = p_middle_nodes[2];
        sub_element = p_element->Create(mLastElemId++, sub_element_nodes, p_element->pGetProperties());
        

        // Encontrar el lugar para setear NodalDataBase, ElementDataBase and SubModelParts
    }

}

/// Get the middle node on an edge defined by two nodes
template< unsigned int TDim>
Node<3>::Pointer NestedRefinementUtility<TDim>::GetNodeBetween(
    const Node<3>::Pointer pNode0,
    const Node<3>::Pointer pNode1
    )
{
    // Initialize the output
    Node<3>::Pointer middle_node;
    
    // Get the middle node key
    std::pair<int, int> node_key;
    node_key = std::minmax(pNode0->Id(), pNode1->Id());

    // Check if the node exist
    auto search = mNodesMap.find(node_key);
    if (search != mNodesMap.end() )
    {
        middle_node = mrModelPart.Nodes()(search->second);
    }
    else
    {
        // Create the new node
        double new_x = pNode0->X() - pNode1->X();
        double new_y = pNode0->Y() - pNode1->Y();
        double new_z = pNode0->Z() - pNode1->Z();
        middle_node = mrModelPart.CreateNewNode(mLastNodeId++, new_x, new_y, new_z);

        // Store the node in the map
        //std::pair< std::pair<int, int>, int > node_map = (node_key, middle_node->Id());
        mNodesMap.insert( std::pair< std::pair<int, int>, int > (node_key, middle_node->Id()) );
    }

    return middle_node;
}

template class NestedRefinementUtility<2>;

}  // namespace Kratos.


