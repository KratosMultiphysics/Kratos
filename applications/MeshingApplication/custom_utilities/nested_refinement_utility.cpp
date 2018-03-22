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
NestedRefinementUtility::NestedRefinementUtility(ModelPart& rModelPart) :
    mrModelPart(rModelPart) {}

/// Destructor
NestedRefinementUtility::~NestedRefinementUtility() {}

/// Turn back information as a string.
std::string NestedRefinementUtility::Info() const {
    return "Nested refinement utility.";
}

/// Print information about this object.
void NestedRefinementUtility::PrintInfo(std::ostream& rOStream) const {
    rOStream << "Nested refinement utility.";
}

/// Print object's data.
void NestedRefinementUtility::PrintData(std::ostream& rOStream) const {
    rOStream << "Nested refinement utility constructed with:\n";
    rOStream << "   Model part: " << mrModelPart.Info() << "\n";
}

/// Execute the refinement
void NestedRefinementUtility::Refine() 
{
    // Initialize the entities Id lists
    std::vector<int> elements_id;
    std::vector<int> conditions_id;

    // Get the nodes id
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
        elements_id.push_back(ielement->Id());
        if (ielement->Id() > mLastElemId)
            mLastElemId = ielement->Id();
    }

    // Get the conditions id
    const int n_conditions = mrModelPart.Conditions().size();
    for (int i = 0; i < n_conditions; i++)
    {
        ModelPart::ConditionsContainerType::iterator icondition = mrModelPart.ConditionsBegin() + i;
        conditions_id.push_back(icondition->Id());
        if (icondition->Id() > mLastCondId)
            mLastCondId = icondition->Id();
    }

    // Loop the origin elements. Get the middle node on each edge and create the nodes
    for (const int id : elements_id)
    {
        // Get the element
        Element::Pointer p_element = mrModelPart.Elements()(id);
        // Initialize the vector containing the middle nodes Id
        array_1d<Node<3>::Pointer, 3> p_middle_nodes;

    }

}

/// Get the middle node on an edge
Node<3>::Pointer NestedRefinementUtility::GetNodeBetween(Node<3>::Pointer node_a, Node<3>::Pointer node_b)
{
    // Initialize the output
    Node<3>::Pointer middle_node;
    
    // Get the middle node key
    std::pair<int, int> node_key;
    node_key = std::minmax(node_a->Id(), node_b->Id());

    // Check if the node exist
    auto search = mNodesMap.find(node_key);
    if (search != mNodesMap.end() )
    {
        middle_node = mrModelPart.Nodes()(search->second);
    }
    else
    {
        // Create the new node
        double new_x = node_a->X() - node_b->X();
        double new_y = node_a->Y() - node_b->Y();
        double new_z = node_a->Z() - node_b->Z();
        middle_node = mrModelPart.CreateNewNode(mLastNodeId++, new_x, new_y, new_z);

        // Store the node in the map
        //std::pair< std::pair<int, int>, int > node_map = (node_key, middle_node->Id());
        mNodesMap.insert( std::pair< std::pair<int, int>, int > (node_key, middle_node->Id()) );
    }

    return middle_node;
}


}  // namespace Kratos.


