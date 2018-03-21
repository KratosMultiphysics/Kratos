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
    unsigned int max_node_id;
    unsigned int max_elem_id;
    unsigned int max_cond_id;

    // Get the nodes id
    const int nnodes = mrModelPart.Nodes().size();
    for (int i = 0; i < nnodes; i++)
    {
        ModelPart::NodesContainerType::iterator inode = mrModelPart.NodesBegin() + i;
        if (inode->Id() > max_node_id)
            max_node_id = inode->Id();
    }

    // Get the elements id
    const int n_elements = mrModelPart.Elements().size();
    for (int i = 0; i < n_elements; i++)
    {
        ModelPart::ElementsContainerType::iterator ielement = mrModelPart.ElementsBegin() + i;
        elements_id.push_back(ielement->Id());
        if (ielement->Id() > max_elem_id)
            max_elem_id = ielement->Id();
    }

    // Get the conditions id
    const int n_conditions = mrModelPart.Conditions().size();
    for (int i = 0; i < n_conditions; i++)
    {
        ModelPart::ConditionsContainerType::iterator icondition = mrModelPart.ConditionsBegin() + i;
        conditions_id.push_back(icondition->Id());
        if (icondition->Id() > max_cond_id)
            max_elem_id = icondition->Id();
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

    // Get the middle node hash
    std::pair<int, int> node_key;
    node_key = std::minmax(node_a->Id(), node_b->Id());

    // Check if the node exist
    // auto search = mNodesHash.find(node_key);
    // if (search != mNodesHash.end() )
    // {
    //     middle_node = mrModelPart.Nodes()(search->second);
    // }
    // else
    // {
    //     //
    // }
    return middle_node;
}


}  // namespace Kratos.


