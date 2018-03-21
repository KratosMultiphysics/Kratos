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

///
void NestedRefinementUtility::Refine() 
{
    // Initialize the nodes hash
    std::unordered_map<array_1d<int,2>, int, KeyHasherRange<array_1d<int,2>>, KeyComparorRange<array_1d<int,2>> > nodes_hash;

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
    const int nelements = mrModelPart.Elements().size();
    for (int i = 0; i < nelements; i++)
    {
        ModelPart::ElementsContainerType::iterator ielement = mrModelPart.ElementsBegin() + i;
        elements_id.push_back(ielement->Id());
        if (ielement->Id() > max_elem_id)
            max_elem_id = ielement->Id();
    }

    // Get the conditions id
    const int nconditions = mrModelPart.Conditions().size();
    for (int i = 0; i < nconditions; i++)
    {
        ModelPart::ConditionsContainerType::iterator icondition = mrModelPart.ConditionsBegin() + i;
        conditions_id.push_back(icondition->Id());
        if (icondition->Id() > max_cond_id)
            max_elem_id = icondition->Id();
    }
}

}  // namespace Kratos.


