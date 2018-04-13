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
#include "includes/variables.h"
#include "uniform_refine_utility.h"


namespace Kratos
{
/// Default constructor
template< unsigned int TDim>
UniformRefineUtility<TDim>::UniformRefineUtility(ModelPart& rModelPart, int RefinementLevel) :
    mrModelPart(rModelPart),
    mFinalRefinementLevel(RefinementLevel)
{
    // Initialize the member variables storing the Id
    mLastNodeId = 0;
    mLastElemId = 0;
    mLastCondId = 0;

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

    mStepDataSize = mrModelPart.GetNodalSolutionStepDataSize();
    mBufferSize = mrModelPart.GetBufferSize();
}


/// Destructor
template< unsigned int TDim>
UniformRefineUtility<TDim>::~UniformRefineUtility() {}


/// Turn back information as a string.
template< unsigned int TDim>
std::string UniformRefineUtility<TDim>::Info() const {
    return "Uniform refine utility.";
}


/// Print information about this object.
template< unsigned int TDim>
void UniformRefineUtility<TDim>::PrintInfo(std::ostream& rOStream) const {
    rOStream << "Uniform refine utility.";
}


/// Print object's data.
template< unsigned int TDim>
void UniformRefineUtility<TDim>::PrintData(std::ostream& rOStream) const {
    rOStream << "Uniform refine utility constructed with:\n";
    rOStream << "   Model part: " << mrModelPart.Info() << "\n";
    rOStream << "   Final refinement level: " << mFinalRefinementLevel << "\n";
}


/// Execute the refinement until the final refinement level is reached
template< unsigned int TDim>
void UniformRefineUtility<TDim>::Refine()
{
    // Get the lowest refinement level
    int minimum_refinement_level = 1e6;
    const int n_elements = mrModelPart.Elements().size();
    for (int i = 0; i < n_elements; i++)
    {
        ModelPart::ElementsContainerType::iterator ielement = mrModelPart.ElementsBegin() + i;
        if (ielement->GetValue(REFINEMENT_LEVEL) < minimum_refinement_level)
            minimum_refinement_level = ielement->GetValue(REFINEMENT_LEVEL);
    }

    for (int level = minimum_refinement_level; level < mFinalRefinementLevel; level++)
    {
        RefineLevel(level);
    }
}


/// Execute the refinement once
template< unsigned int TDim>
void UniformRefineUtility<TDim>::RefineLevel(const int& rThisLevel)
{
    // Initialize the entities Id lists
    std::vector<int> elements_id;
    std::vector<int> conditions_id;


    // Get the elements id
    const int n_elements = mrModelPart.Elements().size();
    for (int i = 0; i < n_elements; i++)
    {
        ModelPart::ElementsContainerType::iterator ielement = mrModelPart.ElementsBegin() + i;

        // Check the refinement level of the origin elements
        int step_refine_level = ielement->GetValue(REFINEMENT_LEVEL);
        if (step_refine_level == rThisLevel)
            elements_id.push_back(ielement->Id());
    }

    // Get the conditions id
    const int n_conditions = mrModelPart.Conditions().size();
    for (int i = 0; i < n_conditions; i++)
    {
        ModelPart::ConditionsContainerType::iterator icondition = mrModelPart.ConditionsBegin() + i;

        // Check the refinement level of the origin conditions
        int step_refine_level = icondition->GetValue(REFINEMENT_LEVEL);
        if (step_refine_level == rThisLevel)
            conditions_id.push_back(icondition->Id());
    }

    // Loop the origin elements. Get the middle node on each edge and create the nodes
    for (const int id : elements_id)
    {
        // Get the element
        Element::Pointer p_element = mrModelPart.Elements()(id);

        // Get the refinement level of the origin element
        int step_refine_level = rThisLevel + 1;

        // Get the geometry
        Geometry<NodeType>& geom = p_element->GetGeometry();

        if (geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle2D3)
        {
            // FIRST: Create the nodes
            // Loop the edges of the father element and get the nodes
            std::array<NodeType::Pointer, 3> middle_nodes;
            middle_nodes[0] = GetNodeBetween( geom.pGetPoint(0), geom.pGetPoint(1), step_refine_level );
            middle_nodes[1] = GetNodeBetween( geom.pGetPoint(1), geom.pGetPoint(2), step_refine_level );
            middle_nodes[2] = GetNodeBetween( geom.pGetPoint(0), geom.pGetPoint(2), step_refine_level );

            // SECOND: create the sub elements
            std::vector<NodeType::Pointer> sub_element_nodes(3);

            for (int position = 0; position < 4; position++)
            {
                sub_element_nodes = GetSubTriangleNodes(position, geom, middle_nodes);
                CreateElement(p_element, sub_element_nodes, step_refine_level);
            }
        }
        else if (geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4)
        {
            // FIRST: Create the nodes
            // Loop the edges of the father element and get the nodes
            std::array<NodeType::Pointer, 5> middle_nodes;
            middle_nodes[0] = GetNodeBetween( geom.pGetPoint(0), geom.pGetPoint(1), step_refine_level );
            middle_nodes[1] = GetNodeBetween( geom.pGetPoint(1), geom.pGetPoint(2), step_refine_level );
            middle_nodes[2] = GetNodeBetween( geom.pGetPoint(2), geom.pGetPoint(3), step_refine_level );
            middle_nodes[3] = GetNodeBetween( geom.pGetPoint(0), geom.pGetPoint(3), step_refine_level );
            middle_nodes[4] = GetNodeInFace( geom.pGetPoint(0), geom.pGetPoint(1), geom.pGetPoint(2), geom.pGetPoint(3), step_refine_level );

            // SECOND: create the sub elements
            std::vector<NodeType::Pointer> sub_element_nodes(4);

            for (int position = 0; position < 4; position++)
            {
                sub_element_nodes = GetSubQuadrilateralNodes(position, geom, middle_nodes);
                CreateElement(p_element, sub_element_nodes, step_refine_level);
            }

        }
        else
        {
            KRATOS_WARNING("UniformRefineUtility") << "WARNING: YOUR GEOMETRY CONTAINS " << geom.PointsNumber() <<" NODES CAN NOT BE REMESHED" << std::endl;
        }
        // Encontrar el lugar para ejecutar SubModelPartsColors

        // Once we have created all the sub elements
        p_element->Set(TO_ERASE, true);
    }
    
    mrModelPart.RemoveElementsFromAllLevels(TO_ERASE);

    // Loop the origin conditions
    for (const int id : conditions_id)
    {
        // Get the condition
        Condition::Pointer p_condition = mrModelPart.Conditions()(id);

        // Check the refinement level of the origin condition
        int step_refine_level = rThisLevel + 1;
        // THIRD: Create the conditions

        /* Do some stuff here */

    }

}


/// Get the middle node on an edge defined by two nodes. If the node does not exist, it creates one
template< unsigned int TDim>
Node<3>::Pointer UniformRefineUtility<TDim>::GetNodeBetween(
    const NodeType::Pointer pNode0,
    const NodeType::Pointer pNode1,
    const int& rRefinementLevel
    )
{
    // Initialize the output
    NodeType::Pointer middle_node;
    
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
        double new_x = 0.5*pNode0->X() + 0.5*pNode1->X();
        double new_y = 0.5*pNode0->Y() + 0.5*pNode1->Y();
        double new_z = 0.5*pNode0->Z() + 0.5*pNode1->Z();
        middle_node = mrModelPart.CreateNewNode(++mLastNodeId, new_x, new_y, new_z);

        // interpolate the variables
        CalculateNodalStepData(middle_node, pNode0, pNode1);

        // Set the refinement level
        int& this_node_level = middle_node->GetValue(REFINEMENT_LEVEL);
        this_node_level = rRefinementLevel;

        // Store the node in the map
        //std::pair< std::pair<int, int>, int > node_map = (node_key, middle_node->Id());
        mNodesMap.insert( std::pair< std::pair<int, int>, int > (node_key, middle_node->Id()) );
    }

    return middle_node;
}


/// Get the middle node on a face defined by four nodes. If the node does not exist, it creates one
template< unsigned int TDim>
Node<3>::Pointer UniformRefineUtility<TDim>::GetNodeInFace(
    const NodeType::Pointer pNode0,
    const NodeType::Pointer pNode1,
    const NodeType::Pointer pNode2,
    const NodeType::Pointer pNode3,
    const int& rRefinementLevel
)
{
    // Initialize the output
    NodeType::Pointer middle_node;
    
    // WARNING: I am working on a 2D space, so I am assuming that the node inside the face only belongs to ONE quadrilateral
    // TODO: develop the 3D and check the existance of the node

    // // Get the middle node key
    // std::pair<int, int> node_key;
    // node_key = std::minmax(pNode0->Id(), pNode1->Id());

    // // Check if the node exist
    // auto search = mNodesMap.find(node_key);
    // if (search != mNodesMap.end() )
    // {
    //     middle_node = mrModelPart.Nodes()(search->second);
    // }
    // else
    // {

    // Create the new node
    double new_x = 0.25*pNode0->X() + 0.25*pNode1->X() + 0.25*pNode2->X() + 0.25*pNode3->X();
    double new_y = 0.25*pNode0->Y() + 0.25*pNode1->Y() + 0.25*pNode2->Y() + 0.25*pNode3->Y();
    double new_z = 0.25*pNode0->Z() + 0.25*pNode1->Z() + 0.25*pNode2->Z() + 0.25*pNode3->Z();
    middle_node = mrModelPart.CreateNewNode(++mLastNodeId, new_x, new_y, new_z);

    // interpolate the variables
    CalculateNodalStepData(middle_node, pNode0, pNode1);

    // Set the refinement level
    int& this_node_level = middle_node->GetValue(REFINEMENT_LEVEL);
    this_node_level = rRefinementLevel;

        // // Store the node in the map
        // mNodesMap.insert( std::pair< std::pair<int, int>, int > (node_key, middle_node->Id()) );
    // }

    return middle_node;
}


/// Compute the nodal data of a node
template< unsigned int TDim >
void UniformRefineUtility<TDim>::CalculateNodalStepData(
    NodeType::Pointer pNewNode, 
    const NodeType::Pointer pNode0, 
    const NodeType::Pointer pNode1
    )
{
    for (unsigned int step = 0; step < mBufferSize; step++)
    {
        double* new_node_data = pNewNode->SolutionStepData().Data(step);

        const double* node_data_0 = pNode0->SolutionStepData().Data(step);
        const double* node_data_1 = pNode1->SolutionStepData().Data(step);

        for (unsigned int variable = 0; variable < mStepDataSize; variable++)
            new_node_data[variable] = 0.5 * node_data_0[variable] + 0.5 * node_data_1[variable];
    }
}


/// Create a sub element
template<unsigned int TDim>
void UniformRefineUtility<TDim>::CreateElement(
    Element::Pointer pOriginElement,
    std::vector<NodeType::Pointer> ThisNodes,
    const int& rRefinementLevel
    )
{
    Element::Pointer sub_element = pOriginElement->Create(++mLastElemId, ThisNodes, pOriginElement->pGetProperties());
    
    if (sub_element != nullptr) 
    {
        mrModelPart.AddElement(sub_element);
        int& this_elem_level = sub_element->GetValue(REFINEMENT_LEVEL);
        this_elem_level = rRefinementLevel;
    }

}


/// Return the nodes defining the i-subtriangle
template<unsigned int TDim>
std::vector<Node<3>::Pointer> UniformRefineUtility<TDim>::GetSubTriangleNodes(
    int Position,
    Geometry<NodeType>& rGeom,
    std::array<NodeType::Pointer, 3>& rMiddleNodes
    )
{
    std::vector<NodeType::Pointer> sub_element_nodes(3);

    if (Position == 0)
    {
        // First sub element
        sub_element_nodes[0] = rGeom.pGetPoint(0);
        sub_element_nodes[1] = rMiddleNodes[0];
        sub_element_nodes[2] = rMiddleNodes[2];
    }
    else if (Position == 1)
    {
        // Second sub element
        sub_element_nodes[0] = rGeom.pGetPoint(1);
        sub_element_nodes[1] = rMiddleNodes[1];
        sub_element_nodes[2] = rMiddleNodes[0];
    }
    else if (Position == 2)
    {
        // Third sub element
        sub_element_nodes[0] = rGeom.pGetPoint(2);
        sub_element_nodes[1] = rMiddleNodes[2];
        sub_element_nodes[2] = rMiddleNodes[1];
    }
    else if (Position == 3)
    {
        // Fourth sub element (inner element)
        sub_element_nodes[0] = rMiddleNodes[0];
        sub_element_nodes[1] = rMiddleNodes[1];
        sub_element_nodes[2] = rMiddleNodes[2];
    }
    else
    {
        KRATOS_ERROR << "Attempting to get " << Position << " sub element inside a triangle" << std::endl;
    }

    return sub_element_nodes;
}

/// Return the nodes defining the i-subquadrilateral
template<unsigned int TDim>
std::vector<Node<3>::Pointer> UniformRefineUtility<TDim>::GetSubQuadrilateralNodes(
    int Position,
    Geometry<NodeType>& rGeom,
    std::array<NodeType::Pointer, 5>& rMiddleNodes
    )
{
    std::vector<NodeType::Pointer> sub_element_nodes(4);

    if (Position == 0)
    {
        // First sub element
        sub_element_nodes[0] = rGeom.pGetPoint(0);
        sub_element_nodes[1] = rMiddleNodes[0];
        sub_element_nodes[2] = rMiddleNodes[4];
        sub_element_nodes[3] = rMiddleNodes[3];
    }
    else if (Position == 1)
    {
        // Second sub element
        sub_element_nodes[0] = rGeom.pGetPoint(1);
        sub_element_nodes[1] = rMiddleNodes[1];
        sub_element_nodes[2] = rMiddleNodes[4];
        sub_element_nodes[3] = rMiddleNodes[0];
    }
    else if (Position == 2)
    {
        // Third sub element
        sub_element_nodes[0] = rGeom.pGetPoint(2);
        sub_element_nodes[1] = rMiddleNodes[2];
        sub_element_nodes[2] = rMiddleNodes[4];
        sub_element_nodes[3] = rMiddleNodes[1];
    }
    else if (Position == 3)
    {
        // Fourth sub element
        sub_element_nodes[0] = rGeom.pGetPoint(3);
        sub_element_nodes[1] = rMiddleNodes[3];
        sub_element_nodes[2] = rMiddleNodes[4];
        sub_element_nodes[3] = rMiddleNodes[2];
    }
    else
    {
        KRATOS_ERROR << "Attempting to get " << Position << " sub element inside a quadrilateral" << std::endl;
    }

    return sub_element_nodes;
}



template class UniformRefineUtility<2>;

}  // namespace Kratos.


