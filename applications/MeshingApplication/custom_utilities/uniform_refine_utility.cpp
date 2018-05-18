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
#include "utilities/sub_model_parts_list_utility.h"


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
    const IndexType nnodes = mrModelPart.Nodes().size();
    for (IndexType i = 0; i < nnodes; i++)
    {
        ModelPart::NodesContainerType::iterator inode = mrModelPart.NodesBegin() + i;
        if (inode->Id() > mLastNodeId)
            mLastNodeId = inode->Id();
    }

    // Get the elements id
    const IndexType n_elements = mrModelPart.Elements().size();
    for (IndexType i = 0; i < n_elements; i++)
    {
        ModelPart::ElementsContainerType::iterator ielement = mrModelPart.ElementsBegin() + i;
        if (ielement->Id() > mLastElemId)
            mLastElemId = ielement->Id();
    }

    // Get the conditions id
    const IndexType n_conditions = mrModelPart.Conditions().size();
    for (IndexType i = 0; i < n_conditions; i++)
    {
        ModelPart::ConditionsContainerType::iterator icondition = mrModelPart.ConditionsBegin() + i;
        if (icondition->Id() > mLastCondId)
            mLastCondId = icondition->Id();
    }

    mStepDataSize = mrModelPart.GetNodalSolutionStepDataSize();
    mBufferSize = mrModelPart.GetBufferSize();
    mDofs = mrModelPart.NodesBegin()->GetDofs();

    // Compute the sub model part maps
    std::unordered_map<int, std::vector<std::string>> colors;
    SubModelPartsListUtility colors_utility(mrModelPart);
    colors_utility.ComputeSubModelPartsList(mNodesColorMap, mCondColorMap, mElemColorMap, colors);
    SubModelPartsListUtility::IntersectColors(colors, mIntersections);
    mColorsPointers = SubModelPartsListUtility::GetModelPartColorsPointers(mrModelPart, colors);
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
    const IndexType n_elements = mrModelPart.Elements().size();
    for (IndexType i = 0; i < n_elements; i++)
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
    std::vector<IndexType> elements_id;
    std::vector<IndexType> conditions_id;

    // Get the elements id
    const IndexType n_elements = mrModelPart.Elements().size();
    for (IndexType i = 0; i < n_elements; i++)
    {
        ModelPart::ElementsContainerType::iterator ielement = mrModelPart.ElementsBegin() + i;

        // Check the refinement level of the origin elements
        int step_refine_level = ielement->GetValue(REFINEMENT_LEVEL);
        if (step_refine_level == rThisLevel)
            elements_id.push_back(ielement->Id());
    }

    // Get the conditions id
    const IndexType n_conditions = mrModelPart.Conditions().size();
    for (IndexType i = 0; i < n_conditions; i++)
    {
        ModelPart::ConditionsContainerType::iterator icondition = mrModelPart.ConditionsBegin() + i;

        // Check the refinement level of the origin conditions
        int step_refine_level = icondition->GetValue(REFINEMENT_LEVEL);
        if (step_refine_level == rThisLevel)
            conditions_id.push_back(icondition->Id());
    }

    // Loop the origin elements to create the middle nodes
    for (auto id : elements_id)
    {
        // Get the element
        Element::Pointer p_element = mrModelPart.Elements()(id);

        // Get the refinement level of the origin element
        int step_refine_level = rThisLevel + 1;

        // Get the geometry
        Geometry<NodeType>& geom = p_element->GetGeometry();

        // Loop the edges of the father element and get the nodes
        for (auto edge : geom.Edges())
            CreateNodeInEdge(edge, step_refine_level);

        if (geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4)
            CreateNodeInFace( geom, step_refine_level );
    }

    // TODO: add OMP
    // Create the elements
    for (auto id : elements_id)
    {
        // Get the element
        Element::Pointer p_element = mrModelPart.Elements()(id);

        // Get the refinement level of the origin element
        int step_refine_level = rThisLevel + 1;

        // Get the geometry
        Geometry<NodeType>& geom = p_element->GetGeometry();

        if (geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle2D3)
        {
            // Loop the edges of the father element and get the nodes
            int i_edge = 0;
            std::array<NodeType::Pointer, 3> middle_nodes;
            for (auto edge : geom.Edges())
                middle_nodes[i_edge++] = GetNodeInEdge(edge);

            // Create the sub elements
            std::vector<NodeType::Pointer> sub_element_nodes(3);
            for (int position = 0; position < 4; position++)
            {
                sub_element_nodes = GetSubTriangleNodes(position, geom, middle_nodes);
                CreateElement(p_element, sub_element_nodes, step_refine_level);
            }
        }
        else if (geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4)
        {
            // Loop the edges of the father element and get the nodes
            int i_edge = 0;
            std::array<NodeType::Pointer, 5> middle_nodes;
            for (auto edge : geom.Edges())
                middle_nodes[i_edge++] = GetNodeInEdge(edge);
            middle_nodes[4] = GetNodeInFace( geom );

            // Create the sub elements
            std::vector<NodeType::Pointer> sub_element_nodes(4);
            for (int position = 0; position < 4; position++)
            {
                sub_element_nodes = GetSubQuadrilateralNodes(position, geom, middle_nodes);
                CreateElement(p_element, sub_element_nodes, step_refine_level);
            }
        }
        else
        {
            KRATOS_ERROR << "Your geometry contains " << geom.GetGeometryType() << " which cannot be refined" << std::endl;
        }

        // Once we have created all the sub elements, the origin element must be deleted
        p_element->Set(TO_ERASE, true);
    }

    mrModelPart.RemoveElementsFromAllLevels(TO_ERASE);

    // Loop the origin conditions
    for (auto id : conditions_id)
    {
        // Get the condition
        Condition::Pointer p_condition = mrModelPart.Conditions()(id);

        // Get the refinement level of the origin condition
        int step_refine_level = rThisLevel + 1;

        // Get the geometry
        Geometry<NodeType>& geom = p_condition->GetGeometry();

        if (geom.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Line2D2)
        {
            NodeType::Pointer middle_node = GetNodeInEdge(geom);

            // Create the sub conditions
            std::vector<NodeType::Pointer> sub_condition_nodes(2);
            for (int position = 0; position < 2; position++)
            {
                sub_condition_nodes = GetSubLineNodes(position, geom, middle_node);
                CreateCondition(p_condition, sub_condition_nodes, step_refine_level);
            }
        }
        else
        {
            KRATOS_ERROR << "Your geometry contains " << geom.GetGeometryType() << " which cannot be refined" << std::endl;
        }

        // Once we have created all the sub conditions, the origin conditions must be deleted
        p_condition->Set(TO_ERASE, true);
    }

    mrModelPart.RemoveConditionsFromAllLevels(TO_ERASE);

}


/// Create a middle node on an edge. If the node does not exist, it creates one
template <unsigned int TDim>
void UniformRefineUtility<TDim>::CreateNodeInEdge(
    const EdgeType& rEdge,
    const int& rRefinementLevel
    )
{
    // Get the middle node key
    std::pair<IndexType, IndexType> node_key;
    node_key = std::minmax(rEdge(0)->Id(), rEdge(1)->Id());

    // Check if the node is not yet created
    auto search = mNodesMap.find(node_key);
    if (search == mNodesMap.end() )
    {
        // Create the new node
        const double new_x = 0.5*rEdge(0)->X() + 0.5*rEdge(1)->X();
        const double new_y = 0.5*rEdge(0)->Y() + 0.5*rEdge(1)->Y();
        const double new_z = 0.5*rEdge(0)->Z() + 0.5*rEdge(1)->Z();
        NodeType::Pointer middle_node = mrModelPart.CreateNewNode(++mLastNodeId, new_x, new_y, new_z);

        // Store the node key in the map
        mNodesMap[node_key] = middle_node->Id();

        // interpolate the variables
        CalculateNodalStepData(middle_node, rEdge(0), rEdge(1));

        // Set the refinement level
        int& this_node_level = middle_node->GetValue(REFINEMENT_LEVEL);
        this_node_level = rRefinementLevel;

        // Set the DoF's
        for (typename NodeType::DofsContainerType::const_iterator it_dof = mDofs.begin(); it_dof != mDofs.end(); ++it_dof)
            middle_node->pAddDof(*it_dof);

        // Add the node to the sub model parts
        const int key0 = mNodesColorMap[rEdge(0)->Id()];
        const int key1 = mNodesColorMap[rEdge(1)->Id()];
        const int key = mIntersections[std::minmax(key0, key1)];
        if (key != 0)  // NOTE: key==0 is the main model part
        {
            for (auto sub_model_part : mColorsPointers[key])
            {
                sub_model_part->AddNode(middle_node);
            }
        }
        mNodesColorMap[middle_node->Id()] = key;
    }
}


/// Get the middle node on an edge
template <unsigned int TDim>
typename NodeType::Pointer UniformRefineUtility<TDim>::GetNodeInEdge(const EdgeType& rEdge)
{
    // Initialize the output
    NodeType::Pointer middle_node;

    // Get the middle node key
    std::pair<IndexType, IndexType> node_key;
    node_key = std::minmax(rEdge(0)->Id(), rEdge(1)->Id());

    // Check if the node exist
    auto search = mNodesMap.find(node_key);
    if (search != mNodesMap.end() )
    {
        middle_node = mrModelPart.Nodes()(search->second);
    }
    else
    {
        KRATOS_WARNING("UniformRefineProcess") << "Middle node not found in edge" << rEdge << std::endl;
    }

    return middle_node;
}


/// Get the middle node on a face defined by four nodes. If the node does not exist, it creates one
template< unsigned int TDim>
void UniformRefineUtility<TDim>::CreateNodeInFace(
    const FaceType& rFace,
    const int& rRefinementLevel
    )
{
    // Get the middle node key
    std::array<IndexType, 4> node_key = {rFace(0)->Id(), rFace(1)->Id(), rFace(2)->Id(), rFace(3)->Id()};
    std::sort(node_key.begin(), node_key.end());

    // Check if the node is not yet created
    auto search = mNodesInFaceMap.find(node_key);
    if (search == mNodesInFaceMap.end() )
    {
        // Create the new node
        const double new_x = 0.25*rFace(0)->X() + 0.25*rFace(1)->X() + 0.25*rFace(2)->X() + 0.25*rFace(3)->X();
        const double new_y = 0.25*rFace(0)->Y() + 0.25*rFace(1)->Y() + 0.25*rFace(2)->Y() + 0.25*rFace(3)->Y();
        const double new_z = 0.25*rFace(0)->Z() + 0.25*rFace(1)->Z() + 0.25*rFace(2)->Z() + 0.25*rFace(3)->Z();
        NodeType::Pointer middle_node = mrModelPart.CreateNewNode(++mLastNodeId, new_x, new_y, new_z);

        // Store the node key in the map
        mNodesInFaceMap[node_key] = middle_node->Id();

        // interpolate the variables
        CalculateNodalStepData(middle_node, rFace(0), rFace(1), rFace(2), rFace(3));

        // Set the refinement level
        int& this_node_level = middle_node->GetValue(REFINEMENT_LEVEL);
        this_node_level = rRefinementLevel;

        // Set the DoF's
        for (typename NodeType::DofsContainerType::const_iterator it_dof = mDofs.begin(); it_dof != mDofs.end(); ++it_dof)
            middle_node->pAddDof(*it_dof);

        // Add the node to the sub model parts
        int key0 = mNodesColorMap[rFace(0)->Id()];
        int key1 = mNodesColorMap[rFace(1)->Id()];
        int key2 = mNodesColorMap[rFace(2)->Id()];
        int key3 = mNodesColorMap[rFace(3)->Id()];
        key0 = mIntersections[std::minmax(key0, key1)];
        key1 = mIntersections[std::minmax(key2, key3)];
        key0 = mIntersections[std::minmax(key0, key1)];
        if (key0 != 0)  // NOTE: key==0 is the main model part
        {
            for (auto sub_model_part : mColorsPointers[key0])
            {
                sub_model_part->AddNode(middle_node);
            }
        }
        mNodesColorMap[middle_node->Id()] = key0;
    }
}


/// Get the middle node on a face defined by four nodes. If the node does not exist, it creates one
template< unsigned int TDim>
typename NodeType::Pointer UniformRefineUtility<TDim>::GetNodeInFace(const FaceType& rFace)
{
    // Initialize the output
    NodeType::Pointer middle_node;

    // Get the middle node key
    std::array<IndexType, 4> node_key = {rFace(0)->Id(), rFace(1)->Id(), rFace(2)->Id(), rFace(3)->Id()};
    std::sort(node_key.begin(), node_key.end());

    // Check if the node exist
    auto search = mNodesInFaceMap.find(node_key);
    if (search != mNodesInFaceMap.end() )
    {
        middle_node = mrModelPart.Nodes()(search->second);
    }
    else
    {
        KRATOS_WARNING("UniformRefineProcess") << "Middle node not found in edge" << rFace << std::endl;
    }

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
    for (IndexType step = 0; step < mBufferSize; step++)
    {
        double* new_node_data = pNewNode->SolutionStepData().Data(step);

        const double* node_data_0 = pNode0->SolutionStepData().Data(step);
        const double* node_data_1 = pNode1->SolutionStepData().Data(step);

        for (IndexType variable = 0; variable < mStepDataSize; variable++)
            new_node_data[variable] = 0.5 * node_data_0[variable] + 0.5 * node_data_1[variable];
    }
}


/// Compute the nodal data of a node
template< unsigned int TDim >
void UniformRefineUtility<TDim>::CalculateNodalStepData(
    NodeType::Pointer pNewNode,
    const NodeType::Pointer pNode0,
    const NodeType::Pointer pNode1,
    const NodeType::Pointer pNode2,
    const NodeType::Pointer pNode3
    )
{
    for (IndexType step = 0; step < mBufferSize; step++)
    {
        double* new_node_data = pNewNode->SolutionStepData().Data(step);

        const double* node_data_0 = pNode0->SolutionStepData().Data(step);
        const double* node_data_1 = pNode1->SolutionStepData().Data(step);
        const double* node_data_2 = pNode2->SolutionStepData().Data(step);
        const double* node_data_3 = pNode3->SolutionStepData().Data(step);

        for (IndexType variable = 0; variable < mStepDataSize; variable++)
            new_node_data[variable] = 0.25 * node_data_0[variable] + 0.25 * node_data_1[variable] +
                                      0.25 * node_data_2[variable] + 0.25 * node_data_3[variable];
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
        // Add the element to the origin model part
        mrModelPart.AddElement(sub_element);

        // Set the refinement level
        int& this_elem_level = sub_element->GetValue(REFINEMENT_LEVEL);
        this_elem_level = rRefinementLevel;

        // Add the element to the sub model parts
        int key = mElemColorMap[pOriginElement->Id()];
        if (key != 0)  // NOTE: key==0 is the main model part
        {
            for (auto sub_model_part : mColorsPointers[key])
            {
                sub_model_part->AddElement(sub_element);
            }
        }
        mElemColorMap[sub_element->Id()] = key;
    }
}


/// Create a sub condition
template<unsigned int TDim>
void UniformRefineUtility<TDim>::CreateCondition(
    Condition::Pointer pOriginCondition,
    std::vector<NodeType::Pointer> ThisNodes,
    const int& rRefinementLevel
    )
{
    Condition::Pointer sub_condition = pOriginCondition->Create(++mLastElemId, ThisNodes, pOriginCondition->pGetProperties());

    if (sub_condition != nullptr)
    {
        // Add the element to the origin model part
        mrModelPart.AddCondition(sub_condition);

        // Set the refinement level
        int& this_cond_level = sub_condition->GetValue(REFINEMENT_LEVEL);
        this_cond_level = rRefinementLevel;

        // Add the element to the sub model parts
        int key = mCondColorMap[pOriginCondition->Id()];
        if (key != 0)  // NOTE: key==0 is the main model part
        {
            for (auto sub_model_part : mColorsPointers[key])
            {
                sub_model_part->AddCondition(sub_condition);
            }
        }
        mCondColorMap[sub_condition->Id()] = key;
    }
}


/// Return the nodes defining the i-subline
template<unsigned int TDim>
std::vector<typename NodeType::Pointer> UniformRefineUtility<TDim>::GetSubLineNodes(
    const int Position,
    const Geometry<NodeType>& rGeom,
    NodeType::Pointer& rMiddleNode
    )
{
    std::vector<NodeType::Pointer> sub_line_nodes(2);

    if (Position == 0)
    {
        // First sub line
        sub_line_nodes[0] = rGeom.pGetPoint(0);
        sub_line_nodes[1] = rMiddleNode;
    }
    else if (Position == 1)
    {
        // second sub line
        sub_line_nodes[0] = rMiddleNode;
        sub_line_nodes[1] = rGeom.pGetPoint(1);
    }
    else
    {
        KRATOS_ERROR << "Attempting to get " << Position << " sub-line inside a line" << std::endl;
    }

    return sub_line_nodes;
}


/// Return the nodes defining the i-subtriangle
template<unsigned int TDim>
std::vector<typename NodeType::Pointer> UniformRefineUtility<TDim>::GetSubTriangleNodes(
    const int Position,
    const Geometry<NodeType>& rGeom,
    std::array<NodeType::Pointer, 3>& rMiddleNodes
    )
{
    std::vector<NodeType::Pointer> sub_triangle_nodes(3);

    if (Position == 0)
    {
        // First sub triangle
        sub_triangle_nodes[0] = rGeom.pGetPoint(0);
        sub_triangle_nodes[1] = rMiddleNodes[0];
        sub_triangle_nodes[2] = rMiddleNodes[2];
    }
    else if (Position == 1)
    {
        // Second sub triangle
        sub_triangle_nodes[0] = rGeom.pGetPoint(1);
        sub_triangle_nodes[1] = rMiddleNodes[1];
        sub_triangle_nodes[2] = rMiddleNodes[0];
    }
    else if (Position == 2)
    {
        // Third sub triangle
        sub_triangle_nodes[0] = rGeom.pGetPoint(2);
        sub_triangle_nodes[1] = rMiddleNodes[2];
        sub_triangle_nodes[2] = rMiddleNodes[1];
    }
    else if (Position == 3)
    {
        // Fourth sub triangle (inner triangle)
        sub_triangle_nodes[0] = rMiddleNodes[0];
        sub_triangle_nodes[1] = rMiddleNodes[1];
        sub_triangle_nodes[2] = rMiddleNodes[2];
    }
    else
    {
        KRATOS_ERROR << "Attempting to get " << Position << " sub-triangle inside a triangle" << std::endl;
    }

    return sub_triangle_nodes;
}

/// Return the nodes defining the i-subquadrilateral
template<unsigned int TDim>
std::vector<typename NodeType::Pointer> UniformRefineUtility<TDim>::GetSubQuadrilateralNodes(
    const int Position,
    const Geometry<NodeType>& rGeom,
    std::array<NodeType::Pointer, 5>& rMiddleNodes
    )
{
    std::vector<NodeType::Pointer> sub_triangle_nodes(4);

    if (Position == 0)
    {
        // First sub element
        sub_triangle_nodes[0] = rGeom.pGetPoint(0);
        sub_triangle_nodes[1] = rMiddleNodes[0];
        sub_triangle_nodes[2] = rMiddleNodes[4];
        sub_triangle_nodes[3] = rMiddleNodes[3];
    }
    else if (Position == 1)
    {
        // Second sub element
        sub_triangle_nodes[0] = rGeom.pGetPoint(1);
        sub_triangle_nodes[1] = rMiddleNodes[1];
        sub_triangle_nodes[2] = rMiddleNodes[4];
        sub_triangle_nodes[3] = rMiddleNodes[0];
    }
    else if (Position == 2)
    {
        // Third sub element
        sub_triangle_nodes[0] = rGeom.pGetPoint(2);
        sub_triangle_nodes[1] = rMiddleNodes[2];
        sub_triangle_nodes[2] = rMiddleNodes[4];
        sub_triangle_nodes[3] = rMiddleNodes[1];
    }
    else if (Position == 3)
    {
        // Fourth sub element
        sub_triangle_nodes[0] = rGeom.pGetPoint(3);
        sub_triangle_nodes[1] = rMiddleNodes[3];
        sub_triangle_nodes[2] = rMiddleNodes[4];
        sub_triangle_nodes[3] = rMiddleNodes[2];
    }
    else
    {
        KRATOS_ERROR << "Attempting to get " << Position << " sub-quadrilateral inside a quadrilateral" << std::endl;
    }

    return sub_triangle_nodes;
}



template class UniformRefineUtility<2>;

}  // namespace Kratos.
