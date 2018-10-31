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
#include "uniform_refinement_utility.h"
#include "utilities/assign_unique_model_part_collection_tag_utility.h"


namespace Kratos
{
/// Default constructor
UniformRefinementUtility::UniformRefinementUtility(ModelPart& rModelPart) :
    mrModelPart(rModelPart)
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
    mDimension = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];
}


/// Destructor
UniformRefinementUtility::~UniformRefinementUtility() {}


/// Turn back information as a string.
std::string UniformRefinementUtility::Info() const {
    return "Uniform refine utility.";
}


/// Print information about this object.
void UniformRefinementUtility::PrintInfo(std::ostream& rOStream) const {
    rOStream << "Uniform refine utility.";
}


/// Print object's data.
void UniformRefinementUtility::PrintData(std::ostream& rOStream) const {
    rOStream << "Uniform refine utility constructed with:\n";
    rOStream << "   Model part: " << mrModelPart.Info() << "\n";
}


/// Execute the refinement until the final number of divisions level is reached
void UniformRefinementUtility::Refine(int& rFinalRefinementLevel)
{
    if (mrModelPart.Nodes().size() == 0)
        KRATOS_WARNING("UniformRefinementUtility") << "Attempting to refine an empty model part" << std::endl;
    else
        mDofs = mrModelPart.NodesBegin()->GetDofs();

    // Get the lowest refinement level
    int minimum_divisions_level = 1e6;
    const IndexType n_elements = mrModelPart.Elements().size();
    for (IndexType i = 0; i < n_elements; i++)
    {
        ModelPart::ElementsContainerType::iterator ielement = mrModelPart.ElementsBegin() + i;
        if (ielement->GetValue(NUMBER_OF_DIVISIONS) < minimum_divisions_level)
            minimum_divisions_level = ielement->GetValue(NUMBER_OF_DIVISIONS);
    }

    // Restart the model part collections utility
    mNodesTags.clear();
    mElementsTags.clear();
    mConditionsTags.clear();
    AssignUniqueModelPartCollectionTagUtility collections_utility(mrModelPart);
    collections_utility.ComputeTags(mNodesTags, mConditionsTags, mElementsTags, mCollections);

    IndexIndexVectorMapType tag_nodes, tag_elements, tag_conditions;

    for (int divisions = minimum_divisions_level; divisions < rFinalRefinementLevel; divisions++)
        ExecuteDivision(divisions, tag_nodes, tag_elements, tag_conditions);

    // Finally, add the new entities to the sub model parts
    for (auto& collection : mCollections)
    {
        const IndexType tag = collection.first;
        if (tag != 0) // NOTE: tag == 0 is the root model part
        {
            for (auto model_part_name : collection.second)
            {
                ModelPart& sub_model_part = mrModelPart.GetSubModelPart(model_part_name);
                sub_model_part.AddNodes(tag_nodes[tag]);
                sub_model_part.AddElements(tag_elements[tag]);
                sub_model_part.AddConditions(tag_conditions[tag]);
            }
        }
    }
}


/// Set the custom ids which will be used to create new entities
void UniformRefinementUtility::SetCustomIds(IndexType& rNodeId, IndexType& rElemId, IndexType& rCondId)
{
    // Set the id
    mLastNodeId = rNodeId;
    mLastElemId = rElemId;
    mLastCondId = rCondId;
}


/// Get the last id of the created nodes, elements and conditions
void UniformRefinementUtility::GetLastCreatedIds(IndexType& rNodeId, IndexType& rElemId, IndexType& rCondId)
{
    // Get the id
    rNodeId = mLastNodeId;
    rElemId = mLastElemId;
    rCondId = mLastCondId;
}


/// Remove the refined entities
void UniformRefinementUtility::RemoveRefinedEntities(Flags ThisFlag)
{
    // Clear the maps
    for (ModelPart::NodeIterator node = mrModelPart.NodesBegin(); node < mrModelPart.NodesEnd(); node++)
    {
        if (node->Is(ThisFlag))
        {
            for (NodesInEdgeMapType::iterator pair = mNodesMap.begin(); pair != mNodesMap.end(); )
            {
                if (node->Id() == pair->second)
                    pair = mNodesMap.erase(pair);
                else
                    pair++;
            }

            for (NodesInFaceMapType::iterator pair = mNodesInFaceMap.begin(); pair != mNodesInFaceMap.end(); )
            {
                if (node->Id() == pair->second)
                    pair = mNodesInFaceMap.erase(pair);
                else
                    pair++;
            }
        }
    }

    // Remove the entities
    mrModelPart.RemoveNodesFromAllLevels(ThisFlag);
    mrModelPart.RemoveElementsFromAllLevels(ThisFlag);
    mrModelPart.RemoveConditionsFromAllLevels(ThisFlag);
}


/// Execute the refinement once
void UniformRefinementUtility::ExecuteDivision(
    const int& rDivision,
    IndexIndexVectorMapType& rTagNodes,
    IndexIndexVectorMapType& rTagElems,
    IndexIndexVectorMapType& rTagConds
)
{
    // Initialize the auxiliary arrays for the elements and conditions to refine
    ElementsArrayType elements_to_refine;
    ConditionsArrayType conditions_to_refine;

    // Restart the elements and conditions maps (since the preexisting elements will be deleted)
    rTagElems.clear();
    rTagConds.clear();

    // Fill the auxiliary array with the elements to refine
    auto all_elem_begin = mrModelPart.ElementsBegin();
    for (int i = 0; i < static_cast<int>(mrModelPart.Elements().size()); i++)
    {
        auto i_element = all_elem_begin + i;

        // Check the divisions level of the origin elements
        if (i_element->GetValue(NUMBER_OF_DIVISIONS) == rDivision)
            elements_to_refine.push_back(*i_element.base());
    }

    // Fill the auxiliary array with the conditions to refine
    auto all_cond_begin = mrModelPart.ConditionsBegin();
    for (int i = 0; i < static_cast<int>(mrModelPart.Conditions().size()); i++)
    {
        auto i_condition = all_cond_begin + i;

        // Check the refinement level of the origin conditions
        if (i_condition->GetValue(NUMBER_OF_DIVISIONS) == rDivision)
            conditions_to_refine.push_back(*i_condition.base());
    }

    // Loop the origin elements
    ElementsArrayType::iterator elements_begin = elements_to_refine.begin();
    for (int i = 0; i < static_cast<int>(elements_to_refine.size()); i++)
    {
        // Get the element iterator
        auto i_element = elements_begin + i;

        // Get the refinement level of the origin element
        int step_divisions_level = rDivision + 1;

        // Get the geometry
        Geometry<NodeType>& geom = i_element->GetGeometry();

        if (geom.GetGeometryType() == GeometryData::Kratos_Triangle2D3)
        {
            // Initialize the vector of middle nodes
            IndexType i_node = 0;
            std::vector<NodeType::Pointer> middle_nodes(3); // 3 edges

            // Loop the edges to get or create the middle nodes
            for (auto edge : geom.Edges())
                middle_nodes[i_node++] = GetNodeInEdge(edge, step_divisions_level, rTagNodes);
            
            // Split the triangle
            PointerVector<NodeType> sub_element_nodes(3);    // a triangle is defined by 3 nodes
            for (int position = 0; position < 4; position++) // there are 4 sub triangles
            {
                sub_element_nodes = GetSubTriangleNodes(position, geom, middle_nodes);
                CreateElement(i_element, sub_element_nodes, step_divisions_level, rTagElems);
            }
        }
        else if (geom.GetGeometryType() == GeometryData::Kratos_Quadrilateral2D4)
        {
            // Initialize the vector of middle nodes
            IndexType i_node = 0;
            std::vector<NodeType::Pointer> middle_nodes(5); // 4 edges and the quadrilateral itself

            // Loop the edges to get or create the middle nodes
            for (auto edge : geom.Edges())
                middle_nodes[i_node++] = GetNodeInEdge(edge, step_divisions_level, rTagNodes);
            middle_nodes[i_node++] = GetNodeInFace(geom, step_divisions_level, rTagNodes);

            // Split the quadrilateral
            PointerVector<NodeType> sub_element_nodes(4);    // a quadrilateral is defined by 4 nodes
            for (int position = 0; position < 4; position++) // there are 4 sub quadrilateral
            {
                sub_element_nodes = GetSubQuadrilateralNodes(position, geom, middle_nodes);
                CreateElement(i_element, sub_element_nodes, step_divisions_level, rTagElems);
            }
        }
        else if (geom.GetGeometryType() == GeometryData::Kratos_Tetrahedra3D4)
        {
            // Initialize the vector of middle nodes
            IndexType i_node = 0;
            std::vector<NodeType::Pointer> middle_nodes(6); // 6 edges

            // Loop the edges to get or create the middle nodes
            for (auto edge : geom.Edges())
                middle_nodes[i_node++] = GetNodeInEdge(edge, step_divisions_level, rTagNodes);

            // Split the tetrahedra
            PointerVector<NodeType> sub_element_nodes(4);    // a tetrahedra is defined by 4 nodes
            for (int position = 0; position < 8; position++) // there are 8 sub tetrahedrons
            {
                sub_element_nodes = GetSubTetrahedraNodes(position, geom, middle_nodes);
                CreateElement(i_element, sub_element_nodes, step_divisions_level, rTagElems);
            }
        }
        else if (geom.GetGeometryType() == GeometryData::Kratos_Hexahedra3D8)
        {
            // Initialize the vector of middle nodes
            IndexType i_node = 0;
            std::vector<NodeType::Pointer> middle_nodes(19); // 12 edges, 6 faces and the hexahedra itself

            // Loop the edges to get or create the middle nodes
            for (auto edge : geom.Edges())
                middle_nodes[i_node++] = GetNodeInEdge(edge, step_divisions_level, rTagNodes);
            for (auto face : geom.Faces())
                middle_nodes[i_node++] = GetNodeInFace(face, step_divisions_level, rTagNodes);
            middle_nodes[i_node++] = GetNodeInBody(geom, step_divisions_level, rTagNodes);

            // Split the hexahedra
            PointerVector<NodeType> sub_element_nodes(8);    // an hexahedra is defined by 8 nodes
            for (int position = 0; position < 8; position++) // there are 8 sub hexahedrons
            {
                sub_element_nodes = GetSubHexahedraNodes(position, geom, middle_nodes);
                CreateElement(i_element, sub_element_nodes, step_divisions_level, rTagElems);
            }
        }
        else
        {
            KRATOS_ERROR << "Your geometry contains " << geom.GetGeometryType() << " which cannot be refined" << std::endl;
        }

        // Once we have created all the sub elements, the origin element must be deleted
        i_element->Set(TO_ERASE, true);
    }

    mrModelPart.RemoveElementsFromAllLevels(TO_ERASE);

    // Loop the origin conditions
    ConditionsArrayType::iterator conditions_begin = conditions_to_refine.begin();
    for (int i = 0; i < static_cast<int>(conditions_to_refine.size()); i++)
    {
        // Get the condition iterator
        auto i_condition = conditions_begin + i;

        // Get the refinement level of the origin condition
        int step_divisions_level = rDivision + 1;

        // Get the geometry
        Geometry<NodeType>& geom = i_condition->GetGeometry();

        if (geom.GetGeometryType() == GeometryData::Kratos_Line2D2)
        {
            NodeType::Pointer middle_node = GetNodeInEdge(geom, step_divisions_level, rTagNodes);

            // Create the sub conditions
            PointerVector<NodeType> sub_condition_nodes(2);
            for (int position = 0; position < 2; position++)
            {
                sub_condition_nodes = GetSubLineNodes(position, geom, middle_node);
                CreateCondition(i_condition, sub_condition_nodes, step_divisions_level, rTagConds);
            }
        }
        else if (geom.GetGeometryType() == GeometryData::Kratos_Triangle3D3)
        {
            // Initialize the middle nodes vector
            IndexType i_node = 0;
            std::vector<NodeType::Pointer> middle_nodes(3);
            // Loop the edges to get or create the middle nodes
            for (auto edge : geom.Edges())
                middle_nodes[i_node++] = GetNodeInEdge(edge, step_divisions_level, rTagNodes);

            PointerVector<NodeType> sub_condition_nodes(3);    // a triangle is defined by 3 nodes
            for (int position = 0; position < 4; position++) // there are 4 sub triangles
            {
                sub_condition_nodes = GetSubTriangleNodes(position, geom, middle_nodes);
                CreateCondition(i_condition, sub_condition_nodes, step_divisions_level, rTagConds);
            }
        }
        else if (geom.GetGeometryType() == GeometryData::Kratos_Quadrilateral3D4)
        {
            // Initialize the middle nodes vector
            IndexType i_node = 0;
            std::vector<NodeType::Pointer> middle_nodes(5);
            // Loop the edges to get or create the middle nodes
            for (auto edge : geom.Edges())
                middle_nodes[i_node++] = GetNodeInEdge(edge, step_divisions_level, rTagNodes);
            middle_nodes[i_node++] = GetNodeInFace(geom, step_divisions_level, rTagNodes);

            PointerVector<NodeType> sub_condition_nodes(4);    // a quadrilateral is defined by 4 nodes
            for (int position = 0; position < 4; position++) // there are 4 sub quadrilaterals
            {
                sub_condition_nodes = GetSubQuadrilateralNodes(position, geom, middle_nodes);
                CreateCondition(i_condition, sub_condition_nodes, step_divisions_level, rTagConds);
            }
        }
        else
        {
            KRATOS_ERROR << "Your geometry contains " << geom.GetGeometryType() << " which cannot be refined" << std::endl;
        }

        // Once we have created all the sub conditions, the origin conditions must be deleted
        i_condition->Set(TO_ERASE, true);
    }

    mrModelPart.RemoveConditionsFromAllLevels(TO_ERASE);

}


/// Get the middle node on an edge
typename NodeType::Pointer UniformRefinementUtility::GetNodeInEdge(
    const EdgeType& rEdge,
    const int& rNumberOfDivisions,
    IndexIndexVectorMapType& rTagNodes
)
{
    // Initialize the output
    NodeType::Pointer middle_node;

    // Get the middle node key
    EdgeKeyType node_key;
    node_key = std::minmax(rEdge(0)->Id(), rEdge(1)->Id());

    // Check if the node exist
    auto search = mNodesMap.find(node_key);
    if (search != mNodesMap.end() )
    {
        middle_node = mrModelPart.Nodes()(search->second);
    }
    else
    {
        middle_node = CreateNodeInEdge(rEdge, rNumberOfDivisions, node_key, rTagNodes);
    }

    return middle_node;
}


/// Create a middle node on an edge. If the node does not exist, it creates one
typename NodeType::Pointer UniformRefinementUtility::CreateNodeInEdge(
    const EdgeType& rEdge,
    const int& rNumberOfDivisions,
    const EdgeKeyType& rNodeKey,
    IndexIndexVectorMapType& rTagNodes
)
{
    // Initialize the output
    NodeType::Pointer middle_node;

    // Create the new node
    const double new_x = 0.5*rEdge(0)->X() + 0.5*rEdge(1)->X();
    const double new_y = 0.5*rEdge(0)->Y() + 0.5*rEdge(1)->Y();
    const double new_z = 0.5*rEdge(0)->Z() + 0.5*rEdge(1)->Z();
    middle_node = mrModelPart.CreateNewNode(++mLastNodeId, new_x, new_y, new_z);

    // Store the node key in the map
    mNodesMap[rNodeKey] = middle_node->Id();

    // interpolate the variables
    CalculateNodalStepData(middle_node, rEdge(0), rEdge(1));

    // Set the number of divisions level
    int& this_node_level = middle_node->GetValue(NUMBER_OF_DIVISIONS);
    this_node_level = rNumberOfDivisions;

    // Set the appropriate flags
    middle_node->Set(NEW_ENTITY, true);

    // Set the DoF's
    for (typename NodeType::DofsContainerType::const_iterator it_dof = mDofs.begin(); it_dof != mDofs.end(); ++it_dof)
        middle_node->pAddDof(*it_dof);

    // Store the created node on the taps map in order to later add it to the sub model parts
    IndexType tag = mNodesTags[rEdge(0)->Id()];
    rTagNodes[tag].push_back(middle_node->Id());
    mNodesTags[middle_node->Id()] = tag;

    return middle_node;
}


/// Get the middle node on a face defined by four nodes. If the node does not exist, it creates one
typename NodeType::Pointer UniformRefinementUtility::GetNodeInFace(
    const FaceType& rFace,
    const int& rNumberOfDivisions,
    IndexIndexVectorMapType& rTagNodes
)
{
    // Initialize the output
    NodeType::Pointer middle_node;

    // Get the middle node key
    FaceKeyType node_key = {{rFace(0)->Id(), rFace(1)->Id(), rFace(2)->Id(), rFace(3)->Id()}};
    std::sort(node_key.begin(), node_key.end());

    // Check if the node exist
    auto search = mNodesInFaceMap.find(node_key);
    if (search != mNodesInFaceMap.end() )
    {
        middle_node = mrModelPart.Nodes()(search->second);
    }
    else
    {
        middle_node = CreateNodeInFace(rFace, rNumberOfDivisions, node_key, rTagNodes);
    }

    return middle_node;
}


/// Get the middle node on a face defined by four nodes. If the node does not exist, it creates one
typename NodeType::Pointer UniformRefinementUtility::CreateNodeInFace(
    const FaceType& rFace,
    const int& rNumberOfDivisions,
    const FaceKeyType& rNodeKey,
    IndexIndexVectorMapType& rTagNodes
)
{
    // Initialize the output
    NodeType::Pointer middle_node;

    // Create the new node
    const double new_x = 0.25*rFace(0)->X() + 0.25*rFace(1)->X() + 0.25*rFace(2)->X() + 0.25*rFace(3)->X();
    const double new_y = 0.25*rFace(0)->Y() + 0.25*rFace(1)->Y() + 0.25*rFace(2)->Y() + 0.25*rFace(3)->Y();
    const double new_z = 0.25*rFace(0)->Z() + 0.25*rFace(1)->Z() + 0.25*rFace(2)->Z() + 0.25*rFace(3)->Z();
    middle_node = mrModelPart.CreateNewNode(++mLastNodeId, new_x, new_y, new_z);

    // Store the node key in the map
    mNodesInFaceMap[rNodeKey] = middle_node->Id();

    // interpolate the variables
    CalculateNodalStepData(middle_node, rFace(0), rFace(1), rFace(2), rFace(3));

    // Set the refinement level
    int& this_node_level = middle_node->GetValue(NUMBER_OF_DIVISIONS);
    this_node_level = rNumberOfDivisions;

    // Set the appropriate flags
    middle_node->Set(NEW_ENTITY, true);

    // Set the DoF's
    for (typename NodeType::DofsContainerType::const_iterator it_dof = mDofs.begin(); it_dof != mDofs.end(); ++it_dof)
        middle_node->pAddDof(*it_dof);

    // Store the created node on the tags map in order to later add it to the sub model parts
    IndexType tag = mNodesTags[rFace(0)->Id()];
    rTagNodes[tag].push_back(middle_node->Id());
    mNodesTags[middle_node->Id()] = tag;

    return middle_node;
}


/// Create the middle node on a body defined by eight nodes.
typename NodeType::Pointer UniformRefinementUtility::GetNodeInBody(
    const BodyType& rBody,
    const int& rNumberOfDivisions,
    IndexIndexVectorMapType& rTagNodes
)
{
    // Initialize the output
    NodeType::Pointer middle_node;

    // Create the new node
    const double new_x = 0.125*rBody(0)->X() + 0.125*rBody(1)->X() + 0.125*rBody(2)->X() + 
        0.125*rBody(3)->X() + 0.125*rBody(4)->X() + 0.125*rBody(5)->X() + 0.125*rBody(6)->X() + 0.125*rBody(7)->X();
    const double new_y = 0.125*rBody(0)->Y() + 0.125*rBody(1)->Y() + 0.125*rBody(2)->Y() + 
        0.125*rBody(3)->Y() + 0.125*rBody(4)->Y() + 0.125*rBody(5)->Y() + 0.125*rBody(6)->Y() + 0.125*rBody(7)->Y();
    const double new_z = 0.125*rBody(0)->Z() + 0.125*rBody(1)->Z() + 0.125*rBody(2)->Z() + 
        0.125*rBody(3)->Z() + 0.125*rBody(4)->Z() + 0.125*rBody(5)->Z() + 0.125*rBody(6)->Z() + 0.125*rBody(7)->Z();
    middle_node = mrModelPart.CreateNewNode(++mLastNodeId, new_x, new_y, new_z);

    // Interpolate the variables
    CalculateNodalStepData(middle_node, rBody);

    // Set the refinement level
    int& this_node_level = middle_node->GetValue(NUMBER_OF_DIVISIONS);
    this_node_level = rNumberOfDivisions;

    // Set the appropriate flags
    middle_node->Set(NEW_ENTITY, true);

    // Set the DoF's
    for (typename NodeType::DofsContainerType::const_iterator it_dof = mDofs.begin(); it_dof != mDofs.end(); ++it_dof)
        middle_node->pAddDof(*it_dof);

    // Store the created node on the tags map in order to later add it to the sub model parts
    IndexType tag = mNodesTags[rBody(0)->Id()];
    rTagNodes[tag].push_back(middle_node->Id());
    mNodesTags[middle_node->Id()] = tag;

    return middle_node;
}

/// Compute the nodal data of a node
void UniformRefinementUtility::CalculateNodalStepData(
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

    WeakPointerVector<NodeType>& r_new_father_nodes = pNewNode->GetValue(FATHER_NODES);
    r_new_father_nodes.clear();
    r_new_father_nodes = pNode0->GetValue(FATHER_NODES);

    std::vector<double>& r_new_father_nodes_weights = pNewNode->GetValue(FATHER_NODES_WEIGHTS);
    r_new_father_nodes_weights.clear();
    r_new_father_nodes_weights = pNode0->GetValue(FATHER_NODES_WEIGHTS);

    AddOtherFatherNodes(r_new_father_nodes, r_new_father_nodes_weights,
        pNode1->GetValue(FATHER_NODES), pNode1->GetValue(FATHER_NODES_WEIGHTS));
}


/// Compute the nodal data of a node
void UniformRefinementUtility::CalculateNodalStepData(
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

    WeakPointerVector<NodeType>& r_new_father_nodes = pNewNode->GetValue(FATHER_NODES);
    r_new_father_nodes.clear();
    r_new_father_nodes = pNode0->GetValue(FATHER_NODES);

    std::vector<double>& r_new_father_nodes_weights = pNewNode->GetValue(FATHER_NODES_WEIGHTS);
    r_new_father_nodes_weights.clear();
    r_new_father_nodes_weights = pNode0->GetValue(FATHER_NODES_WEIGHTS);

    AddOtherFatherNodes(r_new_father_nodes, r_new_father_nodes_weights,
        pNode1->GetValue(FATHER_NODES), pNode1->GetValue(FATHER_NODES_WEIGHTS), 0.5);
    AddOtherFatherNodes(r_new_father_nodes, r_new_father_nodes_weights,
        pNode1->GetValue(FATHER_NODES), pNode2->GetValue(FATHER_NODES_WEIGHTS), 1/3);
    AddOtherFatherNodes(r_new_father_nodes, r_new_father_nodes_weights,
        pNode1->GetValue(FATHER_NODES), pNode3->GetValue(FATHER_NODES_WEIGHTS), 0.25);
}


/// Compute the nodal data of a node
void UniformRefinementUtility::CalculateNodalStepData(
    NodeType::Pointer pNewNode,
    const BodyType& rBody
)
{
    FaceKeyType key;
    // Get the node in the center of the first face
    key = {rBody(0)->Id(), rBody(1)->Id(), rBody(2)->Id(), rBody(3)->Id()};
    std::sort(key.begin(), key.end());
    NodeType::Pointer node_0 = mrModelPart.pGetNode(mNodesInFaceMap[key]);
    // Get the node in the center of the opposite face
    key = {rBody(4)->Id(), rBody(5)->Id(), rBody(6)->Id(), rBody(7)->Id()};
    std::sort(key.begin(), key.end());
    NodeType::Pointer node_1 = mrModelPart.pGetNode(mNodesInFaceMap[key]);
    // Compute the data as an average of this two nodes
    CalculateNodalStepData(pNewNode, node_0, node_1);
}


/// Add the father nodes which does not exist in the current father nodes
void UniformRefinementUtility::AddOtherFatherNodes(
    WeakPointerVector<NodeType>& rThisFatherNodes,
    std::vector<double>& rThisFatherWeights,
    WeakPointerVector<NodeType>& rOtherFatherNodes,
    const std::vector<double>& rOtherFatherWeights,
    const double& rWeight
)
{
    for (auto& weight : rThisFatherWeights)
        weight *= (1-rWeight);

    WeakPointerVector<NodeType>::iterator other_nodes_begin = rOtherFatherNodes.begin();
    for (IndexType o = 0; o < rOtherFatherNodes.size(); o++)
    {
        auto other_node = other_nodes_begin + o;
        bool other_not_found = true;

        WeakPointerVector<NodeType>::iterator this_nodes_begin = rThisFatherNodes.begin();
        for (IndexType t = 0; (t < rThisFatherNodes.size()) && (other_not_found); t++)
        {
            auto this_node = this_nodes_begin + t;
            if (other_node->Id() == this_node->Id())
            {
                rThisFatherWeights[t] = rOtherFatherWeights[o] * rWeight;
                other_not_found = false;
            }
        }
        if (other_not_found)
        {
            rThisFatherNodes.push_back(*other_node.base());
            rThisFatherWeights.push_back(rOtherFatherWeights[o] * rWeight);
        }
    }
}


/// Create a sub element
void UniformRefinementUtility::CreateElement(
    ElementsArrayType::iterator pOriginElement,
    PointerVector<NodeType>& rThisNodes,
    const int& rNumberOfDivisions,
    IndexIndexVectorMapType& rTagElems
)
{
    Element::Pointer sub_element = pOriginElement->Clone(++mLastElemId, rThisNodes);

    if (sub_element != nullptr)
    {
        // Add the element to the origin model part
        mrModelPart.AddElement(sub_element);

        // Set the refinement level
        int& this_elem_level = sub_element->GetValue(NUMBER_OF_DIVISIONS);
        this_elem_level = rNumberOfDivisions;

        // Store the father element pointer
        sub_element->SetValue(FATHER_ELEMENT, pOriginElement->GetValue(FATHER_ELEMENT));

        // Store the created element on the tags map in order to later add it to the sub model parts
        IndexType tag = mElementsTags[pOriginElement->Id()];
        rTagElems[tag].push_back(sub_element->Id());
        mElementsTags[sub_element->Id()] = tag;
    }
}


/// Create a sub condition
void UniformRefinementUtility::CreateCondition(
    ConditionsArrayType::iterator pOriginCondition,
    PointerVector<NodeType>& rThisNodes,
    const int& rNumberOfDivisions,
    IndexIndexVectorMapType& rTagConds
)
{
    Condition::Pointer sub_condition = pOriginCondition->Clone(++mLastCondId, rThisNodes);

    if (sub_condition != nullptr)
    {
        // Add the condition to the origin model part
        mrModelPart.AddCondition(sub_condition);

        // Set the refinement level
        int& this_cond_level = sub_condition->GetValue(NUMBER_OF_DIVISIONS);
        this_cond_level = rNumberOfDivisions;

        // Store the created condition on the tags map in order to later add it to the sub model parts
        IndexType tag = mConditionsTags[pOriginCondition->Id()];
        rTagConds[tag].push_back(sub_condition->Id());
        mConditionsTags[sub_condition->Id()] = tag;
    }
}


/// Return the nodes defining the i-subline
PointerVector<NodeType> UniformRefinementUtility::GetSubLineNodes(
    const int Position,
    const Geometry<NodeType>& rGeom,
    NodeType::Pointer& rMiddleNode
)
{
    PointerVector<NodeType> sub_line_nodes(2);

    if (Position == 0)
    {
        // First sub line
        sub_line_nodes(0) = rGeom.pGetPoint(0);
        sub_line_nodes(1) = rMiddleNode;
    }
    else if (Position == 1)
    {
        // second sub line
        sub_line_nodes(0) = rMiddleNode;
        sub_line_nodes(1) = rGeom.pGetPoint(1);
    }
    else
    {
        KRATOS_ERROR << "Attempting to get " << Position << " sub-line inside a line" << std::endl;
    }

    return sub_line_nodes;
}


/// Return the nodes defining the i-subtriangle
PointerVector<NodeType> UniformRefinementUtility::GetSubTriangleNodes(
    const int Position,
    const Geometry<NodeType>& rGeom,
    std::vector<NodeType::Pointer>& rMiddleNodes
)
{
    PointerVector<NodeType> sub_triangle_nodes(3);

    if (Position == 0)
    {
        // First sub triangle
        sub_triangle_nodes(0) = rGeom.pGetPoint(0);
        sub_triangle_nodes(1) = rMiddleNodes[0];
        sub_triangle_nodes(2) = rMiddleNodes[2];
    }
    else if (Position == 1)
    {
        // Second sub triangle
        sub_triangle_nodes(0) = rGeom.pGetPoint(1);
        sub_triangle_nodes(1) = rMiddleNodes[1];
        sub_triangle_nodes(2) = rMiddleNodes[0];
    }
    else if (Position == 2)
    {
        // Third sub triangle
        sub_triangle_nodes(0) = rGeom.pGetPoint(2);
        sub_triangle_nodes(1) = rMiddleNodes[2];
        sub_triangle_nodes(2) = rMiddleNodes[1];
    }
    else if (Position == 3)
    {
        // Fourth sub triangle (inner triangle)
        sub_triangle_nodes(0) = rMiddleNodes[0];
        sub_triangle_nodes(1) = rMiddleNodes[1];
        sub_triangle_nodes(2) = rMiddleNodes[2];
    }
    else
    {
        KRATOS_ERROR << "Attempting to get " << Position << " sub-triangle inside a triangle" << std::endl;
    }

    return sub_triangle_nodes;
}

/// Return the nodes defining the i-subquadrilateral
PointerVector<NodeType> UniformRefinementUtility::GetSubQuadrilateralNodes(
    const int Position,
    const Geometry<NodeType>& rGeom,
    std::vector<NodeType::Pointer>& rMiddleNodes
)
{
    PointerVector<NodeType> sub_quadrilateral_nodes(4);

    if (Position == 0)
    {
        // First sub element
        sub_quadrilateral_nodes(0) = rGeom.pGetPoint(0);
        sub_quadrilateral_nodes(1) = rMiddleNodes[0];
        sub_quadrilateral_nodes(2) = rMiddleNodes[4];
        sub_quadrilateral_nodes(3) = rMiddleNodes[3];
    }
    else if (Position == 1)
    {
        // Second sub element
        sub_quadrilateral_nodes(0) = rGeom.pGetPoint(1);
        sub_quadrilateral_nodes(1) = rMiddleNodes[1];
        sub_quadrilateral_nodes(2) = rMiddleNodes[4];
        sub_quadrilateral_nodes(3) = rMiddleNodes[0];
    }
    else if (Position == 2)
    {
        // Third sub element
        sub_quadrilateral_nodes(0) = rGeom.pGetPoint(2);
        sub_quadrilateral_nodes(1) = rMiddleNodes[2];
        sub_quadrilateral_nodes(2) = rMiddleNodes[4];
        sub_quadrilateral_nodes(3) = rMiddleNodes[1];
    }
    else if (Position == 3)
    {
        // Fourth sub element
        sub_quadrilateral_nodes(0) = rGeom.pGetPoint(3);
        sub_quadrilateral_nodes(1) = rMiddleNodes[3];
        sub_quadrilateral_nodes(2) = rMiddleNodes[4];
        sub_quadrilateral_nodes(3) = rMiddleNodes[2];
    }
    else
    {
        KRATOS_ERROR << "Attempting to get " << Position << " sub-quadrilateral inside a quadrilateral" << std::endl;
    }

    return sub_quadrilateral_nodes;
}


/// Return the nodes defining the i-tetrahedra
PointerVector<NodeType> UniformRefinementUtility::GetSubTetrahedraNodes(
    const int Position,
    const Geometry<NodeType>& rGeom,
    std::vector<NodeType::Pointer>& rMiddleNodes
)
{
    PointerVector<NodeType> sub_tetrahedra_nodes(4);

    if (Position == 0)
    {
        // First sub element
        sub_tetrahedra_nodes(0) = rGeom.pGetPoint(0);
        sub_tetrahedra_nodes(1) = rMiddleNodes[0];
        sub_tetrahedra_nodes(2) = rMiddleNodes[2];
        sub_tetrahedra_nodes(3) = rMiddleNodes[3];
    }
    else if (Position == 1)
    {
        // Second sub element
        sub_tetrahedra_nodes(0) = rMiddleNodes[0];
        sub_tetrahedra_nodes(1) = rGeom.pGetPoint(1);
        sub_tetrahedra_nodes(2) = rMiddleNodes[1];
        sub_tetrahedra_nodes(3) = rMiddleNodes[4];
    }
    else if (Position == 2)
    {
        // Third sub element
        sub_tetrahedra_nodes(0) = rMiddleNodes[2];
        sub_tetrahedra_nodes(1) = rMiddleNodes[1];
        sub_tetrahedra_nodes(2) = rGeom.pGetPoint(2);
        sub_tetrahedra_nodes(3) = rMiddleNodes[5];
    }
    else if (Position == 3)
    {
        // Fourth sub element
        sub_tetrahedra_nodes(0) = rMiddleNodes[3];
        sub_tetrahedra_nodes(1) = rMiddleNodes[4];
        sub_tetrahedra_nodes(2) = rMiddleNodes[5];
        sub_tetrahedra_nodes(3) = rGeom.pGetPoint(3);
    }
    else if (Position == 4)
    {
        // Fifth sub element (inner element: Mesh quality is not conserved)
        sub_tetrahedra_nodes(0) = rMiddleNodes[0];
        sub_tetrahedra_nodes(1) = rMiddleNodes[1];
        sub_tetrahedra_nodes(2) = rMiddleNodes[2];
        sub_tetrahedra_nodes(3) = rMiddleNodes[3];
    }
    else if (Position == 5)
    {
        // Sixth sub element (inner element: Mesh quality is not conserved)
        sub_tetrahedra_nodes(0) = rMiddleNodes[2];
        sub_tetrahedra_nodes(1) = rMiddleNodes[3];
        sub_tetrahedra_nodes(2) = rMiddleNodes[5];
        sub_tetrahedra_nodes(3) = rMiddleNodes[4];
    }
    else if (Position == 6)
    {
        // Seventh sub element (inner element: Mesh quality is not conserved)
        sub_tetrahedra_nodes(0) = rMiddleNodes[0];
        sub_tetrahedra_nodes(1) = rMiddleNodes[3];
        sub_tetrahedra_nodes(2) = rMiddleNodes[4];
        sub_tetrahedra_nodes(3) = rMiddleNodes[1];
    }
    else if (Position == 7)
    {
        // Eighth sub element (inner element: Mesh quality is not conserved)
        sub_tetrahedra_nodes(0) = rMiddleNodes[1];
        sub_tetrahedra_nodes(1) = rMiddleNodes[4];
        sub_tetrahedra_nodes(2) = rMiddleNodes[5];
        sub_tetrahedra_nodes(3) = rMiddleNodes[2];
    }
    else
    {
        KRATOS_ERROR << "Attempting to get " << Position << " sub-tetrahedra inside a tetrahedra" << std::endl;
    }

    return sub_tetrahedra_nodes;
}

/// Return the nodes defining the i-tetrahedra
PointerVector<NodeType> UniformRefinementUtility::GetSubHexahedraNodes(
    const int Position,
    const Geometry<NodeType>& rGeom,
    std::vector<NodeType::Pointer>& rMiddleNodes
)
{
    PointerVector<NodeType> sub_hexahedra_nodes(8);

    if (Position == 0)
    {
        // Firs sub element
        sub_hexahedra_nodes(0) = rGeom.pGetPoint(0);
        sub_hexahedra_nodes(1) = rMiddleNodes[0];
        sub_hexahedra_nodes(2) = rMiddleNodes[12];
        sub_hexahedra_nodes(3) = rMiddleNodes[3];
        sub_hexahedra_nodes(4) = rMiddleNodes[8];
        sub_hexahedra_nodes(5) = rMiddleNodes[13];
        sub_hexahedra_nodes(6) = rMiddleNodes[18];
        sub_hexahedra_nodes(7) = rMiddleNodes[16];
    }
    else if (Position == 1)
    {
        // Second sub element
        sub_hexahedra_nodes(0) = rMiddleNodes[0];
        sub_hexahedra_nodes(1) = rGeom.pGetPoint(1);
        sub_hexahedra_nodes(2) = rMiddleNodes[1];
        sub_hexahedra_nodes(3) = rMiddleNodes[12];
        sub_hexahedra_nodes(4) = rMiddleNodes[13];
        sub_hexahedra_nodes(5) = rMiddleNodes[9];
        sub_hexahedra_nodes(6) = rMiddleNodes[14];
        sub_hexahedra_nodes(7) = rMiddleNodes[18];
    }
    else if (Position == 2)
    {
        // Third sub element
        sub_hexahedra_nodes(0) = rMiddleNodes[12];
        sub_hexahedra_nodes(1) = rMiddleNodes[1];
        sub_hexahedra_nodes(2) = rGeom.pGetPoint(2);
        sub_hexahedra_nodes(3) = rMiddleNodes[2];
        sub_hexahedra_nodes(4) = rMiddleNodes[18];
        sub_hexahedra_nodes(5) = rMiddleNodes[14];
        sub_hexahedra_nodes(6) = rMiddleNodes[10];
        sub_hexahedra_nodes(7) = rMiddleNodes[15];
    }
    else if (Position == 3)
    {
        // Fourth sub element
        sub_hexahedra_nodes(0) = rMiddleNodes[3];
        sub_hexahedra_nodes(1) = rMiddleNodes[12];
        sub_hexahedra_nodes(2) = rMiddleNodes[2];
        sub_hexahedra_nodes(3) = rGeom.pGetPoint(3);
        sub_hexahedra_nodes(4) = rMiddleNodes[16];
        sub_hexahedra_nodes(5) = rMiddleNodes[18];
        sub_hexahedra_nodes(6) = rMiddleNodes[15];
        sub_hexahedra_nodes(7) = rMiddleNodes[11];
    }
    else if (Position == 4)
    {
        // Fifth sub element
        sub_hexahedra_nodes(0) = rMiddleNodes[8];
        sub_hexahedra_nodes(1) = rMiddleNodes[13];
        sub_hexahedra_nodes(2) = rMiddleNodes[18];
        sub_hexahedra_nodes(3) = rMiddleNodes[16];
        sub_hexahedra_nodes(4) = rGeom.pGetPoint(4);
        sub_hexahedra_nodes(5) = rMiddleNodes[4];
        sub_hexahedra_nodes(6) = rMiddleNodes[17];
        sub_hexahedra_nodes(7) = rMiddleNodes[7];
    }
    else if (Position == 5)
    {
        // Sixth sub element
        sub_hexahedra_nodes(0) = rMiddleNodes[13];
        sub_hexahedra_nodes(1) = rMiddleNodes[9];
        sub_hexahedra_nodes(2) = rMiddleNodes[14];
        sub_hexahedra_nodes(3) = rMiddleNodes[18];
        sub_hexahedra_nodes(4) = rMiddleNodes[4];
        sub_hexahedra_nodes(5) = rGeom.pGetPoint(5);
        sub_hexahedra_nodes(6) = rMiddleNodes[5];
        sub_hexahedra_nodes(7) = rMiddleNodes[17];
    }
    else if (Position == 6)
    {
        // Seventh sub element
        sub_hexahedra_nodes(0) = rMiddleNodes[18];
        sub_hexahedra_nodes(1) = rMiddleNodes[14];
        sub_hexahedra_nodes(2) = rMiddleNodes[10];
        sub_hexahedra_nodes(3) = rMiddleNodes[15];
        sub_hexahedra_nodes(4) = rMiddleNodes[17];
        sub_hexahedra_nodes(5) = rMiddleNodes[5];
        sub_hexahedra_nodes(6) = rGeom.pGetPoint(6);
        sub_hexahedra_nodes(7) = rMiddleNodes[6];
    }
    else if (Position == 7)
    {
        // Eighth sub element
        sub_hexahedra_nodes(0) = rMiddleNodes[16];
        sub_hexahedra_nodes(1) = rMiddleNodes[18];
        sub_hexahedra_nodes(2) = rMiddleNodes[15];
        sub_hexahedra_nodes(3) = rMiddleNodes[11];
        sub_hexahedra_nodes(4) = rMiddleNodes[7];
        sub_hexahedra_nodes(5) = rMiddleNodes[17];
        sub_hexahedra_nodes(6) = rMiddleNodes[6];
        sub_hexahedra_nodes(7) = rGeom.pGetPoint(7);
    }
    else
    {
        KRATOS_ERROR << "Attempting to get " << Position << " sub-hexahedra inside a hexahedra" << std::endl;
    }

    return sub_hexahedra_nodes;
}


}  // namespace Kratos.
