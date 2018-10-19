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


#if !defined( KRATOS_UNIFORM_REFINEMENT_UTILITY_H_INCLUDED )
#define KRATOS_UNIFORM_REFINEMENT_UTILITY_H_INCLUDED


// System includes
#include <string>
#include <iostream>
#include <unordered_map>


// External includes


// Project includes
#include "includes/key_hash.h"
#include "includes/define.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/model_part.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "meshing_application.h"


namespace Kratos
{
///@addtogroup MeshingApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Divide the elements to uniformly refine a mesh
/**
 *  A node is added on each element edge (an additional node is added inside quadrilaterals
 *  and tetrahedrons) to split the elements.
 *  If a higher refinement is needed, the utility can be called recursively.
 */
class UniformRefinementUtility
{
public:
    ///@name Type Definitions
    ///@{

    /**
     * Node type
     */
    typedef Node<3> NodeType;

    /**
     * Type of edge geometry
     */
    typedef Line3D2<NodeType> EdgeType;

    /**
     * Type of face geometry
     */
    typedef Quadrilateral3D4<NodeType> FaceType;

    /**
     * Type of body geometry
     */
    typedef Hexahedra3D8<NodeType> BodyType;

    /**
     * Type of IDs
     */
    typedef std::size_t IndexType;

    /**
     * Nodes, elements and conditions containers definition
     */
    typedef ModelPart::NodesContainerType NodesArrayType;
    typedef ModelPart::ElementsContainerType ElementsArrayType;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;

    /**
     * Key and map types to locate nodes in the mesh
     */
    typedef std::pair<IndexType, IndexType> EdgeKeyType;
    typedef std::array<IndexType, 4> FaceKeyType;
    typedef std::map<EdgeKeyType, IndexType> NodesInEdgeMapType;
    typedef std::unordered_map<FaceKeyType, IndexType, KeyHasherRange<FaceKeyType>, KeyComparorRange<FaceKeyType>> NodesInFaceMapType;

    /**
     * Map types for the AssignUniqueModelPartCollectionTagUtility
     */
    typedef std::vector<IndexType> IndexVectorType;
    typedef std::vector<std::string> StringVectorType;
    typedef std::unordered_map<IndexType, IndexType> IndexIndexMapType;
    typedef std::unordered_map<IndexType, IndexVectorType> IndexIndexVectorMapType;
    typedef std::unordered_map<IndexType, StringVectorType> IndexStringVectorMapType;

    /// Pointer definition of UniformRefinementUtility
    KRATOS_CLASS_POINTER_DEFINITION(UniformRefinementUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    UniformRefinementUtility(ModelPart& rModelPart);

    /// Destructor.
    virtual ~UniformRefinementUtility();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Execute the refinement until the final refinement is reached
     */
    void Refine(int& rFinalRefinementLevel);

    /**
     * @brief Set the custom ids which will be used to create new entities
     * @param rNodeId
     * @param rElemId
     * @param rCondId
     */
    void SetCustomIds(IndexType& rNodeId, IndexType& rElemId, IndexType& rCondId);

    /**
     * @brief Get the las id of the created nodes, elements and conditions
     * @param rNodeId
     * @param rElemId
     * @param rCondId
     */
    void GetLastCreatedIds(IndexType& rNodeId, IndexType& rElemId, IndexType& rCondId);

    /**
     * @brief RemoveRefinedEntities
     * @param ThisFlag
     */
    void RemoveRefinedEntities(Flags ThisFlag);

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const;


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;             /// The model part to refine
    IndexType mDimension;

    IndexType mLastNodeId;           /// The node Id
    IndexType mLastElemId;           /// The element Id
    IndexType mLastCondId;           /// The condition Id
    IndexType mStepDataSize;         /// The size of the nodal database
    IndexType mBufferSize;           /// The buffer size
    NodeType::DofsContainerType mDofs;  /// Storage for the dof of the node

    NodesInEdgeMapType mNodesMap;              /// Where the father nodes IDs are stored
    NodesInFaceMapType mNodesInFaceMap;        /// Where the father nodes IDs are stored

    IndexIndexMapType mNodesTags;
    IndexIndexMapType mElementsTags;
    IndexIndexMapType mConditionsTags;
    IndexStringVectorMapType mCollections;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief ExecuteDivision executes the refinement once
     * @detail Only the entities with NUMBER_OF_DIVISIONS = rDivision are refined to rDivision+1
     * @param rDivision The flag
     * @param rTagNodes The map with the vector of Id's associated to the tag
     * @param rTagElems The map with the vector of Id's associated to the tag
     * @param rTagConds The map with the vector of Id's associated to the tag
     */
    void ExecuteDivision(
        const int& rDivision,
        IndexIndexVectorMapType& rTagNodes,
        IndexIndexVectorMapType& rTagElems,
        IndexIndexVectorMapType& rTagConds
        );

    /**
     * @brief GetNodeInEdge gets the middle node on an edge and return a pointer to it
     * @detail CreateNodeInEdge should be executed before to ensure the node existance
     * @param rEdge The edge containing the two father nodes
     * @param rNumberOfDivisions
     * @param rTagNodes  The map with the vector of Id's associated to the tag
     */
    typename NodeType::Pointer GetNodeInEdge(
        const EdgeType& rEdge,
        const int& rNumberOfDivisions,
        IndexIndexVectorMapType& rTagNodes
        );

    /**
     * @brief CreateNodeInEdge creates a node at the middle point of an edge
     * @detail If the middle node does not exist, creates a new one and set the
     * nodal values. Otherwise, do nothing
     * @see GetNodeInEdge
     * @param rEdge The edge containing the two father nodes
     * @param rNumberOfDivisions The value to set NUMBER_OF_DIVISIONS flag
     * @param rNodeKey The pair of Id's of the edge nodes
     * @param rTagNodes The map with the vector of Id's associated to the tag
     */
    typename NodeType::Pointer CreateNodeInEdge(
        const EdgeType& rEdge,
        const int& rNumberOfDivisions,
        const EdgeKeyType& rNodeKey,
        IndexIndexVectorMapType& rTagNodes
        );

    /**
     * @brief GetNodeInFace gets the node inside a face
     * @detail If the middle node does not exist, create a new one and returns a pointer to it
     * @param rFace The face containing the father nodes 
     * @param rNumberOfDivisions The value to set NUMBER_OF_DIVISIONS flag
     * @param rTagNodes The map with the vector of Id's associated to the tag
     */
    typename NodeType::Pointer GetNodeInFace(
        const FaceType& rFace,
        const int& rNumberOfDivisions,
        IndexIndexVectorMapType& rTagNodes
        );

    /**
     * @brief CreateNodeInFace creates a node at the middle point of a face
     * @detail If the middle node does not exist, create a new one. Otherwise, do nothing
     * @see GetNodeInFace
     * @param rFace The face containing the father nodes
     * @param rNumberOfDivisions The value to set NUMBER_OF_DIVISIONS flag
     * @param rNodeKey The pair of Id's of the edge nodes
     * @param rTagNodes The map with the vector of Id's associated to the tag
     */
    typename NodeType::Pointer CreateNodeInFace(
        const FaceType& rFace,
        const int& rNumberOfDivisions,
        const FaceKeyType& rNodeKey,
        IndexIndexVectorMapType& rTagNodes
        );
    
    /**
     * @brief GetNodeInBody gets the node inside an hexahedrons
     * @detail The middle node does never exist, so always creates a new one and returns a pointer to it
     * @param rBody The geometry containing the father nodes 
     * @param rNumberOfDivisions The value to set NUMBER_OF_DIVISIONS flag
     * @param rTagNodes The map with the vector of Id's associated to the tag
     */
    typename NodeType::Pointer GetNodeInBody(
        const BodyType& rBody,
        const int& rNumberOfDivisions,
        IndexIndexVectorMapType& rTagNodes
        );

    /**
     * @brief CalculateNodalStepData calculates the nodal data as the mean of the father nodes
     * @detail The destination node is assumed to be at the mid point between the origin nodes
     * @param pNewNode The destination node
     * @param pNode0 The first origin node
     * @param pNode1 The second origin node
     */
    void CalculateNodalStepData(
        NodeType::Pointer pNewNode,
        const NodeType::Pointer pNode0,
        const NodeType::Pointer pNode1
        );

    /**
     * @brief CalculateNodalStepData calculates the nodal data as the mean of the father nodes
     * @detail The destination node is assumed to be at the mid point among the origin nodes
     * @param pNewNode The destination node
     * @param pNode0 The first origin node
     * @param pNode1 The second origin node
     * @param pNode2 The third origin node
     * @param pNode3 The fourth origin node
     */
    void CalculateNodalStepData(
        NodeType::Pointer pNewNode,
        const NodeType::Pointer pNode0,
        const NodeType::Pointer pNode1,
        const NodeType::Pointer pNode2,
        const NodeType::Pointer pNode3
        );

    /**
     * @brief CalculateNodalStepData calculates the nodal data as the mean of the father nodes
     * @detail The destination node is assumed to be at the mid point among the origin nodes
     * @param pNewNode The destination node
     * @param rBody The array with the 8 origin nodes
     */
    void CalculateNodalStepData(
        NodeType::Pointer pNewNode,
        const BodyType& rBody
        );

    /**
     * @brief Add the other father nodes if they does not exist in the current father nodes vector
     * @see CalculateNodalStepData
     * @param rThisFatherNodes the current father nodes
     * @param rOtherFatherNodes the father nodes to insert
     */
    void AddOtherFatherNodes(
        WeakPointerVector<NodeType>& rThisFatherNodes,
        std::vector<double>& rThisFatherWeights,
        WeakPointerVector<NodeType>& rOtherFatherNodes,
        const std::vector<double>& rOtherFatherWeights,
        const double& rWeight = 0.5
        );

    /**
     * @brief CreateElement creates an element from an origin element
     * @param pOriginElement pointer to the father element
     * @ThisNodes vector containing the sub element nodes
     * @param rNumberOfDivisions The value to set NUMBER_OF_DIVISIONS flag
     */
    void CreateElement(
        ElementsArrayType::iterator pOriginElement,
        PointerVector<NodeType>& rThisNodes,
        const int& rNumberOfDivisions,
        IndexIndexVectorMapType& rTagElems
        );

    /**
     * @brief CreateCondition creates a condition from an origin condition
     * @param pOriginCondition pointer to the father condition
     * @ThisNodes vector containing the sub condition nodes
     * @param rNumberOfDivisions The value to set NUMBER_OF_DIVISIONS flag
     */
    void CreateCondition(
        ConditionsArrayType::iterator pOriginCondition,
        PointerVector<NodeType>& rThisNodes,
        const int& rNumberOfDivisions,
        IndexIndexVectorMapType& rTagConds
        );

    /**
     * @brief GetSubLineNodes gets the connectivity of a sub-line inside a geometry and the new node
     * @param Position The index which defines the sub-line
     * @param rGeom The original geometry
     * @param rMiddleNode The node which divides the original geometry
     * @return The nodes defining the sub-line
     */
    PointerVector<NodeType> GetSubLineNodes(
        const int Position,
        const Geometry<NodeType>& rGeom,
        NodeType::Pointer& rMiddleNode
        );

    /**
     * @brief GetSubTriangleNodes gets the connectivity of a sub-triangle inside a
     * geometry and the new nodes
     * @param Position The index which defines the sub-triangle
     * @param rGeom The original geometry
     * @param rMiddleNode The nodes which divides the original geometry
     * @return The nodes defining the sub-triangle
     */
    PointerVector<NodeType> GetSubTriangleNodes(
        const int Position,
        const Geometry<NodeType>& rGeom,
        std::vector<NodeType::Pointer>& rMiddleNodes
        );

    /**
     * @brief GetSubQuadrilateralNodes gets the connectivity of a sub-quadrilateral inside
     * a geometry and the new nodes
     * @param Position The index which defines the sub-quadrilateral
     * @param rGeom The original geometry
     * @param rMiddleNode The node which divides the original geometry
     * @return The nodes defining the sub-quadrilateral
     */
    PointerVector<NodeType> GetSubQuadrilateralNodes(
        const int Position,
        const Geometry<NodeType>& rGeom,
        std::vector<NodeType::Pointer>& rMiddleNodes
        );

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    //UniformRefinementUtility& operator=(UniformRefinementUtility const& rOther);

    /// Copy constructor.
    UniformRefinementUtility(UniformRefinementUtility const& rOther);


    ///@}

}; // Class UniformRefinementUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                UniformRefinementUtility& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const UniformRefinementUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_UNIFORM_REFINEMENT_UTILITY_H_INCLUDED  defined