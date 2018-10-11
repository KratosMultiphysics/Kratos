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


#if !defined( KRATOS_UNIFORM_REFINE_UTILITY_H_INCLUDED )
#define KRATOS_UNIFORM_REFINE_UTILITY_H_INCLUDED


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
#include "geometries/quadrilateral_3d_4.h"
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
template<unsigned int TDim>
class UniformRefineUtility
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
     * Type of IDs
     */
    typedef std::size_t IndexType;

    /**
     * Map types to locate nodes in the mesh
     */
    typedef std::pair<IndexType, IndexType> NodeinEdgeKeyType;
    typedef std::map<NodeinEdgeKeyType, IndexType> NodesInEdgeMapType;
    typedef std::unordered_map<std::array<IndexType, 4>, IndexType, KeyHasherRange<std::array<IndexType, 4>>, KeyComparorRange<std::array<IndexType, 4>>> NodesInFaceMapType;

    /// Pointer definition of UniformRefineUtility
    KRATOS_CLASS_POINTER_DEFINITION(UniformRefineUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    UniformRefineUtility(ModelPart& rModelPart);

    /// Destructor.
    virtual ~UniformRefineUtility();


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

    IndexType mLastNodeId;           /// The node Id
    IndexType mLastElemId;           /// The element Id
    IndexType mLastCondId;           /// The condition Id
    IndexType mStepDataSize;         /// The size of the nodal database
    IndexType mBufferSize;           /// The buffer size
    NodeType::DofsContainerType mDofs;  /// Storage for the dof of the node

    std::unordered_map<IndexType,IndexType> mNodesColorMap;
    std::unordered_map<IndexType,IndexType> mCondColorMap;
    std::unordered_map<IndexType,IndexType> mElemColorMap;
    std::unordered_map<IndexType,std::vector<std::string>> mColors;    /// Where the sub model parts IDs are stored

    NodesInEdgeMapType mNodesMap;              /// Where the father nodes IDs are stored
    NodesInFaceMapType mNodesInFaceMap;        /// Where the father nodes IDs are stored

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
     */
    void ExecuteDivision(const int& rDivision);

    /**
     * @brief CreateNodeInEdge creates a node at the middle point of an edge
     * @detail If the middle node does not exist, creates a new one and set the nodal values. Otherwise, do nothing
     * @param rEdge The edge containing the two father nodes
     * @param rNumberOfDivisions The value to set NUMBER_OF_DIVISIONS flag
     */
    typename NodeType::Pointer CreateNodeInEdge(
        const EdgeType& rEdge,
        const int& rNumberOfDivisions,
        const NodeinEdgeKeyType& rNodeKey);

    /**
     * @brief GetNodeInEdge gets the middle node on an edge and return a pointer to it
     * @detail CreateNodeInEdge should be executed before to ensure the node existance
     * @param rEdge The edge containing the two father nodes
     */
    typename NodeType::Pointer GetNodeInEdge(const EdgeType& rEdge, const int& rNumberOfDivisions);

    /**
     * @brief CreateNodeInFace creates a node at the middle point of a face
     * @detail If the middle node does not exist, create a new one. Otherwise, do nothing
     * @param rFace The face containing the father nodes
     * @param rNumberOfDivisions The value to set NUMBER_OF_DIVISIONS flag
     */
    typename NodeType::Pointer CreateNodeInFace(const FaceType& rFace, const int& rNumberOfDivisions);

    /**
     * @brief GetNodeInFace gets the node inside a face
     * @detail If the middle node does not exist, create a new one and returns a pointer to it
     * @param rFace The face containing the father nodes 
     */
    typename NodeType::Pointer GetNodeInFace(const FaceType& rFace);

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
        Element::Pointer pOriginElement,
        PointerVector<NodeType>& rThisNodes,
        const int& rNumberOfDivisions
        );

    /**
     * @brief CreateCondition creates a condition from an origin condition
     * @param pOriginCondition pointer to the father condition
     * @ThisNodes vector containing the sub condition nodes
     * @param rNumberOfDivisions The value to set NUMBER_OF_DIVISIONS flag
     */
    void CreateCondition(
        Condition::Pointer pOriginCondition,
        PointerVector<NodeType>& rThisNodes,
        const int& rNumberOfDivisions
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
     * @brief GetSubTriangleNodes gets the connectivity of a sub-triangle inside a geometry and the new nodes
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
     * @brief GetSubQuadrilateralNodes gets the connectivity of a sub-quadrilateral inside a geometry and the new nodes
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

    /**
     * @brief AddNodeToSubModelParts adds a node to the sub model parts specified by a tag
     * TODO: improve this function with Model
     */
    void AddNodeToSubModelParts(NodeType::Pointer pNode, IndexType Tag);

    /**
     * @brief AddNodeToSubModelParts adds the nodes to the sub model parts specified by a tag
     * TODO: improve this function with Model
     */
    void AddNodesToSubModelParts(std::vector<NodeType::Pointer>& rThisNodes, IndexType Tag);

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
    //UniformRefineUtility& operator=(UniformRefineUtility const& rOther);

    /// Copy constructor.
    UniformRefineUtility(UniformRefineUtility const& rOther);


    ///@}

}; // Class UniformRefineUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<unsigned int TDim>
inline std::istream& operator >> (std::istream& rIStream,
                UniformRefineUtility<TDim>& rThis);

/// output stream function
template<unsigned int TDim>
inline std::ostream& operator << (std::ostream& rOStream,
                const UniformRefineUtility<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_UNIFORM_REFINE_UTILITY_H_INCLUDED  defined
