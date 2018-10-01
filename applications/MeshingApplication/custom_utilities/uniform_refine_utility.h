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
#include "geometries/line_3d_2.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/quadrilateral_3d_4.h"


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
    typedef std::map<std::pair<IndexType, IndexType>, IndexType> NodesInEdgeMapType;
    typedef std::unordered_map<std::array<IndexType, 4>, IndexType, KeyHasherRange<std::array<IndexType, 4>>, KeyComparorRange<std::array<IndexType, 4>>> NodesInFaceMapType;

    /// Pointer definition of UniformRefineUtility
    KRATOS_CLASS_POINTER_DEFINITION(UniformRefineUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    UniformRefineUtility(ModelPart& rModelPart, int RefinementLevel);

    /// Destructor.
    virtual ~UniformRefineUtility();


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * Execute the refinement until the final refinement is reached
     */
    void Refine();


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
    int mFinalRefinementLevel;          /// The model part will be refined to this level

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
     * Execute the refinement once
     * Only the entities with level = ThisLevel are refined to ThisLevel+1
     */
    void RefineLevel(const int& ThisLevel);

    /**
     * Create a node at the middle point of an edge
     * If the middle node is created before, do nothing
     * If the middle node does not exist, create a new one and set the nodal values
     */
    typename NodeType::Pointer  CreateNodeInEdge(const EdgeType& rEdge, const int& rRefinementLevel);

    /**
     * Get the middle node on an edge and return a pointer to it
     * CreateNodeInEdge should be executed before to ensure the node existance
     */
    typename NodeType::Pointer GetNodeInEdge(const EdgeType& rEdge);

    /**
     * Create a node at the middle point of a face
     * If the middle node is created before, do nothing
     * If the middle node does not exist, create a new one
     */
    typename NodeType::Pointer  CreateNodeInFace(const FaceType& rFace, const int& rRefinementLevel);

    /**
     * Get the node inside a face
     * The four input nodes define an element face
     * TODO: If the middle node exist, returns a pointer to the existing node
     * If the middle node does not exist, create a new one and returns a pointer to it
     */
    typename NodeType::Pointer GetNodeInFace(const FaceType& rFace);

    /**
     * Calculate the nodal data
     * The destination node is assumed to be at the mid point between
     * the origin nodes
     */
    void CalculateNodalStepData(
        NodeType::Pointer pNewNode,
        const NodeType::Pointer pNode0,
        const NodeType::Pointer pNode1
        );

      /**
       * Calculate the nodal data
       * The destination node is assumed to be at the mid point among
       * the origin nodes
       */
      void CalculateNodalStepData(
          NodeType::Pointer pNewNode,
          const NodeType::Pointer pNode0,
          const NodeType::Pointer pNode1,
          const NodeType::Pointer pNode2,
          const NodeType::Pointer pNode3
          );

    /**
     * Create an element from an origin element
     * @param pOriginElement pointer to the father element
     * @ThisNodes vector containing the sub element nodes
     * @param rRefinementLevel To assign to REFINEMENT_LEVEL flag
     */
    void CreateElement(
        Element::Pointer pOriginElement,
        std::vector<NodeType::Pointer>& rThisNodes,
        const int& rRefinementLevel
        );

    /**
     * Create condition inside from an origin condition
     * @param pOriginCondition pointer to the father condition
     * @ThisNodes vector containing the sub condition nodes
     * @param rRefinementLevel To assign to REFINEMENT_LEVEL flag
     */
    void CreateCondition(
        Condition::Pointer pOriginCondition,
        std::vector<NodeType::Pointer>& rThisNodes,
        const int& rRefinementLevel
        );

    /**
     * Return the nodes defining the i-subline
     */
    std::vector<typename NodeType::Pointer> GetSubLineNodes(
        const int Position,
        const Geometry<NodeType>& rGeom,
        NodeType::Pointer& rMiddleNode
        );

    /**
     * Return the nodes defining the i-subtriangle
     */
    std::vector<typename NodeType::Pointer> GetSubTriangleNodes(
        const int Position,
        const Geometry<NodeType>& rGeom,
        std::vector<NodeType::Pointer>& rMiddleNodes
        );

    /**
     * Return the nodes defining the i-subquadrilateral
     */
    std::vector<typename NodeType::Pointer> GetSubQuadrilateralNodes(
        const int Position,
        const Geometry<NodeType>& rGeom,
        std::vector<NodeType::Pointer>& rMiddleNodes
        );

    /**
     * This method adds a node to the sub model parts  
     * specified by a tag
     * TODO: improve this function with Model
     */
    void AddNodeToSubModelParts(NodeType::Pointer pNode, IndexType Tag);

    /**
     * This method adds the nodes to the sub model parts  
     * specified by a tag
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
    //UniformRefineUtility(UniformRefineUtility const& rOther);


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
