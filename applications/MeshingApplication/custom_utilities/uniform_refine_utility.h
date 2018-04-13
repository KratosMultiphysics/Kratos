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
    
    ModelPart& mrModelPart;
    int mFinalRefinementLevel;
    std::map<std::pair<int, int>, int> mNodesMap;
    //std::unordered_map<std::pair<int, int>, int, KeyHasherRange<std::pair<int, int>>, KeyComparorRange<std::pair<int, int>> > mNodesMap;
    unsigned int mLastNodeId;               /// The node Id
    unsigned int mLastElemId;               /// The element Id
    unsigned int mLastCondId;               /// The condition Id
    unsigned int mStepDataSize;             /// The size of the database
    unsigned int mBufferSize;               /// The size of the buffer
    NodeType::DofsContainerType mDofs;      /// Storage for the dof of the node


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
     * If the middle node does not exist, create a new one
     */
    void CreateNodeInEdge(
        const EdgeType& rEdge,
        const int& rRefinementLevel
        );

    /**
     * Get the middle node on an edge and return a pointer to it
     * CreateNodeInEdge should be executed before to ensure the node existance
     */
    Node<3>::Pointer GetNodeInEdge(const EdgeType& rEdge);

    /**
     * Get the node between node_a and node_b
     * The two input nodes define an element edge
     * If the middle node exist, returns a pointer to the existing node
     * If the middle node does not exist, create a new one and returns a pointer to it
     */
    Node<3>::Pointer GetNodeBetween(
        const NodeType::Pointer pNode0,
        const NodeType::Pointer pNode1,
        const int& rRefinementLevel
        );

    /**
     * Create a node at the middle point of a face
     * If the middle node is created before, do nothing
     * If the middle node does not exist, create a new one
     */
    void CreateNodeInFace(
        const FaceType& rFace,
        const int& rRefinementLevel
        );

    /**
     * Get the node inside a face
     * The four input nodes define an element face
     * TODO: If the middle node exist, returns a pointer to the existing node
     * If the middle node does not exist, create a new one and returns a pointer to it
     */
    Node<3>::Pointer GetNodeInFace(
        const NodeType::Pointer pNode0,
        const NodeType::Pointer pNode1,
        const NodeType::Pointer pNode2,
        const NodeType::Pointer pNode3,
        const int& rRefinementLevel
        );

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
     * Create a triangle inside the origin element
     */
    void CreateElement(
        Element::Pointer pOriginElement,
        std::vector<NodeType::Pointer> ThisNodes,
        const int& rRefinementLevel
        );

    /**
     * Return the nodes defining the i-subtriangle
     */
    std::vector<Node<3>::Pointer> GetSubTriangleNodes(
        int Position,
        Geometry<NodeType>& rGeom,
        std::array<NodeType::Pointer, 3>& rMiddleNodes
        );

    /**
     * Return the nodes defining the i-subtriangle
     */
    std::vector<Node<3>::Pointer> GetSubQuadrilateralNodes(
        int Position,
        Geometry<NodeType>& rGeom,
        std::array<NodeType::Pointer, 5>& rMiddleNodes
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
