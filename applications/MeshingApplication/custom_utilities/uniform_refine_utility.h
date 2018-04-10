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
     * Get the node between node_a and node_b
     * The two input nodes define an element edge
     * If the middle node exist, returns a pointer to the existing node
     * If the middle node does not exist, create a new one and returns a pointer to it
     */
    Node<3>::Pointer GetNodeBetween(const Node<3>::Pointer pNode0, const Node<3>::Pointer pNode1, const int& rRefinementLevel);

    /**
     * Calculate the nodal data
     * The destination node is assumed to be at the mid point between 
     * the origin nodes
     */
    void CalculateNodalStepData(Node<3>::Pointer pNewNode, const Node<3>::Pointer pNode0, const Node<3>::Pointer pNode1);


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
    UniformRefineUtility& operator=(UniformRefineUtility const& rOther);

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
