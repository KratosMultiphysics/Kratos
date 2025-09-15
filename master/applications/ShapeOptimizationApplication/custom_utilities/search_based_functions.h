// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef SEARCH_BASED_FUNCTIONS
#define SEARCH_BASED_FUNCTIONS

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "includes/model_part.h"
#include "spatial_containers/spatial_containers.h"

// ==============================================================================

namespace Kratos
{

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

/// Short class definition.
/** Detail class definition.

*/

class KRATOS_API(SHAPE_OPTIMIZATION_APPLICATION) SearchBasedFunctions
{
public:
    ///@name Type Definitions
    ///@{

    // For better reading
    typedef array_1d<double,3> array_3d;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;

    // Type definitions for tree-search
    typedef ModelPart::NodeType NodeType;
    typedef ModelPart::NodeType::Pointer NodeTypePointer;
    typedef std::vector<NodeTypePointer> NodeVector;
    typedef std::vector<NodeTypePointer>::iterator NodeVectorIterator;
    typedef std::vector<double> DoubleVector;
    typedef std::vector<double>::iterator DoubleVectorIterator;
    typedef Bucket< 3, NodeType, NodeVector, NodeTypePointer, NodeVectorIterator, DoubleVectorIterator > BucketType;
    typedef Tree< KDTreePartition<BucketType> > KDTree;

    /// Pointer definition of SearchBasedFunctions
    KRATOS_CLASS_POINTER_DEFINITION(SearchBasedFunctions);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SearchBasedFunctions( ModelPart& modelPart );

    /// Destructor.
    virtual ~SearchBasedFunctions()
    {
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // --------------------------------------------------------------------------
    void FlagNodesInRadius( ModelPart::NodesContainerType& rNodeSet, const Kratos::Flags& rFlag, double filter_radius);

    // --------------------------------------------------------------------------

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
    virtual std::string Info() const
    {
        return "SearchBasedFunctions";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "SearchBasedFunctions";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }


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
    const unsigned int mBucketSize = 0;
    unsigned int mMaxNeighborNodes = 0;

    NodeVector mListOfNodesInModelPart;
    KDTree::Pointer mpSearchTree;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

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
//      SearchBasedFunctions& operator=(SearchBasedFunctions const& rOther);

    /// Copy constructor.
//      SearchBasedFunctions(SearchBasedFunctions const& rOther);


    ///@}

}; // Class SearchBasedFunctions

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // SEARCH_BASED_FUNCTIONS
