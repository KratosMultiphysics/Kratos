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
#include <algorithm>
#include <unordered_map>

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

class SearchBasedFunctions
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
    SearchBasedFunctions( ModelPart& modelPart )
        : mrModelPart( modelPart ),
          mBucketSize(100),
          mMaxNeighborNodes(10000)
    {
        // Create search tree
        mListOfNodesInModelPart.resize(mrModelPart.Nodes().size());
        int counter = 0;
        for (ModelPart::NodesContainerType::iterator node_it = mrModelPart.NodesBegin(); node_it != mrModelPart.NodesEnd(); ++node_it)
        {
            NodeTypePointer pnode = *(node_it.base());
            mListOfNodesInModelPart[counter++] = pnode;
        }

        mpSearchTree = Kratos::shared_ptr<KDTree>(new KDTree(mListOfNodesInModelPart.begin(), mListOfNodesInModelPart.end(), mBucketSize));
    }

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
    void FlagNodesInRadius( ModelPart::NodesContainerType& rNodeSet, const Kratos::Flags& rFlag, double filter_radius)
    {
        // Mark all nodes within one filter radius away from the fixed boundary
        bool is_max_number_too_small = true;

        while(is_max_number_too_small == true)
        {
            is_max_number_too_small = false;
            for(auto& node_i : rNodeSet)
            {
                NodeVector neighbor_nodes( mMaxNeighborNodes );
                DoubleVector resulting_squared_distances( mMaxNeighborNodes, 0.0 );
                unsigned int number_of_neighbors = mpSearchTree->SearchInRadius( node_i,
                                                                                 filter_radius,
                                                                                 neighbor_nodes.begin(),
                                                                                 resulting_squared_distances.begin(),
                                                                                 mMaxNeighborNodes );

                if(number_of_neighbors >= mMaxNeighborNodes)
                {
                    KRATOS_WARNING("ShapeOpt::FlagNodesInRadius") << "For node " << node_i.Id() << " maximum number of neighbor nodes (=" << mMaxNeighborNodes << " nodes) is reached! Increasing maximum number by factor 2. " << std::endl;
                    is_max_number_too_small = true;
                    break;
                }

                // Loop over all nodes in radius (including node on damping region itself)
                for(unsigned int j_itr = 0 ; j_itr<number_of_neighbors ; j_itr++)
                {
                    ModelPart::NodeType& neighbor_node = *neighbor_nodes[j_itr];
                    neighbor_node.Set(rFlag,true);
                }
            }
            mMaxNeighborNodes *=2;
        }
    }

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
