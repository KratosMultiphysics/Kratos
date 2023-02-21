// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Armin Geiser, https://github.com/armingeiser
//
// ==============================================================================

#ifndef DIRECTION_DAMPING_UTILITIES_H
#define DIRECTION_DAMPING_UTILITIES_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/model_part.h"
#include "spatial_containers/spatial_containers.h"
#include "custom_utilities/filter_function.h"

// ==============================================================================

namespace Kratos
{

///@name Kratos Classes
///@{

/// Class for direction damping of the shape update.
/** Class for damping a specified direction from the shape update (and
 * according to chain rule of differentiation also from sensitivities.)
 * Can be used to constrain certain nodes of the design surface to slide on a plane.

*/

class KRATOS_API(SHAPE_OPTIMIZATION_APPLICATION) DirectionDampingUtilities
{
public:
    ///@name Type Definitions
    ///@{

    // Type definitions for better reading later
    typedef array_1d<double,3> array_3d;
    typedef ModelPart::NodeType NodeType;
    typedef ModelPart::NodeType::Pointer NodeTypePointer;
    typedef std::vector<NodeTypePointer> NodeVector;
    typedef std::vector<NodeTypePointer>::iterator NodeVectorIterator;
    typedef std::vector<double> DoubleVector;
    typedef std::vector<double>::iterator DoubleVectorIterator;

    // Type definitions for tree-search
    typedef Bucket< 3, NodeType, NodeVector, NodeTypePointer, NodeVectorIterator, DoubleVectorIterator > BucketType;
    typedef Tree< KDTreePartition<BucketType> > KDTree;

    /// Pointer definition of DirectionDampingUtilities
    KRATOS_CLASS_POINTER_DEFINITION(DirectionDampingUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DirectionDampingUtilities( ModelPart& modelPartToDamp, Parameters DampingSettings );

    /// Destructor.
    virtual ~DirectionDampingUtilities()
    {
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void CreateListOfNodesOfModelPart();

    void CreateSearchTreeWithAllNodesOfModelPart();

    void InitalizeDampingFactorsToHaveNoInfluence();

    void SetDampingFactors();

    FilterFunction::Pointer CreateDampingFunction( std::string damping_type ) const;

    void ThrowWarningIfNodeNeighborsExceedLimit( const ModelPart::NodeType& given_node, const unsigned int number_of_neighbors ) const;

    void DampNodalVariable( const Variable<array_3d> &rNodalVariable );


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "DirectionDampingUtilities";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "DirectionDampingUtilities";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }
    ///@}

private:
    ///@name Member Variables
    ///@{

    // Initialized by class constructor
    ModelPart& mrModelPartToDamp;
    Parameters mDampingSettings;
    array_3d mDirection;
    std::vector<double> mDampingFactors;

    // Variables for spatial search
    unsigned int mBucketSize = 100;
    unsigned int mMaxNeighborNodes = 10000;
    NodeVector mListOfNodesOfModelPart;
    KDTree::Pointer mpSearchTree;

    ///@}

}; // Class DirectionDampingUtilities

///@}

}  // namespace Kratos.

#endif // DIRECTION_DAMPING_UTILITIES_H
