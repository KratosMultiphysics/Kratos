//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
//

#ifndef SYMMETRY_UTILITY_H
#define SYMMETRY_UTILITY_H

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

class KRATOS_API(OPTIMIZATION_APPLICATION) SymmetryUtility
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

    /// Pointer definition of SymmetryUtility
    KRATOS_CLASS_POINTER_DEFINITION(SymmetryUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SymmetryUtility( std::string Name, ModelPart& rModelPart, Parameters SymmetrySettings );

    /// Destructor.
    virtual ~SymmetryUtility()
    {
    }

 	struct PlaneSymmetryData{
		// point on the plane
        array_1d<double,3> Point;
        // normal to the plane
        array_1d<double,3> Normal;
        //Reflection Matrix
        Matrix ReflectionMatrix;
        // map
        std::vector<std::pair <NodeTypePointer,NodeTypePointer>> Map;
	};

 	struct RotationalSymmetryData{
		// point on the axis
        array_1d<double,3> Point;
		// angle
        double Angle;
		// number of rotational operations
        int NumRot;
        // axis
        array_1d<double,3> Axis;
        // pre-computed rotation matrices
        std::vector<Matrix> RotationMatrices;
        // map
        std::vector<std::pair <NodeTypePointer,NodeVector>> Map;
	};


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // ==============================================================================
    void Initialize();
    // --------------------------------------------------------------------------
    void Update();
    // --------------------------------------------------------------------------
    void ApplyOnVectorField( const Variable<array_3d> &rNodalVariable );
    // --------------------------------------------------------------------------
    void ApplyOnScalarField( const Variable<double> &rNodalVariable );
    // ==============================================================================

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
        return "SymmetryUtility";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "SymmetryUtility";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }

    ///@}


private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    // ==============================================================================
    // Initialized by class constructor
    // ==============================================================================
    std::string mUtilName;
    ModelPart& mrModelPart;
    Parameters mSymmetrySettings;
    bool mAxisSymmetry=false;
    RotationalSymmetryData mRotationalSymmetryData;
    bool mPlaneSymmetry=false;
    PlaneSymmetryData mPlaneSymmetryData;

    // ==============================================================================
    // Variables for spatial search
    // ==============================================================================
    NodeVector mListOfNodesOfModelPart;
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
    NodeTypePointer GetRotatedNode(NodeType& rNode, int RotationIndex);
    NodeTypePointer GetReflectedNode(NodeType& rNode);
    void GetRotationMatrix(double Angle, Matrix& rRotMat);


    ///@}


}; // Class SymmetryUtility

///@}


}  // namespace Kratos.

#endif // SYMMETRY_UTILITY_H
