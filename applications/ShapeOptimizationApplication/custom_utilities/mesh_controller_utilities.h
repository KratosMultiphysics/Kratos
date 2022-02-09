// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef MESH_CONTROLLER_UTILITIES_H
#define MESH_CONTROLLER_UTILITIES_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/kratos_flags.h"
#include "utilities/variable_utils.h"
#include "shape_optimization_application.h"

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

class MeshControllerUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MeshControllerUtilities
    KRATOS_CLASS_POINTER_DEFINITION(MeshControllerUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MeshControllerUtilities( ModelPart& modelPart )
        : mrModelPart( modelPart )
    {
    }

    /// Destructor.
    virtual ~MeshControllerUtilities()
    {
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // --------------------------------------------------------------------------
    void UpdateMeshAccordingInputVariable( const Variable<array_1d<double,3>> &rInputVariable )
    {
        for(auto & node_i: mrModelPart.Nodes())
            noalias(node_i.Coordinates()) += node_i.FastGetSolutionStepValue(rInputVariable);
    }

    // --------------------------------------------------------------------------
    void RevertMeshUpdateAccordingInputVariable( const Variable<array_1d<double,3>> &rInputVariable )
    {
        for(auto & node_i: mrModelPart.Nodes())
            noalias(node_i.Coordinates()) -= node_i.FastGetSolutionStepValue(rInputVariable);
    }

    // --------------------------------------------------------------------------
    void LogMeshChangeAccordingInputVariable( Variable<array_1d<double,3>> &rInputVariable )
    {
        for(auto & node_i: mrModelPart.Nodes())
            noalias(node_i.FastGetSolutionStepValue(MESH_CHANGE)) += node_i.FastGetSolutionStepValue(rInputVariable);
    }

    // --------------------------------------------------------------------------
    void SetMeshToReferenceMesh()
    {
        for(auto & node_i: mrModelPart.Nodes())
            noalias(node_i.Coordinates()) = node_i.GetInitialPosition();
    }

    // --------------------------------------------------------------------------
    void SetReferenceMeshToMesh()
    {
        for(auto & node_i: mrModelPart.Nodes())
            noalias(node_i.GetInitialPosition()) = node_i.Coordinates();
    }

    // --------------------------------------------------------------------------
    void SetDeformationVariablesToZero()
    {
        if(mrModelPart.GetNodalSolutionStepVariablesList().Has(DISPLACEMENT))
            VariableUtils().SetHistoricalVariableToZero(DISPLACEMENT,mrModelPart.Nodes());
        if(mrModelPart.GetNodalSolutionStepVariablesList().Has(ROTATION))
            VariableUtils().SetHistoricalVariableToZero(ROTATION,mrModelPart.Nodes());
    }

    // --------------------------------------------------------------------------
    void WriteCoordinatesToVariable( const Variable<array_1d<double,3>> &rVariable )
    {
        for(auto & node_i: mrModelPart.Nodes())
            noalias(node_i.FastGetSolutionStepValue(rVariable)) = node_i.Coordinates();
    }

    // --------------------------------------------------------------------------
    void SubtractCoordinatesFromVariable( const Variable<array_1d<double,3>> &rInputVariable,
        const Variable<array_1d<double,3>> &rDistanceVariable )
    {
        for(auto & node_i: mrModelPart.Nodes()){
            noalias(node_i.FastGetSolutionStepValue(rDistanceVariable)) =
                node_i.FastGetSolutionStepValue(rInputVariable) - node_i.Coordinates();
        }
    }

    // --------------------------------------------------------------------------
    void AddFirstVariableToSecondVariable( const Variable<array_1d<double,3>> &rFirstVariable, const Variable<array_1d<double,3>> &rSecondVariable )
    {
        // TODO this is a copy from optimization_utilities
        for (auto & node_i : mrModelPart.Nodes())
            noalias(node_i.FastGetSolutionStepValue(rSecondVariable)) += node_i.FastGetSolutionStepValue(rFirstVariable);
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
        return "MeshControllerUtilities";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MeshControllerUtilities";
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

    // Initialized by class constructor
    ModelPart& mrModelPart;

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
//      MeshControllerUtilities& operator=(MeshControllerUtilities const& rOther);

    /// Copy constructor.
//      MeshControllerUtilities(MeshControllerUtilities const& rOther);


    ///@}

}; // Class MeshControllerUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // MESH_CONTROLLER_UTILITIES_H
