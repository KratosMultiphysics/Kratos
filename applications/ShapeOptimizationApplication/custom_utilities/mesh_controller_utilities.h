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

    // ==============================================================================
    void UpdateMeshAccordingInputVariable( const Variable<array_1d<double,3>> &rInputVariable )
    {
        KRATOS_TRY;

        for(auto & node_i: mrModelPart.Nodes())
        {
            array_1d<double,3> variable_value = node_i.FastGetSolutionStepValue(rInputVariable);
            node_i.X0() += variable_value[0];
            node_i.Y0() += variable_value[1];
            node_i.Z0() += variable_value[2];
            node_i.X() += variable_value[0];
            node_i.Y() += variable_value[1];
            node_i.Z() += variable_value[2];            
        }

        KRATOS_CATCH("");
    }

    // --------------------------------------------------------------------------
    void ResetMeshToReferenceMesh()
    {
        for(auto & node_i: mrModelPart.Nodes())
        {
            node_i.X() = node_i.X0();
            node_i.Y() = node_i.Y0();
            node_i.Z() = node_i.Z0();
        }
    }

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
