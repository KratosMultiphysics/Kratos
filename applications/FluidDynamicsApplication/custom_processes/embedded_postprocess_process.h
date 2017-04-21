//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

#ifndef KRATOS_DISTANCE_MODIFICATION_PROCESS_H
#define KRATOS_DISTANCE_MODIFICATION_PROCESS_H

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/cfd_variables.h"
#include "processes/find_nodal_h_process.h"

// Application includes


namespace Kratos
{
///@addtogroup FluidDynamicsApplication
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

/// Utility to filter the embedded velocity and pressure values
class EmbeddedPostprocessProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of EmbeddedPostprocessProcess
    KRATOS_CLASS_POINTER_DEFINITION(EmbeddedPostprocessProcess);

    typedef Node<3>                     NodeType;
    typedef Geometry<NodeType>      GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    EmbeddedPostprocessProcess(ModelPart& rModelPart)
    {
        mrModelPart = rModelPart;
    }

    /// Destructor.
    virtual ~EmbeddedPostprocessProcess(){}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    void ExecuteFinalizeSolutionStep() override
    {
        const array_1d<double, 3> aux_zero = ZeroVector(3);

        // Store the positive distance nodes VELOCITY in EMBEDDED_WET_VELOCITY variable. The negative distance
        // nodes EMBEDDED_WET_VELOCITY is set to zero for visualization purposes
        // The same is done for the PRESSURE using the EMBEDDED_WET_PRESSURE variable
        for (ModelPart::NodeIterator itNode=mrModelPart.NodesBegin(); itNode!=mrModelPart.NodesEnd(); ++itNode)
        {
            const double dist = itNode->FastGetSolutionStepValue(DISTANCE);
            double& emb_wet_pres = itNode->FastGetSolutionStepValue(EMBEDDED_WET_PRESSURE);
            array_1d<double, 3>& emb_wet_vel = itNode->FastGetSolutionStepValue(EMBEDDED_WET_VELOCITY);

            if (dist >= 0.0)
            {
                emb_wet_pres = itNode->FastGetSolutionStepValue(PRESSURE);
                emb_wet_vel = itNode->FastGetSolutionStepValue(VELOCITY);
            }
            else
            {
                emb_wet_pres = 0.0;
                emb_wet_vel = aux_zero;
            }
        }
    }

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "EmbeddedPostprocessProcess" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override {rOStream << "EmbeddedPostprocessProcess";}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override {}


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

    ModelPart                                  mrModelPart;

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

    /// Default constructor.
    EmbeddedPostprocessProcess(){}

    /// Assignment operator.
    EmbeddedPostprocessProcess& operator=(EmbeddedPostprocessProcess const& rOther){return *this;}

    /// Copy constructor.
    EmbeddedPostprocessProcess(EmbeddedPostprocessProcess const& rOther){}


    ///@}

}; // Class EmbeddedPostprocessProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

};  // namespace Kratos.

#endif // KRATOS_DISTANCE_MODIFICATION_PROCESS_H
