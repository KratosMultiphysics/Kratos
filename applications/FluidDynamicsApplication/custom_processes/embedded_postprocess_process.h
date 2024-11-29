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

#ifndef KRATOS_EMBEDDED_POSTPROCESS_PROCESS_H
#define KRATOS_EMBEDDED_POSTPROCESS_PROCESS_H

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "processes/process.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/cfd_variables.h"
#include "utilities/openmp_utils.h"

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

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    EmbeddedPostprocessProcess(ModelPart& rModelPart) : mrModelPart(rModelPart)
    {

    }

    /// Destructor.
    ~EmbeddedPostprocessProcess() override{}

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
        ModelPart::NodesContainerType& rNodes = mrModelPart.Nodes();

        // Simple check
        if( mrModelPart.NodesBegin()->SolutionStepsDataHas( DISTANCE ) == false )
            KRATOS_ERROR << "Nodes do not have DISTANCE variable!";
        if( mrModelPart.NodesBegin()->SolutionStepsDataHas( PRESSURE ) == false )
            KRATOS_ERROR << "Nodes do not have PRESSURE variable!";
        if( mrModelPart.NodesBegin()->SolutionStepsDataHas( VELOCITY ) == false )
            KRATOS_ERROR << "Nodes do not have VELOCITY variable!";
        if( mrModelPart.NodesBegin()->SolutionStepsDataHas( EMBEDDED_WET_PRESSURE ) == false )
            KRATOS_ERROR << "Nodes do not have EMBEDDED_WET_PRESSURE variable!";
        if( mrModelPart.NodesBegin()->SolutionStepsDataHas( EMBEDDED_WET_VELOCITY ) == false )
            KRATOS_ERROR << "Nodes do not have EMBEDDED_WET_VELOCITY variable!";

        // Embedded postprocess variables set
        #pragma omp parallel for
        for (int k = 0; k < static_cast<int>(rNodes.size()); ++k)
        {
            ModelPart::NodesContainerType::iterator itNode = rNodes.begin() + k;
            const double dist = itNode->FastGetSolutionStepValue(DISTANCE);
            double& emb_wet_pres = itNode->FastGetSolutionStepValue(EMBEDDED_WET_PRESSURE);
            array_1d<double, 3>& emb_wet_vel = itNode->FastGetSolutionStepValue(EMBEDDED_WET_VELOCITY);

            if (dist >= 0.0)
            {
                emb_wet_pres = itNode->FastGetSolutionStepValue(PRESSURE);      // Store the positive distance nodes PRESSURE in EMBEDDED_WET_PRESSURE variable
                emb_wet_vel = itNode->FastGetSolutionStepValue(VELOCITY);       // Store the positive distance nodes VELOCITY in EMBEDDED_WET_VELOCITY variable
            }
            else
            {
                emb_wet_pres = 0.0;         // The negative distance nodes EMBEDDED_WET_PRESSURE is set to zero for visualization purposes
                emb_wet_vel = aux_zero;     // The negative distance nodes EMBEDDED_WET_VELOCITY is set to zero for visualization purposes
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
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "EmbeddedPostprocessProcess" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "EmbeddedPostprocessProcess";}

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {}


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

    ModelPart&                                  mrModelPart;

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
    EmbeddedPostprocessProcess() = delete;

    /// Assignment operator.
    EmbeddedPostprocessProcess& operator=(EmbeddedPostprocessProcess const& rOther) = delete;

    /// Copy constructor.
    EmbeddedPostprocessProcess(EmbeddedPostprocessProcess const& rOther) = delete;


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

#endif // KRATOS_EMBEDDED_POSTPROCESS_PROCESS_H
