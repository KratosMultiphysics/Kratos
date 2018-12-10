//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Mengjie Zhao
//
//

#ifndef KRATOS_TIME_AVERAGING_PROCESS_H
#define KRATOS_TIME_AVERAGING_PROCESS_H

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/kratos_parameters.h"

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

/// Utility to modify the distances of an embedded object in order to avoid bad intersections
/// Besides, it also deactivate the full negative distance elements
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) TimeAveragingProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TimeAveragingProcess
    KRATOS_CLASS_POINTER_DEFINITION(TimeAveragingProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
	TimeAveragingProcess(ModelPart& rModelPart);

    /// Destructor.
    ~TimeAveragingProcess() override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    void ExecuteInitialize() override;

    void ExecuteFinalizeSolutionStep() override;

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
        buffer << "TimeAveragingProcess" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "TimeAveragingProcess";}

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {}


    ///@}
    ///@name Friends
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ModelPart&                                       mrModelPart;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void AverageVelocity(ModelPart::NodeIterator& inode);
    void AveragePressure(ModelPart::NodeIterator& inode);

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
    TimeAveragingProcess() = delete;

    /// Assignment operator.
    TimeAveragingProcess& operator=(TimeAveragingProcess const& rOther) = delete;

    /// Copy constructor.
    TimeAveragingProcess(TimeAveragingProcess const& rOther) = delete;

    ///@}

}; // Class TimeAveragingProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

};  // namespace Kratos.

#endif // KRATOS_TIME_AVERAGING_PROCESS_H
