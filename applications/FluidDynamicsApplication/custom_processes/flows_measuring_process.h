//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Angel Celigueta
//
//

#ifndef KRATOS_FLOWS_MEASURING_PROCESS_H
#define KRATOS_FLOWS_MEASURING_PROCESS_H

// System includes

// External includes

// Project includes
#include "processes/process.h"
#include "includes/define.h"
#include "includes/model_part.h"

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

class FlowsMeasuringProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FlowsMeasuringProcess
    KRATOS_CLASS_POINTER_DEFINITION(FlowsMeasuringProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    FlowsMeasuringProcess(Model& rModel, Kratos::Parameters::Pointer pParameters);

    /// Destructor.
    ~FlowsMeasuringProcess() override{}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

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
        buffer << "FlowsMeasuringProcess" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "FlowsMeasuringProcess";}

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
    Model& mrModel;
    ModelPart* mpModelPartContainingTime;
    std::vector<ModelPart*> mListOfSubmodelparts;
    std::string mFilename = "flow_data_output.txt";

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
    FlowsMeasuringProcess() = delete;

    /// Assignment operator.
    FlowsMeasuringProcess& operator=(FlowsMeasuringProcess const& rOther) = delete;

    /// Copy constructor.
    FlowsMeasuringProcess(FlowsMeasuringProcess const& rOther) = delete;


    ///@}

}; // Class FlowsMeasuringProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

};  // namespace Kratos.

#endif // KRATOS_FLOWS_MEASURING_PROCESS_H
