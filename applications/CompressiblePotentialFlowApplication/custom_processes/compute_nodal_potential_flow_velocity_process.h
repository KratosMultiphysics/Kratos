//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Marc Nunez, based on R.Rossi and V.Mataix work
//
//

#if !defined(KRATOS_COMPUTE_NODAL_POTENTIAL_FLOW_VELOCITY_PROCESS_INCLUDED)
#define  KRATOS_COMPUTE_NODAL_POTENTIAL_FLOW_VELOCITY_PROCESS_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/model_part.h"

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


class KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) ComputeNodalPotentialFlowVelocityProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ComputeNodalPotentialFlowVelocityProcess
    KRATOS_CLASS_POINTER_DEFINITION(ComputeNodalPotentialFlowVelocityProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor. (double)
    ComputeNodalPotentialFlowVelocityProcess(
        ModelPart& rModelPart
        );

    /// Destructor.
    ~ComputeNodalPotentialFlowVelocityProcess() override
    {
    }


    ///@}
    ///@name Operators
    ///@{

    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()()
    {
        Execute();
    }


    ///@}
    ///@name Operations
    ///@{

    /**
     * Execute method is used to execute the Process algorithms.
     * In this process the gradient of a scalar variable will be computed
     */
    void Execute() override;

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
    std::string Info() const override
    {
        return "ComputeNodalPotentialFlowVelocityProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ComputeNodalPotentialFlowVelocityProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
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

    ModelPart& mrModelPart;                                      /// The main model part

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    // TODO: Try to use enable_if!!!

    /**
     * This clears the gradient
     */
    void ClearGradient();

    /**
     * This gets the gradient value
     * @param rThisGeometry The geometry of the element
     * @param i The node index
     */
    array_1d<double, 3>& GetGradient(
        Element::GeometryType& rThisGeometry,
        unsigned int i
        );

    /**
     * This divides the gradient value by the nodal area
     */
    void PonderateGradient();

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
    ComputeNodalPotentialFlowVelocityProcess& operator=(ComputeNodalPotentialFlowVelocityProcess const& rOther);

    /// Copy constructor.
    //ComputeNodalPotentialFlowVelocityProcess(ComputeNodalPotentialFlowVelocityProcess const& rOther);


    ///@}

}; // Class ComputeNodalPotentialFlowVelocityProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.

#endif // KRATOS_COMPUTE_GRADIENT_PROCESS_INCLUDED  defined


