//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:   Juan I. Camarotti
//
//

#ifndef KRATOS_COMPUTE_AERODYNAMIC_COEFFICIENTS_PROCESS_H
#define KRATOS_COMPUTE_AERODYNAMIC_COEFFICIENTS_PROCESS_H

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
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

/// This process computes the lift coefficient as a function of reference fluid properties
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) ComputeAerodynamicCoefficientsProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{
    using NodeType = ModelPart::NodeType;

    /// Pointer definition of ComputeAerodynamicCoefficientsProcess
    KRATOS_CLASS_POINTER_DEFINITION(ComputeAerodynamicCoefficientsProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with Kratos model
    ComputeAerodynamicCoefficientsProcess(
        Model& rModel,
        Parameters Params);

    /// Destructor.
    ~ComputeAerodynamicCoefficientsProcess() override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    const Parameters GetDefaultParameters() const override;

    void Execute() override;

    void ExecuteBeforeOutputStep() override;

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
        std::stringstream buffer;
        buffer << "ComputeAerodynamicCoefficientsProcess" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "ComputeAerodynamicCoefficientsProcess";}

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

    ModelPart& mrModelPart;
    double mReference_Surface = 0.0;
    double mReference_Chord   = 1.0;
    double mReference_Span    = 1.0;
    double mQInf              = 1.0;
    array_1d<double,3> mMomentReferencePoint = ZeroVector(3);
    double mAngleOfAttack = 0.0;
    double mSideslipAngle = 0.0;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void ReadProcessParameters(const Parameters& Params);

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
    ComputeAerodynamicCoefficientsProcess() = delete;

    /// Assignment operator.
    ComputeAerodynamicCoefficientsProcess& operator=(ComputeAerodynamicCoefficientsProcess const& rOther) = delete;

    /// Copy constructor.
    ComputeAerodynamicCoefficientsProcess(ComputeAerodynamicCoefficientsProcess const& rOther) = delete;

    ///@}

}; // Class ComputeAerodynamicCoefficientsProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

};  // namespace Kratos.

#endif // KRATOS_COMPUTE_AERODYNAMIC_COEFFICIENTS_PROCESS_H
