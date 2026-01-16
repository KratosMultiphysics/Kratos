//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Eduard Gómez
//
//

#ifndef KRATOS_COMPUTE_LIFT_COEFFICIENT_PROCESS_H
#define KRATOS_COMPUTE_LIFT_COEFFICIENT_PROCESS_H

// System includes
#include <string>
#include <iostream>

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
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) ComputeLiftCoefficientProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{
    using NodeType = ModelPart::NodeType;

    /// Pointer definition of ComputeLiftCoefficientProcess
    KRATOS_CLASS_POINTER_DEFINITION(ComputeLiftCoefficientProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with Kratos model
    ComputeLiftCoefficientProcess(
        Model& rModel,
        Parameters Params);

    /// Destructor.
    ~ComputeLiftCoefficientProcess() override {}

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
        buffer << "ComputeLiftCoefficientProcess" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "ComputeLiftCoefficientProcess";}

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
    double mAngle_of_Attack = 0.0;   // Angle of attack
    double mReference_Surface = 0.0;  // Reference surface
    std::function<double(const NodeType&)> mGetPressureCoefficient;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void ReadFreestreamValues(const Parameters& Params);

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
    ComputeLiftCoefficientProcess() = delete;

    /// Assignment operator.
    ComputeLiftCoefficientProcess& operator=(ComputeLiftCoefficientProcess const& rOther) = delete;

    /// Copy constructor.
    ComputeLiftCoefficientProcess(ComputeLiftCoefficientProcess const& rOther) = delete;

    ///@}

}; // Class ComputeLiftCoefficientProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

};  // namespace Kratos.

#endif // KRATOS_COMUTE_LIFT_COEFFICIENT_PROCESS_H
