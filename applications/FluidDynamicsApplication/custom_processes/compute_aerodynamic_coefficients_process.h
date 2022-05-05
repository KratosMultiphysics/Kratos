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

#ifndef KRATOS_AERODYNAMIC_COEFFICIENTS_PROCESS_H
#define KRATOS_AERODYNAMIC_COEFFICIENTS_PROCESS_H

// System includes
#include <string>
#include <iostream>
#include <functional>


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

using NodeType = ModelPart::NodeType;
using DensityGetter = std::function<double(Properties const&, NodeType const&)>;

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
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) ComputeAerodynamicCoefficientsProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ComputeAerodynamicCoefficientsProcess
    KRATOS_CLASS_POINTER_DEFINITION(ComputeAerodynamicCoefficientsProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with Kratos parameters.
    ComputeAerodynamicCoefficientsProcess(
        ModelPart& rModelPart,
        Parameters Params);

    /// Constructor with Kratos model
    ComputeAerodynamicCoefficientsProcess(
        Model& rModel,
        Parameters Params);

    /// Destructor.
    ~ComputeAerodynamicCoefficientsProcess() override {}

    ///@}
    ///@name Operators
    ///@{

    void Execute() override;

    void ExecuteInitialize() override;

    void ExecuteFinalizeSolutionStep() override;

    ///@}
    ///@name Operations
    ///@{

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
    double mFreestreamPressure;
    DensityGetter mGetDensity;
    bool mComputeForces;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void CheckDefaultsAndProcessSettings(Parameters Params);
    static DensityGetter ProcessDensityDatabaseInput(const std::string& RequestedDatabase);

    void ComputePressureCoefficient(const Properties& rProperties, NodeType& rNode) const;
    static array_1d<double, 3> IntegrateOnCondition(const Condition& rCondition);

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

#endif // KRATOS_AERODYNAMIC_COEFFICIENTS_PROCESS_H
