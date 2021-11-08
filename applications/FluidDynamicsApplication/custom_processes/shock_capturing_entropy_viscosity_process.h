//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Eduard Gómez
//

#if !defined(KRATOS_SHOCK_CAPTURING_ENTROPY_VISCOSITY_H_INCLUDED)
#define  KRATOS_SHOCK_CAPTURING_ENTROPY_VISCOSITY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "processes/process.h"

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

class KRATOS_API(FLUID_DYNAMICS_APPLICATION) ShockCapturingEntropyViscosityProcess : public Process
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of ShockCapturingEntropyViscosityProcess
    KRATOS_CLASS_POINTER_DEFINITION(ShockCapturingEntropyViscosityProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with model
    ShockCapturingEntropyViscosityProcess(
        Model& rModel,
        Parameters rParameters)
        : ShockCapturingEntropyViscosityProcess(rModel.GetModelPart(rParameters["model_part_name"].GetString()), rParameters) {};

    /// Constructor with model part
    ShockCapturingEntropyViscosityProcess(
        ModelPart& rModelPart,
        Parameters rParameters)
        : Process()
        , mrModelPart(rModelPart)
    {
        ValidateAndAssignParameters(rParameters);
    };

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    int Check() override;

    const Parameters GetDefaultParameters() const override;

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
        buffer << "ShockCapturingEntropyViscosityProcess";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override { rOStream << "ShockCapturingEntropyViscosityProcess"; }

    /// Print object's data.
    void PrintData(std::ostream &rOStream) const override {}

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

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void ValidateAndAssignParameters(Parameters &rParameters);

    static double ComputeEntropy(const double Density, const double Pressure, const double Gamma);

    static double ComputeH(const Element& rElement);

    double ComputeElementalEntropy();

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
    ShockCapturingEntropyViscosityProcess& operator=(ShockCapturingEntropyViscosityProcess const& rOther);

    /// Copy constructor.
    ShockCapturingEntropyViscosityProcess(ShockCapturingEntropyViscosityProcess const& rOther);

    ///@}

}; // Class ShockCapturingEntropyViscosityProcess

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const ShockCapturingEntropyViscosityProcess& rThis);

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SHOCK_CAPTURING_ENTROPY_VISCOSITY_H_INCLUDED  defined
