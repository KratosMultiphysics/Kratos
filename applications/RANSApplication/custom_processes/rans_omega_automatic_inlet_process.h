//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Dharmin Shah
//                   Bence Rochlitz
//
//  Supervised by:   Jordi Cotela
//                   Suneth Warnakulasuriya
//

#if !defined(KRATOS_RANS_OMEGA_AUTOMATIC_INLET_PROCESS_H_INCLUDED)
#define KRATOS_RANS_OMEGA_AUTOMATIC_INLET_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "containers/array_1d.h"
#include "containers/model.h"
#include "processes/process.h"

namespace Kratos
{
///@addtogroup RANSModellingApplication
///@{

///@name Kratos Classes
///@{

/**
 * @brief Sets omega value best on turbulent mixing length
 *
 * This process sets omega values based on the following formula
 *
 * \[
 *
 *  \omega = \frac{k^{0.5}}{C_\mu * L}
 *
 * \]
 *
 * In here $k$ is turbulent kinetic energy, $\omega$ is turbulent specific
 * energy dissipation rate, and $L$ is turbulent mixing length.
 *
 */

class KRATOS_API(RANS_APPLICATION) RansOmegaAutomaticInletProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RansOmegaAutomaticInletProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansOmegaAutomaticInletProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    RansOmegaAutomaticInletProcess(
        Model& rModel,
        Parameters& rParameters);

    /// Destructor.
    ~RansOmegaAutomaticInletProcess() override = default;

    /// Assignment operator.
    RansOmegaAutomaticInletProcess& operator=(RansOmegaAutomaticInletProcess const& rOther) = delete;

    /// Copy constructor.
    RansOmegaAutomaticInletProcess(RansOmegaAutomaticInletProcess const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    void ExecuteInitializeSolutionStep() override;

    int Check() override;

    const Parameters GetDefaultParameters() const override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

    ///@}

private:
    ///@name Member Variables
    ///@{

    Model& mrModel;

    std::string mModelPartName;
    double mTurbulentMixingLength;
    array_1d<double, 3> mWallLocation;
    array_1d<double, 3> mWallOutwardPointintUnitNormal;
    double mKinematicViscosity;
    double mBeta;
    double mMinValue;
    int mEchoLevel;

    ///@}

}; // Class RansOmegaAutomaticInletProcess

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const RansOmegaAutomaticInletProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_OMEGA_AUTOMATIC_INLET_PROCESS_H_INCLUDED defined