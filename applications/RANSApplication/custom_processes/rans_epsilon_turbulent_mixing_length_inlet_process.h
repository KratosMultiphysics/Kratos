//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(KRATOS_RANS_EPSILON_TURBULENT_MIXING_LENGTH_INLET_PROCESS_H_INCLUDED)
#define KRATOS_RANS_EPSILON_TURBULENT_MIXING_LENGTH_INLET_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "containers/model.h"
#include "processes/process.h"

namespace Kratos
{
///@addtogroup RANSApplication
///@{

///@name Kratos Classes
///@{

/**
 * @brief Sets epsilon value best on turbulent mixing length
 *
 * This process sets epsilon values based on the following formula
 *
 * \[
 *
 *  \epsilon = \frac{C_\mu^{0.75} max\left\lbrace k, 0.0\right\rbrace^{1.5}}{L}
 *
 * \]
 *
 * In here $k$ is turbulent kinetic energy, $\epsilon$ is turbulent energy
 * dissipation rate, and $L$ is turbulent mixing length.
 *
 */

class KRATOS_API(RANS_APPLICATION) RansEpsilonTurbulentMixingLengthInletProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RansEpsilonTurbulentMixingLengthInletProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansEpsilonTurbulentMixingLengthInletProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    RansEpsilonTurbulentMixingLengthInletProcess(
        Model& rModel,
        Parameters rParameters);

    /// Destructor.
    ~RansEpsilonTurbulentMixingLengthInletProcess() override = default;

    /// Assignment operator.
    RansEpsilonTurbulentMixingLengthInletProcess& operator=(RansEpsilonTurbulentMixingLengthInletProcess const& rOther) = delete;

    /// Copy constructor.
    RansEpsilonTurbulentMixingLengthInletProcess(RansEpsilonTurbulentMixingLengthInletProcess const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    void ExecuteInitialize() override;

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
    double mCmu_75;
    double mMinValue;
    bool mIsConstrained;
    int mEchoLevel;

    ///@}

}; // Class RansEpsilonTurbulentMixingLengthInletProcess

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const RansEpsilonTurbulentMixingLengthInletProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_EPSILON_TURBULENT_MIXING_LENGTH_INLET_PROCESS_H_INCLUDED defined
