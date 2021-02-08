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

#if !defined(KRATOS_RANS_NUT_K_EPSILON_DERIVATIVES_PROCESS_H_INCLUDED)
#define KRATOS_RANS_NUT_K_EPSILON_DERIVATIVES_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "containers/model.h"
#include "processes/process.h"

// Application includes

namespace Kratos
{
///@addtogroup RANSApplication
///@{

///@name Kratos Classes
///@{

/**
 * @brief Calculates turbulent kinematic viscosity
 *
 * This process uses following formula to calculate turbulent kinematic viscosity
 *
 * \[
 *      \nu_t = C_\mu\frac{k^2}{\epsilon}
 * \]
 *
 * $k$ is the turbulent kinetic energy, $\epsilon$ is the turbulent energy dissipation rate
 */

class KRATOS_API(RANS_APPLICATION) RansNutKEpsilonDerivativesProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RansNutKEpsilonDerivativesProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansNutKEpsilonDerivativesProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor

    RansNutKEpsilonDerivativesProcess(
        Model& rModel,
        Parameters rParameters);

    RansNutKEpsilonDerivativesProcess(
        Model& rModel,
        const std::string& rModelPartName,
        const double MinValue,
        const int EchoLevel);

    /// Destructor.
    ~RansNutKEpsilonDerivativesProcess() override = default;

    /// Assignment operator.
    RansNutKEpsilonDerivativesProcess& operator=(RansNutKEpsilonDerivativesProcess const& rOther) = delete;

    /// Copy constructor.
    RansNutKEpsilonDerivativesProcess(RansNutKEpsilonDerivativesProcess const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    int Check() override;

    void ExecuteInitializeSolutionStep() override;

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
    double mMinValue;
    int mEchoLevel;
    bool mIsInitialized = false;

    ///@}

}; // Class RansNutKEpsilonDerivativesProcess

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const RansNutKEpsilonDerivativesProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_NUT_K_EPSILON_DERIVATIVES_PROCESS_H_INCLUDED defined
