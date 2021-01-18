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

#if !defined(KRATOS_RANS_NUT_K_EPSILON_HIGH_RE_UPDATE_PROCESS_H_INCLUDED)
#define KRATOS_RANS_NUT_K_EPSILON_HIGH_RE_UPDATE_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "containers/model.h"

// Application includes
#include "rans_formulation_process.h"

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

class KRATOS_API(RANS_APPLICATION) RansNutKEpsilonUpdateProcess
: public RansFormulationProcess
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RansNutKEpsilonUpdateProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansNutKEpsilonUpdateProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor

    RansNutKEpsilonUpdateProcess(
        Model& rModel,
        Parameters rParameters);

    RansNutKEpsilonUpdateProcess(
        Model& rModel,
        const std::string& rModelPartName,
        const double Cmu,
        const double MinValue,
        const int EchoLevel);

    /// Destructor.
    ~RansNutKEpsilonUpdateProcess() override = default;

    /// Assignment operator.
    RansNutKEpsilonUpdateProcess& operator=(RansNutKEpsilonUpdateProcess const& rOther) = delete;

    /// Copy constructor.
    RansNutKEpsilonUpdateProcess(RansNutKEpsilonUpdateProcess const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    int Check() override;

    void ExecuteInitializeSolutionStep() override;

    void ExecuteAfterCouplingSolveStep() override;

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
    double mCmu;
    double mMinValue;
    int mEchoLevel;
    bool mIsInitialized = false;

    ///@}

}; // Class RansNutKEpsilonUpdateProcess

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const RansNutKEpsilonUpdateProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_NUT_K_EPSILON_HIGH_RE_UPDATE_PROCESS_H_INCLUDED defined
