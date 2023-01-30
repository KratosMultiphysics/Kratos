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

#if !defined(KRATOS_RANS_K_TURBULENT_INTENSITY_INLET_PROCESS_H_INCLUDED)
#define KRATOS_RANS_K_TURBULENT_INTENSITY_INLET_PROCESS_H_INCLUDED

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
 * @brief Set turbulent kinetic energy value based on the given turbulent intensity
 *
 * This process sets turbulent kinetic energy of a given model part based on the
 * following equation
 *
 * \[
 *
 *     k = \frac{3}{2}\left(I||\underline{u}||\right)^2
 *
 * \]
 *
 * $k$ is the turbulent kinetic energy, $||\underline{u}||$ is the velocity magnitude,
 * $I$ is the turbulent intensity. If the velocity magnitude is zero, then $k_{min}$ is
 * assigned as the turbulent kinetic energy.
 *
 */

class KRATOS_API(RANS_APPLICATION) RansKTurbulentIntensityInletProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RansKTurbulentIntensityInletProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansKTurbulentIntensityInletProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    RansKTurbulentIntensityInletProcess(
        Model& rModel,
        Parameters rParameters);

    /// Destructor.
    ~RansKTurbulentIntensityInletProcess() override = default;

    /// Assignment operator.
    RansKTurbulentIntensityInletProcess& operator=(RansKTurbulentIntensityInletProcess const& rOther) = delete;

    /// Copy constructor.
    RansKTurbulentIntensityInletProcess(RansKTurbulentIntensityInletProcess const& rOther) = delete;

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
    double mTurbulentIntensity;
    double mMinValue;
    int mEchoLevel;
    bool mIsConstrained;

    ///@}

}; // Class RansKTurbulentIntensityInletProcess

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const RansKTurbulentIntensityInletProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_K_TURBULENT_INTENSITY_INLET_PROCESS_H_INCLUDED defined
