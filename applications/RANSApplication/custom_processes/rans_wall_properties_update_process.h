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

#if !defined(KRATOS_RANS_WALL_PROPERTIES_UPDATE_PROCESS_H_INCLUDED)
#define KRATOS_RANS_WALL_PROPERTIES_UPDATE_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "containers/model.h"

// Application includes
#include "rans_point_execution_formulation_process.h"

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

class KRATOS_API(RANS_APPLICATION) RansWallPropertiesUpdateProcess : public RansPointExecutionFormulationProcess
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RansWallPropertiesUpdateProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansWallPropertiesUpdateProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    RansWallPropertiesUpdateProcess(
        Model& rModel,
        Parameters& rParameters);

    RansWallPropertiesUpdateProcess(
        Model& rModel,
        const std::string& rModelPartName,
        const bool UpdateWallNormals,
        const bool UpdateConditionWallHeights,
        const std::vector<std::string>& rUpdateExecutionPoints,
        const int EchoLevel);

    /// Destructor.
    ~RansWallPropertiesUpdateProcess() override = default;

    /// Assignment operator.
    RansWallPropertiesUpdateProcess& operator=(RansWallPropertiesUpdateProcess const& rOther) = delete;

    /// Copy constructor.
    RansWallPropertiesUpdateProcess(RansWallPropertiesUpdateProcess const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

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
    int mEchoLevel;

    bool mUpdateWallNormals;
    bool mUpdateConditionWallHeights;

    ///@}
    ///@name Private Operations
    ///@{

    void ExecuteOperation() override;

    ///@}

}; // Class RansWallPropertiesUpdateProcess

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const RansWallPropertiesUpdateProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_WALL_PROPERTIES_UPDATE_PROCESS_H_INCLUDED defined