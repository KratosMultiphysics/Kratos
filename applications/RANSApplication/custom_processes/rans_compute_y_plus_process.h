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

#if !defined(KRATOS_RANS_COMPUTE_Y_PLUS_PROCESS_H_INCLUDED)
#define KRATOS_RANS_COMPUTE_Y_PLUS_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "containers/model.h"

// Application includes
#include "custom_processes/rans_point_execution_formulation_process.h"

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

class KRATOS_API(RANS_APPLICATION) RansComputeYPlusProcess : public RansPointExecutionFormulationProcess
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = RansPointExecutionFormulationProcess;

    /// Pointer definition of RansComputeYPlusProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansComputeYPlusProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    RansComputeYPlusProcess(
        Model& rModel,
        Parameters& rParameters);

    /// Destructor.
    ~RansComputeYPlusProcess() override = default;

    /// Assignment operator.
    RansComputeYPlusProcess& operator=(RansComputeYPlusProcess const& rOther) = delete;

    /// Copy constructor.
    RansComputeYPlusProcess(RansComputeYPlusProcess const& rOther) = delete;

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
    std::string mOutputVariableName;

    double mKappa;
    double mBeta;

    int mEchoLevel;

    ///@}
    ///@name Private Operations
    ///@{

    void ExecuteOperation() override;

    ///@}

}; // Class RansComputeYPlusProcess

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const RansComputeYPlusProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_COMPUTE_Y_PLUS_PROCESS_H_INCLUDED defined