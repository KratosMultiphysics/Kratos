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

#if !defined(KRATOS_CALCULATE_SENSITIVITY_METRIC_PROCESS_H_INCLUDED)
#define KRATOS_CALCULATE_SENSITIVITY_METRIC_PROCESS_H_INCLUDED

// System includes
#include <string>
#include <vector>
#include <tuple>

// External includes

// Project includes
#include "containers/model.h"
#include "containers/variable.h"
#include "includes/model_part.h"
#include "processes/process.h"

// Application incldues

namespace Kratos
{
///@addtogroup RANSApplication
///@{

///@name Kratos Classes
///@{

/**
 * @brief Apply a specific flag for nodes and conditions
 *
 * This process apply a given flag to nodes of the modelpart.
 * Then, if preferred, applies to all conditions in the given model, which has
 * the given flag applied to all the nodes in the specific condition.
 *
 */

class KRATOS_API(FLUID_DYNAMICS_APPLICATION) CalculateSensitivityMetricProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    using NodeType = ModelPart::NodeType;

    using Array3D = array_1d<double, 3>;

    /// Pointer definition of CalculateSensitivityMetricProcess
    KRATOS_CLASS_POINTER_DEFINITION(CalculateSensitivityMetricProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    CalculateSensitivityMetricProcess(
        Model& rModel,
        Parameters rParameters);

    /// Destructor.
    ~CalculateSensitivityMetricProcess() override = default;

    /// Assignment operator.
    CalculateSensitivityMetricProcess& operator=(CalculateSensitivityMetricProcess const& rOther) = delete;

    /// Copy constructor.
    CalculateSensitivityMetricProcess(CalculateSensitivityMetricProcess const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    void Execute() override;

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
    ///@name Private Member Variables
    ///@{

    Model& mrModel;
    int mEchoLevel;

    std::string mSensitivityVariableName;
    std::string mModelPartName;

    bool mIsHistoricalVariable;

    ///@}

}; // Class CalculateSensitivityMetricProcess

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const CalculateSensitivityMetricProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_CALCULATE_SENSITIVITY_METRIC_PROCESS_H_INCLUDED defined
