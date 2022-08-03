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

#if !defined(KRATOS_COMPUTE_Y_PLUS_PROCESS_H_INCLUDED)
#define KRATOS_COMPUTE_Y_PLUS_PROCESS_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/model_part.h"
#include "includes/condition.h"
#include "processes/process.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Classes
///@{

class KRATOS_API(FLUID_DYNAMICS_APPLICATION) ComputeYPlusProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    using NodeType = ModelPart::NodeType;

    using ConditionType = ModelPart::ConditionType;

    using ConditionsContainerType = ModelPart::ConditionsContainerType;

    /// Pointer definition of ComputeYPlusProcess
    KRATOS_CLASS_POINTER_DEFINITION(ComputeYPlusProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    ComputeYPlusProcess(
        Model& rModel,
        Parameters rParameters);

    /// Destructor.
    ~ComputeYPlusProcess() override = default;

    /// Assignment operator.
    ComputeYPlusProcess& operator=(ComputeYPlusProcess const& rOther) = delete;

    /// Copy constructor.
    ComputeYPlusProcess(ComputeYPlusProcess const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    int Check() override;

    void ExecuteInitializeSolutionStep() override;

    void ExecuteFinalizeSolutionStep() override;

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

    std::string mModelPartName;
    std::string mOutputVariableName;
    bool mIsOutputStoredInElements;
    bool mIsCalculatedEveryTimeStep;
    int mEchoLevel;

    bool mIsNormalsCalculated;

    ///@}

}; // Class ComputeYPlusProcess

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const ComputeYPlusProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_COMPUTE_Y_PLUS_PROCESS_H_INCLUDED defined
