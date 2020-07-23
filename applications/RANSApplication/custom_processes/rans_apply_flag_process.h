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

#if !defined(KRATOS_RANS_APPLY_FLAG_PROCESS_H_INCLUDED)
#define KRATOS_RANS_APPLY_FLAG_PROCESS_H_INCLUDED

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
 * @brief Apply a specific flag for nodes and conditions
 *
 * This process apply a given flag to nodes of the modelpart.
 * Then, if preferred, applies to all conditions in the given model, which has
 * the given flag applied to all the nodes in the specific condition.
 *
 */

class KRATOS_API(RANS_APPLICATION) RansApplyFlagProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RansApplyFlagProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansApplyFlagProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    RansApplyFlagProcess(
        Model& rModel,
        Parameters rParameters);

    /// Destructor.
    ~RansApplyFlagProcess() override = default;

    /// Assignment operator.
    RansApplyFlagProcess& operator=(RansApplyFlagProcess const& rOther) = delete;

    /// Copy constructor.
    RansApplyFlagProcess(RansApplyFlagProcess const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    int Check() override;

    void ExecuteInitialize() override;

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
    Parameters mrParameters;
    int mEchoLevel;

    std::string mModelPartName;
    std::string mFlagVariableName;

    bool mFlagVariableValue;

    std::vector<std::string> mModelPartsForConditionFlags;

    ///@}
    ///@name Private Operations
    ///@{

    void ApplyNodeFlags();

    void ApplyConditionFlags(ModelPart& rModelPart);

    ///@}

}; // Class RansApplyFlagProcess

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream, const RansApplyFlagProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_APPLY_FLAG_PROCESS_H_INCLUDED defined
