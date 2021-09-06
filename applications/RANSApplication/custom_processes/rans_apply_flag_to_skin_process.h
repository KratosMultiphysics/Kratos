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

class KRATOS_API(RANS_APPLICATION) RansApplyFlagToSkinProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RansApplyFlagToSkinProcess
    KRATOS_CLASS_POINTER_DEFINITION(RansApplyFlagToSkinProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    RansApplyFlagToSkinProcess(
        Model& rModel,
        Parameters rParameters);

    /// Destructor.
    ~RansApplyFlagToSkinProcess() override = default;

    /// Assignment operator.
    RansApplyFlagToSkinProcess& operator=(RansApplyFlagToSkinProcess const& rOther) = delete;

    /// Copy constructor.
    RansApplyFlagToSkinProcess(RansApplyFlagToSkinProcess const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    void ExecuteInitialize() override;

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
    int mEchoLevel;

    std::string mModelPartName;
    std::string mFlagVariableName;

    bool mFlagVariableValue;

    std::vector<std::string> mModelPartsForConditionFlags;

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Applies flags to nodes
     *
     * This method applies given mFlagVariableValue to mFlagVariableName in all
     * nodes of mModelPartName.
     *
     */
    void ApplyNodeFlags();

    /**
     * @brief Applies flags to conditions
     *
     * This method applies mFlagVariableValue for mFlagVariableName in conditions for given rModelPart,
     * if given condition has all of its nodes with mFlagVariableValue for flag mFlagVariableName. Otherwise, it
     * applies inverse mFlagVariableValue for mFlagVariableName flag in the same condition.
     *
     * @param rModelPart    Model part to look for conditions
     */
    void ApplyConditionFlags(
        ModelPart& rModelPart);

    ///@}

}; // Class RansApplyFlagToSkinProcess

///@}
///@name Input and output
///@{

/// output stream function
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const RansApplyFlagToSkinProcess& rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_APPLY_FLAG_PROCESS_H_INCLUDED defined
