//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jonathan Nuttall
//

#pragma once

#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "utilities/function_parser_utility.h"

#include "../../StructuralMechanicsApplication/custom_processes/set_moving_load_process.h"

namespace Kratos
{

/**
 * @class SetMultipleMovingLoadsProcess
 * @ingroup GeoMechanicsApplication
 * @brief Process to set and manage multiple moving loads offset according to a configuration variable
 * @details This process applies multiple moving loads at intervals according to a set configuration pattern.
 * @author Jonathan Nuttall
 */

class KRATOS_API(GEO_MECHANICS_APPLICATION) SetMultipleMovingLoadsProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of SetMultipleMovingLoadsProcess
    KRATOS_CLASS_POINTER_DEFINITION(SetMultipleMovingLoadsProcess);

    using SizeType = std::size_t;

    ///@}
    ///@name Life Cycle
    ///@{

    SetMultipleMovingLoadsProcess(ModelPart& rModelPart, const Parameters& rProcessSettings);

    SetMultipleMovingLoadsProcess(const SetMultipleMovingLoadsProcess&)            = delete;
    SetMultipleMovingLoadsProcess& operator=(const SetMultipleMovingLoadsProcess&) = delete;

    ///@}

    ///@}

    ///@name Operations
    ///@{

    /**
     * \brief  Initializes the set moving load process. Check if load functions and a velocity function are present in the parameters.
     * Sort vector of conditions, and find the start position of the moving load, within the conditions vector.
     */
    void ExecuteInitialize() override;

    /**
     * \brief Initialize solution step. Calculate the load based on the load functions if present, else retrieve the load from the input parameters.
     * Loop over the conditions and find, on which condition the load is located. Then set the load on the condition element, if the load is located
     * within the element. If the moving load is not located on the condition element, set the load to zero.
     */
    void ExecuteInitializeSolutionStep() override;

    /**
     * \brief Finalizes solution step. Sets load velocity based on load velocity function if present, else load velocity is retrieved from the input values.
     * Then move the load based on the current position and the load velocity.
     */
    void ExecuteFinalizeSolutionStep() override;

    ///@}

private:
    ///@name Member Variables
    ///@{
    ModelPart&                                            mrModelPart;
    Parameters                                            mParameters;
    std::vector<Kratos::unique_ptr<SetMovingLoadProcess>> mMovingPointLoadsProcesses;
    ///@}
    ///
    ///@name Operations
    ///@{

    /**
     * \brief Clones condition into a new sub body part of the compute model part
     */
    ModelPart& CloneMovingConditionInComputeModelPart(const std::string& NewBodyPartName);

    /**
     * \brief Get maximum index of current conditions in root
     */
    [[nodiscard]] int GetMaxConditionsIndex() const;

    /**
     * \brief Remove cloned conditions as they are not executed
     */
    void RemoveClonedConditions();

    ///@}
};
} // namespace Kratos