//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

#pragma once

// Core includes
#include "includes/define.h"
#include "processes/process.h"


namespace Kratos {


/** @brief Utility process that generates elements from @ref MasterSlaveConstraint "multipoint constraints" for visualization purposes.
 *  @details Default parameters:
 *           @code
 *           {
 *              "model_part_name" : "",
 *              "interval" : [0.0, "End"]
 *           }
 *           @endcode
 *  @note This process requires exclusive access to its @ref ModelPart, which it will manage while active.
 *        Elements are added and deleted to reflect active @ref MasterSlaveConstraint "multipoint constraints".
 *  @ingroup KratosCore
 */
class KRATOS_API(KRATOS_CORE) MultipointConstraintToElementProcess : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(MultipointConstraintToElementProcess);

    MultipointConstraintToElementProcess() noexcept;

    MultipointConstraintToElementProcess(Model& rModel, Parameters Settings);

    ~MultipointConstraintToElementProcess() override;

    void Execute() override;

    void ExecuteInitialize() override;

    void ExecuteInitializeSolutionStep() override;

    const Parameters GetDefaultParameters() const override;

private:
    struct Impl;
    std::unique_ptr<Impl> mpImpl;
}; // class MultipointConstraintToElementProcess


} // namespace Kratos
