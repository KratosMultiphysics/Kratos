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
#include "includes/define.h" // KRATOS_API, KRATOS_CLASS_POINTER_DEFINITION
#include "processes/process.h" // Process, Model, Parameters

// System includes
#include <memory> // unique_ptr


namespace Kratos {


/** @brief Utility process that generates elements from @ref MasterSlaveConstraint "multipoint constraints" for visualization purposes.
 *  @details Default parameters:
 *           @code
 *           {
 *              "input_model_part_name" : "",
 *              "output_model_part_name" : "",
 *              "interval" : [0.0, "End"]
 *           }
 *           @endcode
 *           Multipoint constraints act on pairs of @ref Dof "DoFs", not on @ref Node "nodes", so elements are further
 *           partitioned into sub model parts, based on the pair of @ref Variable "variables" they constrain. For example,
 *           if @a output_model_part_name is @a root and a multipoint constraint relates @a DISPLACEMENT_X to @a PRESSURE,
 *           then @ref Element2D2N "Element2D2Ns" will be constructed in the @a root.DISPLACEMENT_X_PRESSURE sub model part.
 *
 *  @param input_model_part_name full name of the model part to scane multipoint constraints in.
 *  @param output_model_part_name full name of the model part to generate elements in. It is created if it does not exist yet.
 *  @param interval time interval to manage the generated elements in.
 *
 *  @note The coefficient of the relationships is stored in elements' @a SCALAR_LAGRANGE_MULTIPLIER variable.
 *  @note This process requires exclusive access to its output @ref ModelPart, which it will manage while active.
 *        Elements are added and deleted to reflect active @ref MasterSlaveConstraint "multipoint constraints".
 *  @todo This process would be better off generating @ref Geometry "geometries", but no output process writes those yet @matekelemen.
 *
 *  @ingroup KratosCore
 */
class KRATOS_API(KRATOS_CORE) MultipointConstraintToElementProcess final : public Process
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
