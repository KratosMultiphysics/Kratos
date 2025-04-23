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

// Project Includes
#include "solving_strategies/builder_and_solvers/p_multigrid/constraint_assembler.hpp" // ConstraintAssembler
#include "includes/kratos_parameters.h" // Parameters


namespace Kratos {


/** @brief Construct a @ref ConstraintAssembler instance.
 *  @param Settings Parameters configuring the requested constraint assembler.
 *  @param rInstanceName Name of the constraint assembler instance to be constructed.
 *  @details This factory operates similary to @ref LinearSolverFactory.
 *           @p "method" setting controls what type of constraint assembler to construct,
 *           while the rest of the settings are passed on to the instance's constructor.
 *
 *           Options for @p "method":
 *           - @p "master_slave" @ref MasterSlaveConstraintAssembler
 *           - @p "augmented_lagrange" @ref AugmentedLagrangeConstraintAssembler
 *           - @p "none" @ref NoOpConstraintAssembler
 */
template <class TSparseSpace, class TDenseSpace>
typename ConstraintAssembler<TSparseSpace,TDenseSpace>::Pointer
ConstraintAssemblerFactory(Parameters Settings,
                           std::string&& rInstanceName);


} // namespace Kratos
