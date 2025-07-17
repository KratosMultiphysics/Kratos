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

// Project includes
#include "solving_strategies/builder_and_solvers/p_multigrid/constraint_assembler_factory.hpp" // ConstraintAssemblerFactory
#include "solving_strategies/builder_and_solvers/p_multigrid/noop_constraint_assembler.hpp" // NoOpConstraintAssembler
#include "solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler.hpp" // MasterSlaveConstraintAssembler
#include "solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler.hpp" // AugmentedLagrangeConstraintAssembler
#include "spaces/ublas_space.h" // TUblasSparseSpace, TUblasDenseSpace


namespace Kratos {


template <class TSparse, class TDense>
typename ConstraintAssembler<TSparse,TDense>::Pointer
ConstraintAssemblerFactory(Parameters Settings,
                           std::string&& rInstanceName)
{
    KRATOS_TRY
    const std::string imposition_name = Settings["method"].Get<std::string>();
    if (imposition_name == "none") {
        return std::make_shared<NoOpConstraintAssembler<TSparse,TDense>>(Settings, std::move(rInstanceName));
    } else if (imposition_name == "master_slave") {
        return std::make_shared<MasterSlaveConstraintAssembler<TSparse,TDense>>(Settings, std::move(rInstanceName));
    } else if (imposition_name == "augmented_lagrange") {
        return std::make_shared<AugmentedLagrangeConstraintAssembler<TSparse,TDense>>(Settings, std::move(rInstanceName));
    } else {
        std::stringstream message;
        message << "Unsupported constraint imposition \"" << imposition_name << "\". Options are:\n"
                << "\t\"none\"\n"
                << "\t\"master_slave\"\n"
                << "\t\"augmented_lagrange\"";
        KRATOS_ERROR << message.str();
    }
    KRATOS_CATCH("")
}


#define KRATOS_INSTANTIATE_CONSTRAINT_ASSEMBLER_FACTORY(TSparse, TDense)                                                        \
    template ConstraintAssembler<TSparse,TDense>::Pointer ConstraintAssemblerFactory<TSparse,TDense>(Parameters, std::string&&)

KRATOS_INSTANTIATE_CONSTRAINT_ASSEMBLER_FACTORY(TUblasSparseSpace<double>, TUblasDenseSpace<double>);

KRATOS_INSTANTIATE_CONSTRAINT_ASSEMBLER_FACTORY(TUblasSparseSpace<float>, TUblasDenseSpace<double>);

#undef KRATOS_INSTANTIATE_CONSTRAINT_ASSEMBLER_FACTORY


} // namespace Kratos
