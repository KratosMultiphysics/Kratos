//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

#if !defined(KRATOS_STANDARD_BUILDER_AND_SOLVER_FACTORY_H_INCLUDED )
#define  KRATOS_STANDARD_BUILDER_AND_SOLVER_FACTORY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "factories/builder_and_solver_factory.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class StandardBuilderAndSolverFactory
 * @ingroup KratosCore
 * @brief Here we add the functions needed for the registration of convergence criterias
 * @details Defines the standard convergence criteria factory
 * @author Vicente Mataix Ferrandiz
 * @tparam TSparseSpace The sparse space definition
 * @tparam TLocalSpace The dense space definition
 * @tparam TLinearSolver The linear solver type considered
 * @tparam TBuilderAndSolverType The convergence criteria type
 */
template <typename TSparseSpace, typename TLocalSpace, class TLinearSolver, typename TBuilderAndSolverType>
class StandardBuilderAndSolverFactory
    : public BuilderAndSolverFactory<TSparseSpace,TLocalSpace, TLinearSolver>
{
    ///@name Type Definitions
    ///@{

    // The definition of the linear solver type
    typedef LinearSolver<TSparseSpace,TLocalSpace> LinearSolverType;

    /// The definition of the convergence criteria
    typedef BuilderAndSolver<TSparseSpace,TLocalSpace, TLinearSolver> BuilderAndSolverType;

    ///@}
protected:

    ///@name Operations
    ///@{

    /**
     * @brief This method is an auxiliar method to create a new convergence criteria
     * @return The pointer to the convergence criteria of interest
     */
    typename BuilderAndSolverType::Pointer CreateBuilderAndSolver(typename LinearSolverType::Pointer pLinearSolver, Kratos::Parameters Settings) const override
    {
        return typename BuilderAndSolverType::Pointer(new TBuilderAndSolverType(pLinearSolver, Settings));
    }

    ///@}
};

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
template <typename TSparseSpace, typename TLocalSpace, class TLinearSolver, typename TBuilderAndSolverType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const StandardBuilderAndSolverFactory<TSparseSpace,TLocalSpace, TLinearSolver, TBuilderAndSolverType>& rThis)
{
    rOStream << "StandardBuilderAndSolverFactory" << std::endl;
    return rOStream;
}

///@}
///@name Input and output

void RegisterBuilderAndSolvers();

///@}

}  // namespace Kratos.

#endif // KRATOS_STANDARD_BUILDER_AND_SOLVER_FACTORY_H_INCLUDED  defined
