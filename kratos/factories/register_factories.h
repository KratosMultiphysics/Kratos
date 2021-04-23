//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_REGISTER_FACTORIES_H_INCLUDED )
#define  KRATOS_REGISTER_FACTORIES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "factories/factory.h"
#include "includes/kratos_components.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/explicit_solving_strategy.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "solving_strategies/builder_and_solvers/explicit_builder.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

namespace Kratos
{
///@name Type Definitions
///@{

typedef TUblasSparseSpace<double> SparseSpaceType;
typedef TUblasDenseSpace<double> LocalSpaceType;
typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;

///@}

typedef SolvingStrategy<SparseSpaceType, LocalSpaceType, LinearSolverType> SolvingStrategyType;
typedef ExplicitSolvingStrategy<SparseSpaceType, LocalSpaceType> ExplicitSolvingStrategyType;
typedef BuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType> BuilderAndSolverType;
typedef ExplicitBuilder<SparseSpaceType, LocalSpaceType> ExplicitBuilderType;
typedef Scheme<SparseSpaceType,LocalSpaceType> SchemeType;
typedef ConvergenceCriteria<SparseSpaceType,LocalSpaceType> ConvergenceCriteriaType;

KRATOS_API_EXTERN template class KRATOS_API(KRATOS_CORE) KratosComponents<SolvingStrategyType>;

void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, SolvingStrategyType const& ThisComponent);
  
KRATOS_API_EXTERN template class KRATOS_API(KRATOS_CORE) KratosComponents<ExplicitSolvingStrategyType>;

void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, ExplicitSolvingStrategyType const& ThisComponent);
  
KRATOS_API_EXTERN template class KRATOS_API(KRATOS_CORE) KratosComponents<BuilderAndSolverType>;

void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, BuilderAndSolverType const& ThisComponent);
  
KRATOS_API_EXTERN template class KRATOS_API(KRATOS_CORE) KratosComponents<ExplicitBuilderType>;

void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, ExplicitBuilderType const& ThisComponent);
  
KRATOS_API_EXTERN template class KRATOS_API(KRATOS_CORE) KratosComponents<SchemeType>;
  
void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, SchemeType const& ThisComponent);
  
KRATOS_API_EXTERN template class KRATOS_API(KRATOS_CORE) KratosComponents<ConvergenceCriteriaType>;
  
void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, ConvergenceCriteriaType const& ThisComponent);

}  // namespace Kratos.

#endif // KRATOS_REGISTER_FACTORIES_H_INCLUDED  defined
