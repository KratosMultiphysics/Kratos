//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:             BSD License
//                                       Kratos default license: kratos/license.txt
//
//  Main authors:   Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_TRILINOS_FACTORY_H_INCLUDED )
#define  KRATOS_TRILINOS_FACTORY_H_INCLUDED

// System includes

// External includes
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"

// Project includes
#include "includes/define.h"
#include "trilinos_space.h"
#include "factories/factory.h"

namespace Kratos
{
///@name Kratos Classes
///@{

void RegisterTrilinosStrategies();
void RegisterTrilinosBuilderAndSolvers();
void RegisterTrilinosSchemes();
void RegisterTrilinosConvergenceCriterias();

typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
typedef UblasSpace<double, Matrix, Vector> TrilinosLocalSpaceType;
typedef LinearSolver<TrilinosSparseSpaceType,  TrilinosLocalSpaceType> TrilinosLinearSolverType;

typedef SolvingStrategy<TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType> TrilinosSolvingStrategyType;
typedef BuilderAndSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType> TrilinosBuilderAndSolverType;
typedef Scheme<TrilinosSparseSpaceType,TrilinosLocalSpaceType> TrilinosSchemeType;
typedef ConvergenceCriteria<TrilinosSparseSpaceType,TrilinosLocalSpaceType> TrilinosConvergenceCriteriaType;

KRATOS_API_EXTERN template class KRATOS_API(TRILINOS_APPLICATION) KratosComponents<TrilinosSolvingStrategyType>;

#ifdef KRATOS_REGISTER_TRILINOS_STRATEGY
#undef KRATOS_REGISTER_TRILINOS_STRATEGY
#endif
#define KRATOS_REGISTER_TRILINOS_STRATEGY(name, reference) \
    KratosComponents<TrilinosSolvingStrategyType>::Add(name, reference);

KRATOS_API_EXTERN template class KRATOS_API(TRILINOS_APPLICATION) KratosComponents<TrilinosBuilderAndSolverType>;

#ifdef KRATOS_REGISTER_TRILINOS_BUILDER_AND_SOLVER
#undef KRATOS_REGISTER_TRILINOS_BUILDER_AND_SOLVER
#endif
#define KRATOS_REGISTER_TRILINOS_BUILDER_AND_SOLVER(name, reference) \
    KratosComponents<TrilinosBuilderAndSolverType>::Add(name, reference);

KRATOS_API_EXTERN template class KRATOS_API(TRILINOS_APPLICATION) KratosComponents<TrilinosSchemeType>;

#ifdef KRATOS_REGISTER_TRILINOS_SCHEME
#undef KRATOS_REGISTER_TRILINOS_SCHEME
#endif
#define KRATOS_REGISTER_TRILINOS_SCHEME(name, reference) \
    KratosComponents<TrilinosSchemeType>::Add(name, reference);

KRATOS_API_EXTERN template class KRATOS_API(TRILINOS_APPLICATION) KratosComponents<TrilinosConvergenceCriteriaType>;

#ifdef KRATOS_REGISTER_TRILINOS_CONVERGENCE_CRITERIA
#undef KRATOS_REGISTER_TRILINOS_CONVERGENCE_CRITERIA
#endif
#define KRATOS_REGISTER_TRILINOS_CONVERGENCE_CRITERIA(name, reference) \
    KratosComponents<TrilinosConvergenceCriteriaType>::Add(name, reference);

///@}

}  // namespace Kratos.

#endif // KRATOS_TRILINOS_FACTORY_H_INCLUDED  defined
