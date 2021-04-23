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

#if !defined(KRATOS_REGISTER_TRILINOS_FACTORY_H_INCLUDED )
#define  KRATOS_REGISTER_TRILINOS_FACTORY_H_INCLUDED

// System includes

// External includes
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"

// Project includes
#include "includes/define.h"
#include "trilinos_space.h"
#include "factories/register_factories.h"

namespace Kratos
{
///@name Kratos Classes
///@{

typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
typedef UblasSpace<double, Matrix, Vector> TrilinosLocalSpaceType;
typedef LinearSolver<TrilinosSparseSpaceType,  TrilinosLocalSpaceType> TrilinosLinearSolverType;

typedef SolvingStrategy<TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType> TrilinosSolvingStrategyType;
typedef BuilderAndSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType> TrilinosBuilderAndSolverType;
typedef Scheme<TrilinosSparseSpaceType,TrilinosLocalSpaceType> TrilinosSchemeType;
typedef ConvergenceCriteria<TrilinosSparseSpaceType,TrilinosLocalSpaceType> TrilinosConvergenceCriteriaType;

KRATOS_API_EXTERN template class KRATOS_API(TRILINOS_APPLICATION) KratosComponents<TrilinosSolvingStrategyType>;

void KRATOS_API(TRILINOS_APPLICATION) AddKratosComponent(std::string const& Name, TrilinosSolvingStrategyType const& ThisComponent);

KRATOS_API_EXTERN template class KRATOS_API(TRILINOS_APPLICATION) KratosComponents<TrilinosBuilderAndSolverType>;

void KRATOS_API(TRILINOS_APPLICATION) AddKratosComponent(std::string const& Name, TrilinosBuilderAndSolverType const& ThisComponent);

KRATOS_API_EXTERN template class KRATOS_API(TRILINOS_APPLICATION) KratosComponents<TrilinosSchemeType>;

void KRATOS_API(TRILINOS_APPLICATION) AddKratosComponent(std::string const& Name, TrilinosSchemeType const& ThisComponent);

KRATOS_API_EXTERN template class KRATOS_API(TRILINOS_APPLICATION) KratosComponents<TrilinosConvergenceCriteriaType>;

void KRATOS_API(TRILINOS_APPLICATION) AddKratosComponent(std::string const& Name, TrilinosConvergenceCriteriaType const& ThisComponent);

///@}

}  // namespace Kratos.

#endif // KRATOS_REGISTER_TRILINOS_FACTORY_H_INCLUDED  defined
