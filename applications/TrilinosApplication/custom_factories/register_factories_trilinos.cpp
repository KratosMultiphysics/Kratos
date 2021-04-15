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

// System includes

// External includes

// Project includes
#include "custom_factories/register_factories_trilinos.h"

namespace Kratos
{

void AddKratosComponent(std::string const& Name, TrilinosSolvingStrategyType const& ThisComponent)
{
    KratosComponents<TrilinosSolvingStrategyType>::Add(Name, ThisComponent);
}

void AddKratosComponent(std::string const& Name, TrilinosBuilderAndSolverType const& ThisComponent)
{
    KratosComponents<TrilinosBuilderAndSolverType>::Add(Name, ThisComponent);
}

void AddKratosComponent(std::string const& Name, TrilinosSchemeType const& ThisComponent)
{
    KratosComponents<TrilinosSchemeType>::Add(Name, ThisComponent);
}

void AddKratosComponent(std::string const& Name, TrilinosConvergenceCriteriaType const& ThisComponent)
{
    KratosComponents<TrilinosConvergenceCriteriaType>::Add(Name, ThisComponent);
}

}  // namespace Kratos.
