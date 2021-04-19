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

// System includes

// External includes

// Project includes
#include "factories/register_factories.h"

namespace Kratos
{

void AddKratosComponent(std::string const& Name, SolvingStrategyType const& ThisComponent)
{
    KratosComponents<SolvingStrategyType>::Add(Name, ThisComponent);
}

void AddKratosComponent(std::string const& Name, ExplicitSolvingStrategyType const& ThisComponent)
{
    KratosComponents<ExplicitSolvingStrategyType>::Add(Name, ThisComponent);
}

void AddKratosComponent(std::string const& Name, BuilderAndSolverType const& ThisComponent)
{
    KratosComponents<BuilderAndSolverType>::Add(Name, ThisComponent);
}

void AddKratosComponent(std::string const& Name, ExplicitBuilderType const& ThisComponent)
{
    KratosComponents<ExplicitBuilderType>::Add(Name, ThisComponent);
}

void AddKratosComponent(std::string const& Name, SchemeType const& ThisComponent)
{
    KratosComponents<SchemeType>::Add(Name, ThisComponent);
}

void AddKratosComponent(std::string const& Name, ConvergenceCriteriaType const& ThisComponent)
{
    KratosComponents<ConvergenceCriteriaType>::Add(Name, ThisComponent);
}

}  // namespace Kratos.
