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

#pragma once

// System includes

// External includes

// Project includes
#include "factories/factory.h"
#include "includes/kratos_components.h"
#include "solving_strategies/builder_and_solvers/explicit_builder.h"
#include "spaces/default_spaces.h"

namespace Kratos
{
///@name Type Definitions
///@{

typedef DefaultSparseSpaceType SparseSpaceType;
typedef DefaultLocalSpaceType LocalSpaceType;

///@}

typedef ExplicitBuilder<SparseSpaceType, LocalSpaceType> ExplicitBuilderType;

KRATOS_API_EXTERN template class KRATOS_API(KRATOS_CORE) KratosComponents<ExplicitBuilderType>;

void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, ExplicitBuilderType const& ThisComponent);

}  // namespace Kratos.
