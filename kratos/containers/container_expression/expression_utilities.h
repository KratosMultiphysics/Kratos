//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//                   Máté Kelemen
//

#pragma once

// Project includes
#include "containers/container_expression/container_expression.h"
#include "includes/define.h"


namespace Kratos {


struct KRATOS_API(KRATOS_CORE) ExpressionUtilities
{
    /// @name Copy operations
    /// @{

    template <class TContainer>
    static void Clone(const ContainerExpression<TContainer>& rSource,
                      ContainerExpression<TContainer>& rTarget);

    /// @}
    /// @name Views
    /// @{

    /// @}
    /// @name Arithmetic Operations
    /// @{

    /// @}
}; // struct ExpressionUtilities


} // namespace Kratos
