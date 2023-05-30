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

    /** @brief Clones the existing data container.
     *  @param rSource Source expression to copy from.
     *  @param rTarget Target expression to copy to.
     *  @details Clone an existing expression. This is light weight operation
     *           since this just clones the expression pointer. No underlying
     *           data is copied.
     */
    template <class TContainer>
    static void Clone(const ContainerExpression<TContainer>& rSource,
                      ContainerExpression<TContainer>& rTarget);

    /// @}
    /// @name Views
    /// @{

    /// @}
    /// @name Arithmetic Operations
    /// @{

    template <class TContainer>
    static void Pow(ContainerExpression<TContainer>& rBase, double Exponent);

    template <class TContainer>
    static void Pow(ContainerExpression<TContainer>& rBase,
                    const ContainerExpression<TContainer>& rExponent);

    /// @}
}; // struct ExpressionUtilities


} // namespace Kratos
