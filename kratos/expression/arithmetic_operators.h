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
#include "expression/expression.h"
#include "includes/define.h"


namespace Kratos {


/// @name Arithmetic Operators
/// @{

KRATOS_API(KRATOS_CORE) Expression::Pointer operator+(const Expression::ConstPointer& rpLeft, const double Right);

KRATOS_API(KRATOS_CORE) Expression::Pointer operator+(const double Left, const Expression::ConstPointer& rpRight);

KRATOS_API(KRATOS_CORE) Expression::Pointer operator+(const Expression::ConstPointer& rpLeft, const Expression::ConstPointer& rpRight);

KRATOS_API(KRATOS_CORE) Expression::Pointer operator-(const Expression::ConstPointer& rpLeft, const double Right);

KRATOS_API(KRATOS_CORE) Expression::Pointer operator-(const double Left, const Expression::ConstPointer& rpRight);

KRATOS_API(KRATOS_CORE) Expression::Pointer operator-(const Expression::ConstPointer& rpLeft, const Expression::ConstPointer& rpRight);

KRATOS_API(KRATOS_CORE) Expression::Pointer operator*(const Expression::ConstPointer& rpLeft, const double Right);

KRATOS_API(KRATOS_CORE) Expression::Pointer operator*(const double Left, const Expression::ConstPointer& rpRight);

KRATOS_API(KRATOS_CORE) Expression::Pointer operator*(const Expression::ConstPointer& rpLeft, const Expression::ConstPointer& rpRight);

KRATOS_API(KRATOS_CORE) Expression::Pointer operator/(const Expression::ConstPointer& rpLeft, const double Right);

KRATOS_API(KRATOS_CORE) Expression::Pointer operator/(const double Left, const Expression::ConstPointer& rpRight);

KRATOS_API(KRATOS_CORE) Expression::Pointer operator/(const Expression::ConstPointer& rpLeft, const Expression::ConstPointer& rpRight);

/// @}


} // namespace Kratos
