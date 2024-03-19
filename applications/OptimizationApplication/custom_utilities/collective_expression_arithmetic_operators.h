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
//

#pragma once

// Project includes
#include "includes/define.h"

namespace Kratos {


/// @name Arithmetic Operators
/// @{

class CollectiveExpression;

KRATOS_API(OPTIMIZATION_APPLICATION) CollectiveExpression operator+(const CollectiveExpression& rLeft, const double Right);

KRATOS_API(OPTIMIZATION_APPLICATION) CollectiveExpression operator+(const double Left, const CollectiveExpression& rRight);

KRATOS_API(OPTIMIZATION_APPLICATION) CollectiveExpression operator+(const CollectiveExpression& rLeft, const CollectiveExpression& rRight);

KRATOS_API(OPTIMIZATION_APPLICATION) CollectiveExpression operator-(const CollectiveExpression& rLeft, const double Right);

KRATOS_API(OPTIMIZATION_APPLICATION) CollectiveExpression operator-(const double Left, const CollectiveExpression& rRight);

KRATOS_API(OPTIMIZATION_APPLICATION) CollectiveExpression operator-(const CollectiveExpression& rLeft, const CollectiveExpression& rRight);

KRATOS_API(OPTIMIZATION_APPLICATION) CollectiveExpression operator*(const CollectiveExpression& rLeft, const double Right);

KRATOS_API(OPTIMIZATION_APPLICATION) CollectiveExpression operator*(const double Left, const CollectiveExpression& rRight);

KRATOS_API(OPTIMIZATION_APPLICATION) CollectiveExpression operator*(const CollectiveExpression& rLeft, const CollectiveExpression& rRight);

KRATOS_API(OPTIMIZATION_APPLICATION) CollectiveExpression operator/(const CollectiveExpression& rLeft, const double Right);

KRATOS_API(OPTIMIZATION_APPLICATION) CollectiveExpression operator/(const double Left, const CollectiveExpression& rRight);

KRATOS_API(OPTIMIZATION_APPLICATION) CollectiveExpression operator/(const CollectiveExpression& rLeft, const CollectiveExpression& rRight);

/// @}


} // namespace Kratos
