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
#include "expression/traits.h"
#include "includes/define.h"

namespace Kratos {


/// @name Arithmetic Operators
/// @{

template<class TContainerType, MeshType TMeshType>
class ContainerExpression;

template<class TContainerType, MeshType TMeshType>
KRATOS_API(KRATOS_CORE) ContainerExpression<TContainerType, TMeshType> operator+(const ContainerExpression<TContainerType, TMeshType>& rLeft, const double Right);

template<class TContainerType, MeshType TMeshType>
KRATOS_API(KRATOS_CORE) ContainerExpression<TContainerType, TMeshType> operator+(const double Left, const ContainerExpression<TContainerType, TMeshType>& rRight);

template<class TContainerType, MeshType TMeshType>
KRATOS_API(KRATOS_CORE) ContainerExpression<TContainerType, TMeshType> operator+(const ContainerExpression<TContainerType, TMeshType>& rLeft, const ContainerExpression<TContainerType, TMeshType>& rRight);

template<class TContainerType, MeshType TMeshType>
KRATOS_API(KRATOS_CORE) ContainerExpression<TContainerType, TMeshType> operator-(const ContainerExpression<TContainerType, TMeshType>& rLeft, const double Right);

template<class TContainerType, MeshType TMeshType>
KRATOS_API(KRATOS_CORE) ContainerExpression<TContainerType, TMeshType> operator-(const double Left, const ContainerExpression<TContainerType, TMeshType>& rRight);

template<class TContainerType, MeshType TMeshType>
KRATOS_API(KRATOS_CORE) ContainerExpression<TContainerType, TMeshType> operator-(const ContainerExpression<TContainerType, TMeshType>& rLeft, const ContainerExpression<TContainerType, TMeshType>& rRight);

template<class TContainerType, MeshType TMeshType>
KRATOS_API(KRATOS_CORE) ContainerExpression<TContainerType, TMeshType> operator*(const ContainerExpression<TContainerType, TMeshType>& rLeft, const double Right);

template<class TContainerType, MeshType TMeshType>
KRATOS_API(KRATOS_CORE) ContainerExpression<TContainerType, TMeshType> operator*(const double Left, const ContainerExpression<TContainerType, TMeshType>& rRight);

template<class TContainerType, MeshType TMeshType>
KRATOS_API(KRATOS_CORE) ContainerExpression<TContainerType, TMeshType> operator*(const ContainerExpression<TContainerType, TMeshType>& rLeft, const ContainerExpression<TContainerType, TMeshType>& rRight);

template<class TContainerType, MeshType TMeshType>
KRATOS_API(KRATOS_CORE) ContainerExpression<TContainerType, TMeshType> operator/(const ContainerExpression<TContainerType, TMeshType>& rLeft, const double Right);

template<class TContainerType, MeshType TMeshType>
KRATOS_API(KRATOS_CORE) ContainerExpression<TContainerType, TMeshType> operator/(const double Left, const ContainerExpression<TContainerType, TMeshType>& rRight);

template<class TContainerType, MeshType TMeshType>
KRATOS_API(KRATOS_CORE) ContainerExpression<TContainerType, TMeshType> operator/(const ContainerExpression<TContainerType, TMeshType>& rLeft, const ContainerExpression<TContainerType, TMeshType>& rRight);

/// @}


} // namespace Kratos
