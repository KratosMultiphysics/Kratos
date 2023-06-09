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
#include "containers/container_expression/traits.h"
#include "includes/define.h"

namespace Kratos {


/// @name Arithmetic Operators
/// @{

template<class TContainerType, MeshType TMeshType>
class ContainerExpression;

template<class TContainerType, MeshType TMeshType>
KRATOS_API(KRATOS_CORE) ContainerExpression<TContainerType, TMeshType> operator+(ContainerExpression<TContainerType, TMeshType> Left, const double Right);

template<class TContainerType, MeshType TMeshType>
KRATOS_API(KRATOS_CORE) ContainerExpression<TContainerType, TMeshType> operator+(const double Left, ContainerExpression<TContainerType, TMeshType> Right);

template<class TContainerType, MeshType TMeshType>
KRATOS_API(KRATOS_CORE) ContainerExpression<TContainerType, TMeshType> operator+(ContainerExpression<TContainerType, TMeshType> Left, ContainerExpression<TContainerType, TMeshType> Right);

template<class TContainerType, MeshType TMeshType>
KRATOS_API(KRATOS_CORE) ContainerExpression<TContainerType, TMeshType> operator-(ContainerExpression<TContainerType, TMeshType> Left, const double Right);

template<class TContainerType, MeshType TMeshType>
KRATOS_API(KRATOS_CORE) ContainerExpression<TContainerType, TMeshType> operator-(const double Left, ContainerExpression<TContainerType, TMeshType> Right);

template<class TContainerType, MeshType TMeshType>
KRATOS_API(KRATOS_CORE) ContainerExpression<TContainerType, TMeshType> operator-(ContainerExpression<TContainerType, TMeshType> Left, ContainerExpression<TContainerType, TMeshType> Right);

template<class TContainerType, MeshType TMeshType>
KRATOS_API(KRATOS_CORE) ContainerExpression<TContainerType, TMeshType> operator*(ContainerExpression<TContainerType, TMeshType> Left, const double Right);

template<class TContainerType, MeshType TMeshType>
KRATOS_API(KRATOS_CORE) ContainerExpression<TContainerType, TMeshType> operator*(const double Left, ContainerExpression<TContainerType, TMeshType> Right);

template<class TContainerType, MeshType TMeshType>
KRATOS_API(KRATOS_CORE) ContainerExpression<TContainerType, TMeshType> operator*(ContainerExpression<TContainerType, TMeshType> Left, ContainerExpression<TContainerType, TMeshType> Right);

template<class TContainerType, MeshType TMeshType>
KRATOS_API(KRATOS_CORE) ContainerExpression<TContainerType, TMeshType> operator/(ContainerExpression<TContainerType, TMeshType> Left, const double Right);

template<class TContainerType, MeshType TMeshType>
KRATOS_API(KRATOS_CORE) ContainerExpression<TContainerType, TMeshType> operator/(const double Left, ContainerExpression<TContainerType, TMeshType> Right);

template<class TContainerType, MeshType TMeshType>
KRATOS_API(KRATOS_CORE) ContainerExpression<TContainerType, TMeshType> operator/(ContainerExpression<TContainerType, TMeshType> Left, ContainerExpression<TContainerType, TMeshType> Right);

template<class TContainerType, MeshType TMeshType>
KRATOS_API(KRATOS_CORE) ContainerExpression<TContainerType, TMeshType> Power(const double Base, const ContainerExpression<TContainerType, TMeshType>& rExponent);

template<class TContainerType, MeshType TMeshType>
KRATOS_API(KRATOS_CORE) ContainerExpression<TContainerType, TMeshType> Power(ContainerExpression<TContainerType, TMeshType> Base, const double Exponent);

template<class TContainerType, MeshType TMeshType>
KRATOS_API(KRATOS_CORE) ContainerExpression<TContainerType, TMeshType> Power(ContainerExpression<TContainerType, TMeshType> Base, ContainerExpression<TContainerType, TMeshType> Exponent);

/// @}


} // namespace Kratos
