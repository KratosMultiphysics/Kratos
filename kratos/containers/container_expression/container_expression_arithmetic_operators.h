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

namespace Kratos {


/// @name Arithmetic Operators
/// @{

template<class TContainerType, class TMeshType>
class ContainerExpression;

template<class TContainerType, class TMeshType>
ContainerExpression<TContainerType, TMeshType> operator+(const ContainerExpression<TContainerType, TMeshType>& rLeft, const double Right);

template<class TContainerType, class TMeshType>
ContainerExpression<TContainerType, TMeshType> operator+(const double Left, const ContainerExpression<TContainerType, TMeshType>& rRight);

template<class TContainerType, class TMeshType>
ContainerExpression<TContainerType, TMeshType> operator+(const ContainerExpression<TContainerType, TMeshType>& rLeft, const ContainerExpression<TContainerType, TMeshType>& rRight);

template<class TContainerType, class TMeshType>
ContainerExpression<TContainerType, TMeshType> operator-(const ContainerExpression<TContainerType, TMeshType>& rLeft, const double Right);

template<class TContainerType, class TMeshType>
ContainerExpression<TContainerType, TMeshType> operator-(const double Left, const ContainerExpression<TContainerType, TMeshType>& rRight);

template<class TContainerType, class TMeshType>
ContainerExpression<TContainerType, TMeshType> operator-(const ContainerExpression<TContainerType, TMeshType>& rLeft, const ContainerExpression<TContainerType, TMeshType>& rRight);

template<class TContainerType, class TMeshType>
ContainerExpression<TContainerType, TMeshType> operator*(const ContainerExpression<TContainerType, TMeshType>& rLeft, const double Right);

template<class TContainerType, class TMeshType>
ContainerExpression<TContainerType, TMeshType> operator*(const double Left, const ContainerExpression<TContainerType, TMeshType>& rRight);

template<class TContainerType, class TMeshType>
ContainerExpression<TContainerType, TMeshType> operator*(const ContainerExpression<TContainerType, TMeshType>& rLeft, const ContainerExpression<TContainerType, TMeshType>& rRight);

template<class TContainerType, class TMeshType>
ContainerExpression<TContainerType, TMeshType> operator/(const ContainerExpression<TContainerType, TMeshType>& rLeft, const double Right);

template<class TContainerType, class TMeshType>
ContainerExpression<TContainerType, TMeshType> operator/(const double Left, const ContainerExpression<TContainerType, TMeshType>& rRight);

template<class TContainerType, class TMeshType>
ContainerExpression<TContainerType, TMeshType> operator/(const ContainerExpression<TContainerType, TMeshType>& rLeft, const ContainerExpression<TContainerType, TMeshType>& rRight);

template<class TContainerType, class TMeshType>
ContainerExpression<TContainerType, TMeshType> Power(const double Base, const ContainerExpression<TContainerType, TMeshType>& rExponent);

template<class TContainerType, class TMeshType>
ContainerExpression<TContainerType, TMeshType> Power(const ContainerExpression<TContainerType, TMeshType>& rBase, const double Exponent);

template<class TContainerType, class TMeshType>
ContainerExpression<TContainerType, TMeshType> Power(const ContainerExpression<TContainerType, TMeshType>& rBase, const ContainerExpression<TContainerType, TMeshType>& rExponent);

/// @}


} // namespace Kratos
