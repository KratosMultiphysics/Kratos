//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:
//
//

// Project includes
#include "utilities/indirect_scalar.h"

namespace Kratos
{
template <class TDataType>
TDataType IndirectScalar<TDataType>::mZero = TDataType{};

template <class TDataType>
TDataType IndirectScalar<TDataType>::mGarbageValue = TDataType{};

template <class TDataType>
IndirectScalar<TDataType>::IndirectScalar(TDataType* pVariableValue)
    : mpGetVariable(pVariableValue), mpSetVariable(pVariableValue)
{
}

template <class TDataType>
IndirectScalar<TDataType>::IndirectScalar()
    : mpGetVariable(&IndirectScalar::mZero), mpSetVariable(&IndirectScalar::mGarbageValue)
{
}

template <class TDataType>
IndirectScalar<TDataType>& IndirectScalar<TDataType>::operator=(const TDataType value)
{
    *mpSetVariable = value;
    return *this;
}

template <class TDataType>
IndirectScalar<TDataType>& IndirectScalar<TDataType>::operator+=(const TDataType value)
{
    *mpSetVariable = *mpGetVariable + value;
    return *this;
}

template <class TDataType>
IndirectScalar<TDataType>& IndirectScalar<TDataType>::operator-=(const TDataType value)
{
    *mpSetVariable = *mpGetVariable - value;
    return *this;
}

template <class TDataType>
IndirectScalar<TDataType>& IndirectScalar<TDataType>::operator*=(const TDataType value)
{
    *mpSetVariable = *mpGetVariable * value;
    return *this;
}

template <class TDataType>
IndirectScalar<TDataType>& IndirectScalar<TDataType>::operator/=(const TDataType value)
{
    *mpSetVariable = *mpGetVariable / value;
    return *this;
}

template <class TDataType>
IndirectScalar<TDataType>::operator TDataType() const
{
    return *mpGetVariable;
}

template <class TDataType>
TDataType* const IndirectScalar<TDataType>::pGetValue() const
{
    return mpGetVariable;
}

template <class TDataType>
std::ostream& IndirectScalar<TDataType>::print(std::ostream& os) const
{
    return os << *mpGetVariable;
}

// Template initialization

template class IndirectScalar<double>;

} // namespace Kratos
