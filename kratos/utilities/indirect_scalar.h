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

#if !defined(KRATOS_INDIRECT_SCALAR_H_INCLUDED)
#define KRATOS_INDIRECT_SCALAR_H_INCLUDED

// System includes

// Project includes
#include "includes/node.h"
#include "utilities/indirect_scalar_fwd.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @class IndirectScalar
 * @ingroup KratosCore
 * @brief Wrapper for a function which behaves like an arithmetic type.
 * @details Example:
 *     std::array<double, 3> a = {1.0, 2.0, 3.0};
 *     IndirectScalar<double> f{&a[1]};
 *     std::cout << f << std::endl;
 *     f = 0.0;
 *     std::cout << f << std::endl;
 */
template <class TDataType>
class IndirectScalar
{
public:
    IndirectScalar();

    IndirectScalar(TDataType* pVariableValue);

    IndirectScalar<TDataType>& operator=(const TDataType value);

    IndirectScalar<TDataType>& operator+=(const TDataType value);

    IndirectScalar<TDataType>& operator-=(const TDataType value);

    IndirectScalar<TDataType>& operator*=(const TDataType value);

    IndirectScalar<TDataType>& operator/=(const TDataType value);

    operator TDataType() const;

    TDataType* const pGetValue() const;

    std::ostream& print(std::ostream& os) const;

private:
    static TDataType mZero;
    static TDataType mGarbageValue;
    TDataType* mpGetVariable;
    TDataType* mpSetVariable;

    friend class Serializer;

    void save(Serializer& rSerializer) const
    {
    }

    void load(Serializer& rSerializer)
    {
    }
};

template <class TDataType>
std::ostream& operator<<(std::ostream& os, const IndirectScalar<TDataType>& s)
{
    return s.print(os);
}

///@} // Kratos Classes

template <class TVariableType>
IndirectScalar<typename TVariableType::Type> MakeIndirectScalar(Node<3>& rNode,
                                                                const TVariableType& rVariable,
                                                                std::size_t Step = 0)
{
    return IndirectScalar<typename TVariableType::Type>(
        &rNode.FastGetSolutionStepValue(rVariable, Step));
}

} // namespace Kratos.

#endif // KRATOS_INDIRECT_SCALAR_H_INCLUDED  defined
