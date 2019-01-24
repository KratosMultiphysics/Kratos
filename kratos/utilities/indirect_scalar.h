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
#define  KRATOS_INDIRECT_SCALAR_H_INCLUDED

// System includes
#include <functional>

// Project includes
#include "utilities/indirect_scalar_fwd.h"
#include "includes/node.h"

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
 *     auto fget = [&a]() -> double { return a[1]; };
 *     auto fset = [&a](double d) { a[1] = d; };
 *     IndirectScalar<double> f{fset, fget};
 *     std::cout << f << std::endl;
 *     f = 0.0;
 *     std::cout << f << std::endl;
 */
template <class T, typename S>
class IndirectScalar
{
public:
    IndirectScalar() : set([](T) {}), get([]() -> T { return T{}; })
    {
    }

    IndirectScalar(std::function<void(T)> set, std::function<T()> get)
        : set(set), get(get)
    {
    }

    IndirectScalar<T, S>& operator=(const T value)
    {
        set(value);
        return *this;
    }

    IndirectScalar<T, S>& operator+=(const T value)
    {
        set(get() + value);
        return *this;
    }

    IndirectScalar<T, S>& operator-=(const T value)
    {
        set(get() - value);
        return *this;
    }

    IndirectScalar<T, S>& operator*=(const T value)
    {
        set(get() * value);
        return *this;
    }

    IndirectScalar<T, S>& operator/=(const T value)
    {
        set(get() / value);
        return *this;
    }

    operator T() const
    {
        return get();
    }

    std::ostream& print(std::ostream& os) const
    {
        return os << get();
    }

private:
    std::function<void(T)> set;
    std::function<T()> get;

    friend class Serializer;

    void save(Serializer& rSerializer) const
    {
    }

    void load(Serializer& rSerializer)
    {
    }
};

template <class T>
std::ostream& operator<<(std::ostream& os, const IndirectScalar<T>& s)
{
    return s.print(os);
}

///@} // Kratos Classes

template <class TVariableType>
IndirectScalar<typename TVariableType::Type> MakeIndirectScalar(Node<3>& rNode,
                                                                const TVariableType& rVariable)
{
    auto fset = [&rNode, &rVariable](typename TVariableType::Type s) {
        rNode.FastGetSolutionStepValue(rVariable) = s;
    };

    auto fget = [&rNode, &rVariable]() -> typename TVariableType::Type {
        return rNode.FastGetSolutionStepValue(rVariable);
    };

    return IndirectScalar<typename TVariableType::Type>(fset, fget);
}

template <class TVariableType>
IndirectScalar<typename TVariableType::Type> MakeIndirectScalar(Node<3>& rNode,
                                                                const TVariableType& rVariable,
                                                                std::size_t Step)
{
    if (Step == 0)
    {
        return MakeIndirectScalar(rNode, rVariable);
    }
    // Here we don't capture Step to avoid allocating memory.
    else if (Step == 1)
    {
        auto fset = [&rNode, &rVariable](typename TVariableType::Type s) {
            rNode.FastGetSolutionStepValue(rVariable, 1) = s;
        };
        auto fget = [&rNode, &rVariable]() -> typename TVariableType::Type {
            return rNode.FastGetSolutionStepValue(rVariable, 1);
        };
        return IndirectScalar<typename TVariableType::Type>(fset, fget);
    }
    else if (Step == 2)
    {
        auto fset = [&rNode, &rVariable](typename TVariableType::Type s) {
            rNode.FastGetSolutionStepValue(rVariable, 2) = s;
        };
        auto fget = [&rNode, &rVariable]() -> typename TVariableType::Type {
            return rNode.FastGetSolutionStepValue(rVariable, 2);
        };
        return IndirectScalar<typename TVariableType::Type>(fset, fget);
    }
    else
    {
        KRATOS_ERROR << "Invalid Step = " << Step << std::endl;
        return IndirectScalar<typename TVariableType::Type>{};
    }
}

}  // namespace Kratos.

#endif // KRATOS_INDIRECT_SCALAR_H_INCLUDED  defined 
