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

// System includes
#include <functional>
#include <type_traits>
#include <utility>

// External includes

// Project includes

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @class TestOperators
 * @brief Operator class to test whether operators "<" and "==" is defined for given TDataType.
 */
struct TestOperators
{
    template<class TDataType>
    static auto TestEqual(TDataType*) -> decltype(std::declval<TDataType>() == std::declval<TDataType>());
    template<typename>
    static auto TestEqual(...) -> std::false_type;

    template<class TDataType>
    static auto TestLess(TDataType*) -> decltype(std::declval<TDataType>() < std::declval<TDataType>());
    template<typename>
    static auto TestLess(...) -> std::false_type;
};

/**
 * @class KeyGenerator
 * @brief A functor that serves as the identity function if "<" and "==" operators are defined, otherwise, returns memory location.
 * @details This class provides two overloaded operator() functions, one for non-const objects
 * and another for const objects. The operator() returns the input objects memory location as it is,
 * effectively acting as an identity function for address if the operators "<" and "==" is not defined. Otherwise,
 * this returns the exact input object.
 *
 * The purpose of this functor is to allow objects' memory locations or exact object to be used as keys in sets or other
 * containers that require comparison. By using this functor, you can avoid defining
 * custom comparison functions or operators when the object itself can be considered
 * for comparison.
 * @tparam TDataType The data type of the object that the functor operates on.
 * @author Suneth Warnakulasuriya
 */
template<class TDataType>
class KeyGenerator
{
private:
    ///@name Private static variables
    ///@{

    static constexpr bool HasOperatorEqual = std::is_same_v<bool, decltype(TestOperators::TestEqual<TDataType>(0))>;

    static constexpr bool HasOperatorLess = std::is_same_v<bool, decltype(TestOperators::TestLess<TDataType>(0))>;

    static constexpr bool HasOperatorsDefined = HasOperatorEqual && HasOperatorLess;

    ///@}

public:
    ///@name Public operators
    ///@{

    /**
     * @brief Operator that returns a non-const reference to the input object if "<" and "==" operators are defined.
     * Otherwise this returns a pointer to the input object.
     * @param data The input object of type TDataType.
     * @return A non-const reference or pointer to the same object as provided in the parameter.
     */
    std::conditional_t<HasOperatorsDefined, TDataType&, TDataType*> operator()(TDataType& rData)
    {
        if constexpr(HasOperatorsDefined) {
            return rData;
        } else {
            return &rData;
        }
    }

    /**
     * @brief Operator that returns a const reference to the input object if "<" and "==" operators are defined.
     * Otherwise this returns a pointer to the const input object.
     * @param data The input object of type TDataType.
     * @return A const reference or pointer to the same object as provided in the parameter.
     */
    std::conditional_t<HasOperatorsDefined, TDataType const&, TDataType const*> operator()(const TDataType& rData) const
    {
        if constexpr(HasOperatorsDefined) {
            return rData;
        } else {
            return &rData;
        }
    }

    ///@}
};

///@}

}  // namespace Kratos.
