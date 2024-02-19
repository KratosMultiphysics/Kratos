//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Suneth Warnakulasuriya
//
//

#pragma once

// System includes
#include <functional>

// External includes

// Project includes

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @class SetAddressIdentityFunction
 * @brief A functor that serves as the identity function for memory locations.
 * @details This class provides two overloaded operator() functions, one for non-const objects
 * and another for const objects. The operator() returns the input objects memory location as it is,
 * effectively acting as an identity function for address.
 * The purpose of this functor is to allow objects' memory locations to be used as keys in sets or other
 * containers that require comparison. By using this functor, you can avoid defining
 * custom comparison functions or operators when the object itself can be considered
 * for comparison.
 * @tparam TDataType The data type of the object that the functor operates on.
 * @author Pooyan Dadvand
 * @author Suneth Warnakulasuriya
 */
template<class TDataType>
class SetAddressIdentityFunction
{
public:
    /**
     * @brief Operator that returns a non-const reference pointer to the input object.
     * @param data The input object of type TDataType.
     * @return A non-const reference pointer to the same object as provided in the parameter.
     */
    TDataType* operator()(TDataType& data)
    {
        return &data;
    }

    /**
     * @brief Operator that returns a const reference pointer to the input object.
     * @param data The input object of type TDataType.
     * @return A const reference pointer to the same object as provided in the parameter.
     */
    TDataType const * operator()(const TDataType& data) const
    {
        return &data;
    }
};

///@}

}  // namespace Kratos.
