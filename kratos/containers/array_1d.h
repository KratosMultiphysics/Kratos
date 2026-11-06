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
//                   Riccardo Rossi
//                   Vicente Mataix Ferrandiz
//

#pragma once

// The fixed-size dense vector of Kratos. The implementation follows the
// linear-algebra backend selected at configure time
// (KRATOS_LINEAR_ALGEBRA_BACKEND): an Eigen fixed-size vector under the Eigen
// backend, the historical boost::numeric::ublas vector_expression otherwise.
// Both implementations honour the same contract (see
// tests/cpp_tests/containers/test_dense_backend_contract.cpp): plain T[N]
// layout starting at the object address, the uBLAS-style member surface
// (operator()/operator[], size, resize, iterators, ...) and the uBLAS text
// format for streaming.

#ifdef KRATOS_USE_EIGEN_BACKEND
#include "containers/eigen_array_1d.h"
#else
#include "containers/ublas_array_1d.h"
#endif

namespace AuxiliaryHashCombine
{
    /**
     * @brief This method creates an "unique" hash for the input value
     * @details It comes from boost, taken from here: https://www.boost.org/doc/libs/1_55_0/doc/html/hash/reference.html#boost.hash_combine
     * @tparam TClassType The type of class to be hashed
     * @param rSeed This is the seed used to create the hash
     * @param rValue This is the value to be hashed
     * @todo Once the hashers and comparors are moved from key_hash.h, include key_hash and remove this. Right now there is a cross inclussion
     */
    template <class TClassType>
    inline void HashCombine(
        std::size_t& rSeed,
        const TClassType& rValue
        )
    {
        std::hash<TClassType> hasher;
        rSeed ^= hasher(rValue) + 0x9e3779b9 + (rSeed<<6) + (rSeed>>2);
    }
} /// namespace AuxiliaryHashCombine

namespace std
{
template<class T, std::size_t N>
struct hash<Kratos::array_1d<T,N>>
{
    std::size_t operator()(const Kratos::array_1d<T,N>& rArray) {
            std::size_t seed = 0;
            for (auto component : rArray) {AuxiliaryHashCombine::HashCombine(seed, component);}
            return seed;
        }
};
} /// namespace std
