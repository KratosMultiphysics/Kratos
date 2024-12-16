//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/dof.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    /// The definition of the index type
    using IndexType = std::size_t;

    /// The definition of the hash type
    using HashType = std::size_t;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

    /**
     * @brief This method creates an "unique" hash for the input value
     * @details It comes from boost, taken from here: https://www.boost.org/doc/libs/1_55_0/doc/html/hash/reference.html#boost.hash_combine
     * @tparam TClassType The type of class to be hashed
     * @param Seed This is the seed used to create the hash
     * @param Value This is the value to be hashed
     */
    template <class TClassType>
    inline void HashCombine(
        HashType& Seed,
        const TClassType& Value
        )
    {
        std::hash<TClassType> hasher;
        Seed ^= hasher(Value) + 0x9e3779b9 + (Seed<<6) + (Seed>>2);
    }

    /**
     * @brief This method combines hash until it reaches the last class in order to obtain a corresponding seed
     * @tparam TClassType The type of class to be hashed
     * @param First The first class to be compared
     * @param Last The last class to be compared
     * @return The resulting seed
     */
    template <class TClassType>
    inline HashType HashRange(
        TClassType First,
        TClassType Last
        )
    {
        HashType seed = 0;

        while (First!=Last) {
            HashCombine(seed, *First);
            ++First;
        }

        return seed;
    }

///@}
///@name Kratos Classes
///@{

    class VariableData; // forward declaration

    /**
     * @brief This is a key comparer of general pourpose between two classes
     * @tparam TClassType The type of class to be hashed
     */
    template<class TClassType>
    struct KeyComparorRange
    {
        /**
         * @brief This is the () operator
         * @param first The first class to be compared
         * @param second The second class to be compared
         */
        bool operator()(
            const TClassType& first,
            const TClassType& second
            ) const
        {
            if(first.size() != second.size()) {
                return false;
            }

            auto it_first = first.begin();
            auto it_second = second.begin();

            while(it_first != first.end()) { // NOTE: We already checked that are same size
                if(*it_first != *it_second) {
                    return false;
                }
                if(it_first != first.end()) {
                    ++it_first;
                    ++it_second;
                }
            }

            return true;
        }
    };

    /**
     * @brief This is a hasher of general pourpose
     * @tparam TClassType The type of class to be hashed
     */
    template<class TClassType>
    struct KeyHasherRange
    {
        /**
         * @brief This is the () operator
         * @param rRange The shared pointer to be hashed
         * @return The corresponding hash
         */
        HashType operator()(const TClassType& rRange) const
        {
            return HashRange(rRange.begin(), rRange.end());
        }
    };

    /**
     * @brief This is a hasher for variables
     */
    struct KRATOS_API(KRATOS_CORE) VariableHasher
    {
        /**
         * @brief This is the () operator
         * @param rVariable The variable to be hashed
         * @return The corresponding hash
         */
        HashType operator()(const VariableData& rVariable) const;
    };

    /**
     * @brief This is a key comparer between two variables
     */
    struct KRATOS_API(KRATOS_CORE) VariableComparator
    {
        /**
         * @brief This is the () operator
         * @param rFirst The first class to be compared
         * @param rSecond The second class to be compared
         */
        bool operator()(
            const VariableData& rFirst,
            const VariableData& rSecond
            ) const;
    };

    /**
     * @brief This is a hasher for variables pointers
     */
    struct KRATOS_API(KRATOS_CORE) pVariableHasher
    {
        /**
         * @brief This is the () operator
         * @param pVariable The variable pointer to be hashed
         * @return The corresponding hash
         */
        HashType operator()(const VariableData* pVariable) const;
    };

    /**
     * @brief This is a key comparer between two variables pointers
     */
    struct KRATOS_API(KRATOS_CORE) pVariableComparator
    {
        /**
         * @brief This is the () operator
         * @param pFirst The first variable to be compared
         * @param pSecond The second variable to be compared
         */
        bool operator()(
            const VariableData* pFirst,
            const VariableData* pSecond
            ) const;
    };

    /**
     * @brief This is a hasher for indexed objects
     * @tparam TIndexedObject The type of indexed object
     */
    template<class TIndexedObject>
    struct IndexedObjectHasher
    {
        /**
         * @brief This is the () operator
         * @param rIndexedObject The indexed object to be hashed
         * @return The corresponding hash
         */
        HashType operator()(const TIndexedObject& rIndexedObject) const
        {
            return rIndexedObject.Id();
        }
    };

    /**
     * @brief This is a key comparer between two indexed objects
     * @tparam TIndexedObject The type of indexed object
     */
    template<class TIndexedObject>
    struct IndexedObjectComparator
    {
        /**
         * @brief This is the () operator
         * @param rFirst The first class to be compared
         * @param rSecond The second class to be compared
         */
        bool operator()(
            const TIndexedObject& rFirst,
            const TIndexedObject& rSecond
            ) const
        {
            return rFirst.Id() == rSecond.Id();
        }
    };

    /**
     * @brief This is a hasher for indexed objects (pointer)
     * @tparam TpIndexedObject Pointer type to indexed object
     * @note Must be tenmplated to take into account the shared, intrusive,etc... pointers
     */
    template<class TpIndexedObject>
    struct IndexedObjectPointerHasher
    {
        /**
         * @brief This is the () operator
         * @param pIndexedObject The indexed object pointer to be hashed
         * @return The corresponding hash
         */
        HashType operator()(const TpIndexedObject pIndexedObject) const
        {
            return pIndexedObject->Id();
        }
    };

    /**
     * @brief This is a key comparer between two indexed objects (pointer)
     * @tparam TpIndexedObject Pointer type to indexed object
     * @note Must be tenmplated to take into account the shared, intrusive,etc... pointers
     */
    template<class TpIndexedObject>
    struct IndexedObjectPointerComparator
    {
        /**
         * @brief This is the () operator
         * @param pFirst The first class to be compared
         * @param pSecond The second class to be compared
         */
        bool operator()(
            const TpIndexedObject pFirst,
            const TpIndexedObject pSecond
            ) const
        {
            return pFirst->Id() == pSecond->Id();
        }
    };

    /**
     * @brief This is a hasher for shared pointers
     * @tparam TSharedPointer The type of shared pointer to be hashed
     */
    template<class TSharedPointer>
    struct SharedPointerHasher
    {
        /**
         * @brief This is the () operator
         * @param pPointer The shared pointer to be hashed
         * @return The corresponding hash
         */
        HashType operator()(const TSharedPointer& pPointer) const
        {
            return reinterpret_cast<HashType>(pPointer.get());
        }
    };

    /**
     * @brief This is a key comparer between two shared pointers
     * @tparam TSharedPointer The type of shared pointer to be compared
     */
    template<class TSharedPointer>
    struct SharedPointerComparator
    {
        /**
         * @brief This is the () operator
         * @param first The first class to be compared
         * @param second The second class to be compared
         */
        bool operator()(
            const TSharedPointer& first,
            const TSharedPointer& second
            ) const
        {
            return first.get() == second.get();
        }
    };

    /**
     * @brief This is a key comparer between two pointers
     * @details This compares two pointers by using its value comparison operators.
     */
    template<class TDataType>
    struct PointerComparor
    {
        // Data type used in the comparison
        using data_type = TDataType;
        /**
         * @brief The () operator
         * @param first  The first pointer
         * @param second The second pointer
         */
        bool operator()(TDataType const * first, TDataType const * second) const
        {
            return (*first) == (*second);
        }
    };    

    /**
     * @brief This is a hasher between two vectors of indexes
     * @tparam TVectorIndex The type of vector indexes to be compared
     */
    template<class TVectorIndex>
    struct VectorIndexHasher
    {
        /**
         * @brief This is the () operator
         * @param k The vector of indexes to be hashed
         * @return The corresponding hash
         */
        HashType operator()(const TVectorIndex& k) const
        {
            return HashRange(k.begin(), k.end());
        }
    };

    /**
     * @brief This is a key comparer between two vectors of indexes
     * @tparam TVectorIndex The type of vector indexes to be compared
     */
    template<class TVectorIndex>
    struct VectorIndexComparor
    {
        /**
         * @brief This is the () operator
         * @param lhs The first class to be compared
         * @param rhs The second class to be compared
         */
        bool operator()(const TVectorIndex& lhs, const TVectorIndex& rhs) const
        {
            if(lhs.size() != rhs.size())
                return false;

            for(IndexType i = 0; i < lhs.size(); i++) {
                if(lhs[i] != rhs[i]) return false;
            }

            return true;
        }
    };

    /**
     * @brief This is a hasher for a dof pointers
     * @details Used for example for the B&S
     */
    struct DofPointerHasher
    {
        /**
         * @brief The () operator
         * @param pDoF The DoF pointer
         */
        HashType operator()(const Dof<double>::Pointer& pDoF) const
        {
            HashType seed = 0;
            HashCombine(seed, pDoF->Id());
            HashCombine(seed, (pDoF->GetVariable()).Key());
            return seed;
        }
    };

    /**
     * @brief This is a hasher for pairs
     * @details Used for example for edges ids
     */
    template<class TType1, class TType2>
    struct PairHasher
    {
        /**
         * @brief The () operator
         * @param rPair The index pair to hash
         */
        HashType operator()(const std::pair<TType1, TType2>& rPair) const
        {
            HashType seed = 0;
            HashCombine<TType1>(seed, rPair.first);
            HashCombine<TType2>(seed, rPair.second);
            return seed;
        }
    };

    /**
     * @brief This is a key comparer between two indexes pairs
     * @details Used for example for the B&S
     */
    template<class TType1, class TType2>
    struct PairComparor
    {
        /**
         * @brief The () operator
         * @param rIndexPair1 The first index pair to hash
         * @param rIndexPair2 The second index pair to hash
         */
        bool operator()(const std::pair<TType1, TType2>& rPair1, const std::pair<TType1, TType2>& rPair2) const
        {
            return ((std::get<0>(rPair1) == std::get<0>(rPair2)) && (std::get<1>(rPair1) == std::get<1>(rPair2)));
        }
    };

///@}
///@name Kratos Classes
///@{

} // namespace Kratos.

/**
 * @brief This defines the missing hashs for the std namespace
*/
namespace std 
{
    /**
     * @brief This is a hasher for pairs
     * @details Used for example for edges ids
     * @tparam T1 The first type of the pair
     * @tparam T2 The second type of the pair
     * @note This is needed to use pairs as keys in unordered maps
     */
    template<typename T1, typename T2>
    struct hash<std::pair<T1, T2>> 
    {
        /**
         * @brief Calculates the hash value of a given pair of values in a way that combines the hash values of the individual elements
         * @param p the pair of values to be hashed.
         * @return the resulting hash value.
         */
        size_t operator()(const std::pair<T1, T2>& p) const 
        {
            size_t seed = 0;
            const size_t h1 = std::hash<T1>()(p.first);
            const size_t h2 = std::hash<T2>()(p.second);
            Kratos::HashCombine(seed, h1);
            Kratos::HashCombine(seed, h2);
            return seed;
        }
    };
} // namespace std.
