//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#ifndef KRATOS_KEY_HASH_H_INCLUDED
#define KRATOS_KEY_HASH_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/dof.h"
#include "utilities/indexed_object.h"

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
    typedef std::size_t IndexType;

    /// The definition of the hash type
    typedef std::size_t HashType;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

    /**
     * @brief This method creates an "unique" hash for the input value
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
     * @tparam TVariable The type of variable to be hashed
     */
    template<class TVariable>
    struct VariableHasher
    {
        /**
         * @brief This is the () operator
         * @param rVariable The variable to be hashed
         * @return The corresponding hash
         */
        HashType operator()(const TVariable& rVariable) const
        {
            return rVariable.Key();
        }
    };

    /**
     * @brief This is a hasher for variables pointers
     * @tparam TVariable The type of variable to be hashed
     */
    template<class TVariable>
    struct pVariableHasher
    {
        /**
         * @brief This is the () operator
         * @param pVariable The variable pointer to be hashed
         * @return The corresponding hash
         */
        HashType operator()(const TVariable* pVariable) const
        {
            return pVariable->Key();
        }
    };

    /**
     * @brief This is a key comparer between two variables
     * @tparam TVariable The type of variable to be compared
     */
    template<class TVariable>
    struct VariableComparator
    {
        /**
         * @brief This is the () operator
         * @param rFirst The first class to be compared
         * @param rSecond The second class to be compared
         */
        bool operator()(
            const TVariable& rFirst,
            const TVariable& rSecond
            ) const
        {
            return rFirst.Key() == rSecond.Key();
        }
    };

    /**
     * @brief This is a key comparer between two variables pointers
     * @tparam TVariable The type of variable to be compared
     */
    template<class TVariable>
    struct pVariableComparator
    {
        /**
         * @brief This is the () operator
         * @param pFirst The first class to be compared
         * @param pSecond The second class to be compared
         */
        bool operator()(
            const TVariable* pFirst,
            const TVariable* pSecond
            ) const
        {
            return pFirst->Key() == pSecond->Key();
        }
    };

    /**
     * @brief This is a hasher for indexed objects
     */
    struct IndexedObjectHasher
    {
        /**
         * @brief This is the () operator
         * @param rIndexedObject The indexed object to be hashed
         * @return The corresponding hash
         */
        HashType operator()(const IndexedObject& rIndexedObject) const
        {
            return rIndexedObject.Id();
        }
    };

    /**
     * @brief This is a key comparer between two indexed objects
     */
    struct IndexedObjectComparator
    {
        /**
         * @brief This is the () operator
         * @param rFirst The first class to be compared
         * @param rSecond The second class to be compared
         */
        bool operator()(
            const IndexedObject& rFirst,
            const IndexedObject& rSecond
            ) const
        {
            return rFirst.Id() == rSecond.Id();
        }
    };

    /**
     * @brief This is a hasher for indexed objects (pointer)
     * @param TpIndexedObject Pointer type to indexed object
     * @note Must be tenmplated to take into account the shared, intrussive,etc... pointers
     */
    template<class TpIndexedObject>
    struct IndexedObjecPointertHasher
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
     * @param TpIndexedObject Pointer type to indexed object
     * @note Must be tenmplated to take into account the shared, intrussive,etc... pointers
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
     * @brief This is a key comparer between two dof pointers
     * @details Used for example for the B&S
     */
    struct DofPointerComparor
    {
        /**
         * @brief The () operator
         * @param pDoF1 The first DoF pointer
         * @param pDoF2 The second DoF pointer
         */
        bool operator()(const Dof<double>::Pointer& pDoF1, const Dof<double>::Pointer& pDoF2) const
        {
            return (((pDoF1->Id() == pDoF2->Id() && (pDoF1->GetVariable()).Key()) == (pDoF2->GetVariable()).Key()));
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
        std::size_t operator()(const std::pair<TType1, TType2>& rPair) const
        {
            std::size_t seed = 0;
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
#endif
