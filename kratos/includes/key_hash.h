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

namespace Kratos 
{
///@addtogroup Kratos Core
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

    /**
     * This method creates an "unique" hash for the input value
     * @param TClassType The type of class to be hashed
     * @param Seed This is the seed used to create the hash
     * @param Value This is the value to be hashed
     */
    template <class TClassType>
    inline void HashCombine(
        std::size_t& Seed, 
        const TClassType& Value
        )
    {
        std::hash<TClassType> hasher;
        Seed ^= hasher(Value) + 0x9e3779b9 + (Seed<<6) + (Seed>>2);
    }
    
    /**
     * This method combines hash until it reaches the last class in order to obtain a corresponding seed
     * @param TClassType The type of class to be hashed
     * @param First The first class to be compared
     * @param Last The last class to be compared
     * @return The resulting seed
     */
    template <class TClassType>
    inline std::size_t HashRange(
        TClassType First, 
        TClassType Last
        )
    {
        std::size_t seed = 0;

        while (First!=Last)
        {
            HashCombine(seed, *First);
            ++First;
        }
        
        return seed;
    }
    
    /**
     * \brief This is a key comparer of general pourpose between two classes 
     * @param TClassType The type of class to be hashed
     */
    template<class TClassType>
    struct KeyComparorRange
    {
        /**
         * This is the () operator
         * @param first The first class to be compared
         * @param second The second class to be compared
         */
        bool operator()(
            const TClassType& first, 
            const TClassType& second
            ) const
        {
            if(first.size() != second.size())
            {
                return false;
            }

            auto it_first = first.begin();
            auto it_second = second.begin();

            while(it_first != first.end()) // NOTE: We already checked that are same size
            {
                if(*it_first != *it_second) 
                {
                    return false;
                }
                if(it_first != first.end())
                {
                    ++it_first;
                    ++it_second;
                }
            }

            return true;
        }
    };
    
    /**
     * \brief This is a hasher of general pourpose 
     * @param TClassType The type of class to be hashed
     */
    template<class TClassType>
    struct KeyHasherRange
    {
        /**
         * This is the () operator
         * @param pPointer The shared pointer to be hashed
         * @return The corresponding hash
         */
        std::size_t operator()(const TClassType& rRange) const
        {
            return HashRange(rRange.begin(), rRange.end());
        }
    };
    
    /**
     * \brief This is a hasher for shared pointers 
     * @param TSharedPointer The type of shared pointer to be hashed
     */
    template<class TSharedPointer>
    struct SharedPointerHasher
    {
        /**
         * This is the () operator
         * @param pPointer The shared pointer to be hashed
         * @return The corresponding hash
         */
        std::size_t operator()(const TSharedPointer& pPointer) const
        {
            return reinterpret_cast<std::size_t>(pPointer.get());
        }
    };
    
    /**
     * \brief This is a key comparer between two shared pointers 
     * @param TClassType The type of shared pointer to be compared
     */
    template<class TSharedPointer>
    struct SharedPointerComparator
    {
        /**
         * This is the () operator
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
    
///@}
///@name Kratos Classes
///@{

} // namespace Kratos.
#endif
