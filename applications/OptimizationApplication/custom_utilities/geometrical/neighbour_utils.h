//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <unordered_map>
#include <vector>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

// Application includes

namespace Kratos {

///@name Kratos Classes
///@{

class KRATOS_API(OPTIMIZATION_APPLICATION) NeighbourUtils {
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    ///@}
    ///@name Static operations
    ///@{

    /**
     * @brief This method initializes NEIGHBOUR_ELEMENTS variable of conditions with its parent element
     *
     * This method initializes the NEIGHBOUR_ELEMENTS variable of each condition with one item, which
     * is the parent element which a condition is sharing a face (in 3D), an edge (in 2D), point (in 1D).
     *
     * The input model part must contain both elements and conditions (including each conditions parents).
     *
     * This should work in MPI as well, because, MPI mesh partitioner always keeps conditions and their parents
     * in the same partition (i.e. rank), hence the methodology does not differ if run in MPI or OpenMP.
     *
     * This is an expensive initialization, so number of calls to this method should be kept minimum.
     *
     * @param rModelPart Input model part which contains conditions and their parent elements.
     */
    static void InitializeParentElementForConditions(ModelPart& rModelPart);

    /**
     * @brief Get the Condition Id And Parent Element Ids Map using the initialized NEIGHBOUR_ELEMENTS varaible of conditions
     *
     * This method gives a map of condition id, condition parent element ids for the given rModelPart. This method
     * should be only called once the NEIGHBOUR_ELEMENTS of the given model part conditions have been properly initialized.
     *
     * @see InitializeParentElementForConditions
     * @param rModelPart
     * @return std::unordered_map<IndexType, std::vector<IndexType>>
     */
    static std::unordered_map<IndexType, std::vector<IndexType>> GetConditionIdAndParentElementIdMap(const ModelPart& rModelPart);

    ///@}
private:
    ///@name Classes
    ///@{

    /**
     * @brief Unordered map reducer for shared memory parallelized loops.
     *
     * @tparam KeyType
     * @tparam ValueType
     */
    template<class KeyType, class ValueType, class Hasher = std::hash<KeyType>, class Comparator = std::equal_to<KeyType>>
    class HashMapReducer
    {
    public:
        using value_type = std::pair<KeyType, ValueType>;
        using return_type = std::unordered_map<KeyType, ValueType, Hasher, Comparator>;

        return_type mValue;

        /// access to reduced value
        return_type GetValue() const
        {
            return mValue;
        }

        /// NON-THREADSAFE (fast) value of reduction, to be used within a single thread
        void LocalReduce(value_type rValue){
            mValue.emplace(rValue);
        }

        /// THREADSAFE (needs some sort of lock guard) reduction, to be used to sync threads
        void ThreadSafeReduce(const HashMapReducer& rOther)
        {
            KRATOS_CRITICAL_SECTION
            for (const auto& it : rOther.mValue) {
                mValue.emplace(it);
            }
        }
    };

    ///@}
};

///@}
} // namespace Kratos