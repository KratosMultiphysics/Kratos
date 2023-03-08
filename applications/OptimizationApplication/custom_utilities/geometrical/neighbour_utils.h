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
     * @brief This method initializes NEIGHBOUR_ELEMENTS variable of conditions with its parent elements
     *
     * This method initializes the NEIGHBOUR_ELEMENTS variable of each condition with their parent elements.
     *
     * In elements with LocalSpaceDimension == 3, this assign parent elements for surface conditions, line conditions
     * and point conditions
     *
     * In elements with LocalSpaceDimension == 2, this assigns parent elements for line conditions and point conditions.
     *
     * The input model part must contain both elements and conditions (including each conditions parents).
     *
     * This should work in MPI as well, because, MPI mesh partitioner always keeps conditions and their parents
     * in the same partition (i.e. rank), hence the methodology does not differ if run in MPI or OpenMP.
     *
     * In the case of Point conditions, this assigns only the local parents to the node (or ghost node). Hence,
     * the computations done on to those nodes, should be assembled to have MPI compatibility.
     *
     * This is an expensive initialization, so number of calls to this method should be kept minimum.
     *
     * @param rModelPart Input model part which contains conditions and their parent elements.
     */
    static void InitializeParentElementsForConditions(ModelPart& rModelPart);

    /**
     * @brief Get the Condition Id And Parent Element Ids Map using the initialized NEIGHBOUR_ELEMENTS varaible of conditions
     *
     * This method gives a map of condition id, condition parent element ids for the given rModelPart. This method
     * should be only called once the NEIGHBOUR_ELEMENTS of the given model part conditions have been properly initialized.
     *
     * @see InitializeParentElementsForConditions
     * @param rModelPart
     * @return std::unordered_map<IndexType, std::vector<IndexType>>
     */
    static std::unordered_map<IndexType, std::vector<IndexType>> GetConditionIdAndParentElementIdsMap(const ModelPart& rModelPart);

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
    template<class MapType>
    class MapReduction
    {
    public:
        using value_type = typename MapType::value_type;
        using return_type = MapType;

        return_type mValue;

        /// access to reduced value
        return_type GetValue() const
        {
            return mValue;
        }

        /// NON-THREADSAFE (fast) value of reduction, to be used within a single thread
        void LocalReduce(const value_type rValue){
            mValue.emplace(rValue);
        }

        /// THREADSAFE (needs some sort of lock guard) reduction, to be used to sync threads
        void ThreadSafeReduce(const MapReduction<MapType>& rOther)
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