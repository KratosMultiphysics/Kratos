//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Denis Demidov
//

#if !defined(KRATOS_PARALLEL_UTILITIES_H_INCLUDED)
#define KRATOS_PARALLEL_UTILITIES_H_INCLUDED

// System includes
#include <iostream>
#include <array>
#include <vector>
#include <tuple>
#include <cmath>
#include <limits>
#include <future>
#include <thread>

// External includes
#include <omp.h>

// Project includes
#include "includes/define.h"
#include "includes/global_variables.h"
#include "utilities/reduction_utilities.h"


namespace Kratos
{
///@addtogroup KratosCore

//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
/** @param TContainerType - the type of the container used in the loop (must provide random access iterators)
 *  @param TIteratorType - type of iterator (by default as provided by the TContainerType)
 *  @param TMaxThreads - maximum number of threads allowed in the partitioning.
 *                       must be known at compile time to avoid heap allocations in the partitioning
 */
template<
        class TContainerType,
        class TIteratorType=typename TContainerType::iterator,
        int TMaxThreads=Globals::MaxAllowedThreads
        >
class BlockPartition
{
public:

    /** @param it_begin - iterator pointing at the beginning of the container
     *  @param it_end - iterator pointing to the end of the container
     *  @param Nchunks - number of threads to be used in the loop (must be lower than TMaxThreads)
     */
    BlockPartition(TIteratorType it_begin,
                   TIteratorType it_end,
                   int Nchunks = omp_get_max_threads())
        : mNchunks(Nchunks)
    {
        ptrdiff_t mBlockPartitionSize = (it_end-it_begin) / mNchunks;
        mBlockPartition[0] = it_begin;
        mBlockPartition[mNchunks] = it_end;
        for (int i=1; i<mNchunks; i++) {
            mBlockPartition[i] = mBlockPartition[i-1] + mBlockPartitionSize;
        }
    }

    /** @param rData - the continer to be iterated upon
     *  @param Nchunks - number of threads to be used in the loop (must be lower than TMaxThreads)
     */
    BlockPartition(TContainerType& rData,
                   int Nchunks = omp_get_max_threads())
        : BlockPartition(rData.begin(), rData.end(), Nchunks)
    {}

    virtual ~BlockPartition() = default;

    /** @brief simple iteration loop. f called on every entry in rData
     * @param f - must be a unary function accepting as input TContainerType::value_type&
     */
    template <class TUnaryFunction>
    inline void for_each(TUnaryFunction&& f)
    {
        #pragma omp parallel for
        for (int i=0; i<mNchunks; ++i) {
            for (auto it = mBlockPartition[i]; it != mBlockPartition[i+1]; ++it) {
                f(*it); //note that we pass the value to the function, not the iterator
            }
        }
    }

    /** @brie loop allowing reductions. f called on every entry in rData
     * the function f needs to return the values to be used by the reducer
     * @param TReducer template parameter specifying the reduction operation to be done
     * @param f - must be a unary function accepting as input TContainerType::value_type&
     */
    template <class TReducer, class TUnaryFunction>
    inline typename TReducer::value_type for_each(TUnaryFunction &&f)
    {
        TReducer global_reducer;
        #pragma omp parallel for
        for (int i=0; i<mNchunks; ++i) {
            TReducer local_reducer;
            for (auto it = mBlockPartition[i]; it != mBlockPartition[i+1]; ++it) {
                local_reducer.LocalReduce(f(*it));
            }
            global_reducer.ThreadSafeReduce(local_reducer);
        }
        return global_reducer.GetValue();
    }

private:
    int mNchunks;
    std::array<TIteratorType, TMaxThreads> mBlockPartition;
};

/** @brief simplified version of the basic loop (without reduction) to enable template type deduction
 * @param v - containers to be looped upon
 * @param func - must be a unary function accepting as input TContainerType::value_type&
 *
 */
template <class TContainerType, class TFunctionType>
void block_for_each(TContainerType &&v, TFunctionType &&func)
{
    BlockPartition<typename std::decay<TContainerType>::type>(std::forward<TContainerType>(v)).for_each(std::forward<TFunctionType>(func));
}

template <class TReducer, class TContainerType, class TFunctionType>
typename TReducer::value_type block_for_each(TContainerType &&v, TFunctionType &&func)
{
    return BlockPartition<typename std::decay<TContainerType>::type>
        (std::forward<TContainerType>(v)).template for_each<TReducer>(std::forward<TFunctionType>(func));
}

//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
/** @brief This class is useful for index iteration over containers
 *  @param TIndexType type of index to be used in the loop
 *  @param TMaxThreads - maximum number of threads allowed in the partitioning.
 *                       must be known at compile time to avoid heap allocations in the partitioning
 */
template<class TIndexType, int TMaxThreads=Globals::MaxAllowedThreads>
class IndexPartition
{
public:

/** @brief constructor using the size of the partition to be used
 *  @param Size - the size of the partition
 *  @param Nchunks - number of threads to be used in the loop (must be lower than TMaxThreads)
 */
    IndexPartition(TIndexType Size,
                   int Nchunks = omp_get_max_threads())
        : mNchunks(Nchunks)
    {
        int mBlockPartitionSize = Size / mNchunks;
        mBlockPartition[0] = 0;
        mBlockPartition[mNchunks] = Size;
        for (int i=1; i<mNchunks; i++) {
            mBlockPartition[i] = mBlockPartition[i-1] + mBlockPartitionSize;
        }

    }

    virtual ~IndexPartition() = default;

    //NOT COMMENTING IN DOXYGEN - THIS SHOULD BE SORT OF HIDDEN UNTIL GIVEN PRIME TIME
    //pure c++11 version (can handle exceptions)
    template <class TUnaryFunction>
    inline void for_pure_c11(TUnaryFunction &&f)
    {
        std::vector< std::future<void> > runners(mNchunks);
        const auto& partition = mBlockPartition;
        for (int i=0; i<mNchunks; ++i) {
            runners[i] = std::async(std::launch::async, [&partition, i,  &f]()
                {
                    for (auto k = partition[i]; k < partition[i+1]; ++k) {
                        f(k);
                    }
                });
        }

        //here we impose a syncronization and we check the exceptions
        for(int i=0; i<mNchunks; ++i) {
            try {
                runners[i].get();
            }
            catch(Exception& e) {
                KRATOS_ERROR << std::endl << "THREAD number: " << i << " caught exception " << e.what() << std::endl;
            } catch(std::exception& e) {
                KRATOS_ERROR << std::endl << "THREAD number: " << i << " caught exception " << e.what() << std::endl;
            } catch(...) {
                KRATOS_ERROR << std::endl << "unknown error" << std::endl;
            }
        }
    }

    /** simple version of for_each (no reduction) to be called for each index in the partition
     * @param f - must be a unary function accepting as input IndexType
     */
    template <class TUnaryFunction>
    inline void for_each(TUnaryFunction &&f)
    {
        #pragma omp parallel for
        for (int i=0; i<mNchunks; ++i) {
            for (auto k = mBlockPartition[i]; k < mBlockPartition[i+1]; ++k) {
                f(k); //note that we pass a reference to the value, not the iterator
            }
        }
    }

    /** version with reduction to be called for each index in the partition
     * function f is expected to return the values to be reduced
     * @param TReducer - template parameter specifying the type of reducer to be applied
     * @param f - must be a unary function accepting as input IndexType
     */
    template <class TReducer, class TUnaryFunction>
    inline typename TReducer::value_type for_each(TUnaryFunction &&f)
    {
        TReducer global_reducer;
        #pragma omp parallel for
        for (int i=0; i<mNchunks; ++i) {
            TReducer local_reducer;
            for (auto k = mBlockPartition[i]; k < mBlockPartition[i+1]; ++k) {
                local_reducer.LocalReduce(f(k));
            }
            global_reducer.ThreadSafeReduce(local_reducer);
        }
        return global_reducer.GetValue();
    }

private:
    int mNchunks;
    std::array<TIndexType, TMaxThreads> mBlockPartition;
};

} // namespace Kratos.

#endif // KRATOS_PARALLEL_UTILITIES_H_INCLUDED defined
