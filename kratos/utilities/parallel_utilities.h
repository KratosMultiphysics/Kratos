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
//                   Philipp Bucher (https://github.com/philbucher)
//                   Vicente Mataix Ferrandiz
//

#pragma once

// System includes
#include <iostream>
#include <array>
#include <iterator>
#include <vector>
#include <tuple>
#include <cmath>
#include <limits>
#include <future>
#include <thread>
#include <mutex>

// External includes
#ifdef KRATOS_SMP_OPENMP
#include <omp.h>
#endif

// Project includes
#include "includes/define.h"
#include "includes/global_variables.h"
#include "includes/lock_object.h"

#define KRATOS_CRITICAL_SECTION const std::lock_guard scope_lock(ParallelUtilities::GetGlobalLock());

#define KRATOS_PREPARE_CATCH_THREAD_EXCEPTION std::stringstream err_stream;

#ifndef KRATOS_NO_TRY_CATCH
    #define KRATOS_CATCH_THREAD_EXCEPTION \
    } catch(Exception& e) { \
        KRATOS_CRITICAL_SECTION \
        err_stream << "Thread #" << i << " caught exception: " << e.what(); \
    } catch(std::exception& e) { \
        KRATOS_CRITICAL_SECTION \
        err_stream << "Thread #" << i << " caught exception: " << e.what(); \
    } catch(...) { \
        KRATOS_CRITICAL_SECTION \
        err_stream << "Thread #" << i << " caught unknown exception:"; \
    }
#else
    #define KRATOS_CATCH_THREAD_EXCEPTION {}
#endif

#define KRATOS_CHECK_AND_THROW_THREAD_EXCEPTION \
const std::string& err_msg = err_stream.str(); \
KRATOS_ERROR_IF_NOT(err_msg.empty()) << "The following errors occured in a parallel region!\n" << err_msg << std::endl;

namespace Kratos
{
///@addtogroup KratosCore

/// Shared memory parallelism related helper class
/** Provides access to functionalities for shared memory parallelism
 * such as the number of threads in usa.
*/
class KRATOS_API(KRATOS_CORE) ParallelUtilities
{
public:
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /** @brief Returns the current number of threads
     * @return number of threads
     */
    [[nodiscard]] static int GetNumThreads();

    /** @brief Sets the current number of threads
     * @param NumThreads - the number of threads to be used
     */
    static void SetNumThreads(const int NumThreads);

    /** @brief Returns the number of processors available to this device
     * This can include the multiple threads per processing unit
     * @return number of processors
     */
    [[nodiscard]] static int GetNumProcs();

    ///@}

    /** @brief Returns the global lock
     * Global lock that can be used for critical sections
     * @return global lock
     */
    [[nodiscard]] static LockObject& GetGlobalLock();

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    static LockObject* mspGlobalLock;

    static int* mspNumThreads;

    ///@}
    ///@name Private Operations
    ///@{

    /// Default constructor.
    ParallelUtilities() = delete;

    /** @brief Initializes the number of threads to be used.
     * @return number of threads
     */
    static int InitializeNumberOfThreads();
    ///@}

    ///@name Private Access
    ///@{

    static int& GetNumberOfThreads();

    ///@}
}; // Class ParallelUtilities


//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
/** @param TIterator - type of iterator (must be a random access iterator)
 *  @param MaxThreads - maximum number of threads allowed in the partitioning.
 *                       must be known at compile time to avoid heap allocations in the partitioning
 */
template<class TIterator, int MaxThreads=Globals::MaxAllowedThreads>
class BlockPartition
{
public:
    /** @param it_begin - iterator pointing at the beginning of the container
     *  @param it_end - iterator pointing to the end of the container
     *  @param Nchunks - number of threads to be used in the loop (must be lower than TMaxThreads)
     */
    BlockPartition(TIterator it_begin,
                   TIterator it_end,
                   int Nchunks = ParallelUtilities::GetNumThreads())
    {
        static_assert(
            std::is_same_v<typename std::iterator_traits<TIterator>::iterator_category, std::random_access_iterator_tag>,
            "BlockPartition requires random access iterators to divide the input range into partitions"
        );
        KRATOS_ERROR_IF(Nchunks < 1) << "Number of chunks must be > 0 (and not " << Nchunks << ")" << std::endl;

        const std::ptrdiff_t size_container = it_end-it_begin;

        if (size_container == 0) {
            mNchunks = Nchunks;
        } else {
            // in case the container is smaller than the number of chunks
            mNchunks = std::min(static_cast<int>(size_container), Nchunks);
        }
        const std::ptrdiff_t block_partition_size = size_container / mNchunks;
        mBlockPartition[0] = it_begin;
        mBlockPartition[mNchunks] = it_end;
        for (int i=1; i<mNchunks; i++) {
            mBlockPartition[i] = mBlockPartition[i-1] + block_partition_size;
        }
    }

    /** @brief simple iteration loop. f called on every entry in rData
     * @param f - must be a unary function accepting as input TContainerType::value_type&
     */
    template <class TUnaryFunction>
    inline void for_each(TUnaryFunction&& f)
    {
        KRATOS_PREPARE_CATCH_THREAD_EXCEPTION

        #pragma omp parallel for schedule(runtime)
        for (int i=0; i<mNchunks; ++i) {
            KRATOS_TRY
            for (auto it = mBlockPartition[i]; it != mBlockPartition[i+1]; ++it) {
                f(*it); //note that we pass the value to the function, not the iterator
            }
            KRATOS_CATCH_THREAD_EXCEPTION
        }

        KRATOS_CHECK_AND_THROW_THREAD_EXCEPTION
    }

    /** @brief loop allowing reductions. f called on every entry in rData
     * the function f needs to return the values to be used by the reducer
     * @param TReducer template parameter specifying the reduction operation to be done
     * @param f - must be a unary function accepting as input TContainerType::value_type&
     */
    template <class TReducer, class TUnaryFunction>
    [[nodiscard]] inline typename TReducer::return_type for_each(TUnaryFunction &&f)
    {
        KRATOS_PREPARE_CATCH_THREAD_EXCEPTION

        TReducer global_reducer;
        #pragma omp parallel for schedule(runtime)
        for (int i=0; i<mNchunks; ++i) {
            KRATOS_TRY
            TReducer local_reducer;
            for (auto it = mBlockPartition[i]; it != mBlockPartition[i+1]; ++it) {
                local_reducer.LocalReduce(f(*it));
            }
            global_reducer.ThreadSafeReduce(local_reducer);
            KRATOS_CATCH_THREAD_EXCEPTION
        }

        KRATOS_CHECK_AND_THROW_THREAD_EXCEPTION

        return global_reducer.GetValue();
    }

    /** @brief loop with thread local storage (TLS). f called on every entry in rData
     * @param TThreadLocalStorage template parameter specifying the thread local storage
     * @param f - must be a function accepting as input TContainerType::value_type& and the thread local storage
     */
    template <class TThreadLocalStorage, class TFunction>
    inline void for_each(const TThreadLocalStorage& rThreadLocalStoragePrototype, TFunction &&f)
    {
        static_assert(std::is_copy_constructible<TThreadLocalStorage>::value, "TThreadLocalStorage must be copy constructible!");

        KRATOS_PREPARE_CATCH_THREAD_EXCEPTION

        #pragma omp parallel
        {
            // copy the prototype to create the thread local storage
            TThreadLocalStorage thread_local_storage(rThreadLocalStoragePrototype);

            #pragma omp for
            for(int i=0; i<mNchunks; ++i){
                KRATOS_TRY
                for (auto it = mBlockPartition[i]; it != mBlockPartition[i+1]; ++it){
                    f(*it, thread_local_storage); // note that we pass the value to the function, not the iterator
                }
                KRATOS_CATCH_THREAD_EXCEPTION
            }
        }
        KRATOS_CHECK_AND_THROW_THREAD_EXCEPTION
    }

    /** @brief loop with thread local storage (TLS) allowing reductions. f called on every entry in rData
     * the function f needs to return the values to be used by the reducer
     * @param TReducer template parameter specifying the reduction operation to be done
     * @param TThreadLocalStorage template parameter specifying the thread local storage
     * @param f - must be a function accepting as input TContainerType::value_type& and the thread local storage
     */
    template <class TReducer, class TThreadLocalStorage, class TFunction>
    [[nodiscard]] inline typename TReducer::return_type for_each(const TThreadLocalStorage& rThreadLocalStoragePrototype, TFunction &&f)
    {
        static_assert(std::is_copy_constructible<TThreadLocalStorage>::value, "TThreadLocalStorage must be copy constructible!");

        KRATOS_PREPARE_CATCH_THREAD_EXCEPTION

        TReducer global_reducer;

        #pragma omp parallel
        {
            // copy the prototype to create the thread local storage
            TThreadLocalStorage thread_local_storage(rThreadLocalStoragePrototype);

            #pragma omp for
            for (int i=0; i<mNchunks; ++i) {
                KRATOS_TRY
                TReducer local_reducer;
                for (auto it = mBlockPartition[i]; it != mBlockPartition[i+1]; ++it) {
                    local_reducer.LocalReduce(f(*it, thread_local_storage));
                }
                global_reducer.ThreadSafeReduce(local_reducer);
                KRATOS_CATCH_THREAD_EXCEPTION
            }
        }
        KRATOS_CHECK_AND_THROW_THREAD_EXCEPTION
        return global_reducer.GetValue();
    }

private:
    int mNchunks;
    std::array<TIterator, MaxThreads> mBlockPartition;
};

/** @brief Execute a functor on all items of a range in parallel.
 *  @tparam TIterator: random access iterator.
 *  @tparam TFunction: functor taking the dereferenced type of @a TIterator.
 *  @param itBegin: iterator to the first item in the container to loop on.
 *  @param itEnd: iterator past the last item in the container.
 *  @param rFunction: function to execute on each item.
 */
template <class TIterator,
          class TFunction,
          std::enable_if_t<std::is_base_of_v<std::input_iterator_tag, typename std::iterator_traits<TIterator>::iterator_category>,bool> = true>
void block_for_each(TIterator itBegin, TIterator itEnd, TFunction&& rFunction)
{
    BlockPartition<TIterator>(itBegin, itEnd).for_each(std::forward<TFunction>(rFunction));
}

/** @brief Execute a functor on all items of a range in parallel, and perform a reduction.
 *  @tparam TReduction: type of reduction to apply. See @ref SumReduction for an example.
 *  @tparam TIterator: random access iterator.
 *  @tparam TFunction: functor taking the dereferenced type of @a TIterator.
 *  @param itBegin: iterator to the first item in the container to loop on.
 *  @param itEnd: iterator past the last item in the container.
 *  @param rFunction: function to execute on each item.
 */
template <class TReduction,
          class TIterator,
          class TFunction,
          std::enable_if_t<std::is_base_of_v<std::input_iterator_tag, typename std::iterator_traits<TIterator>::iterator_category>,bool> = true>
[[nodiscard]] typename TReduction::return_type block_for_each(TIterator itBegin, TIterator itEnd, TFunction&& rFunction)
{
    return  BlockPartition<TIterator>(itBegin, itEnd).template for_each<TReduction>(std::forward<TFunction>(std::forward<TFunction>(rFunction)));
}

/** @brief Execute a functor with thread local storage on all items of a range in parallel.
 *  @tparam TIterator:  random access iterator.
 *  @tparam TTLS: copy constructible thread-local type.
 *  @tparam TFunction: functor taking the dereferenced type of @a TIterator.
 *  @param itBegin: iterator to the first item in the container to loop on.
 *  @param itEnd: iterator past the last item in the container.
 *  @param rTLS: thread local storage
 *  @param rFunction: function to execute on each item.
 */
template <class TIterator,
          class TTLS,
          class TFunction,
          std::enable_if_t<std::is_base_of_v<std::input_iterator_tag, typename std::iterator_traits<TIterator>::iterator_category>,bool> = true>
void block_for_each(TIterator itBegin, TIterator itEnd, const TTLS& rTLS, TFunction &&rFunction)
{
     BlockPartition<TIterator>(itBegin, itEnd).for_each(rTLS, std::forward<TFunction>(rFunction));
}

/** @brief Execute a functor with thread local storage on all items of a range in parallel, and perform a reduction.
 *  @tparam TReduction: type of reduction to apply. See @ref SumReduction for an example.
 *  @tparam TIterator: random access iterator.
 *  @tparam TTLS: copy constructible thread-local type.
 *  @tparam TFunction: functor taking the dereferenced type of @a TIterator.
 *  @param itBegin: iterator to the first item in the container to loop on.
 *  @param itEnd: iterator past the last item in the container.
 *  @param rTLS: thread local storage
 *  @param rFunction: function to execute on each item.
 */
template <class TReduction,
          class TIterator,
          class TTLS,
          class TFunction,
          std::enable_if_t<std::is_base_of_v<std::input_iterator_tag, typename std::iterator_traits<TIterator>::iterator_category>,bool> = true>
[[nodiscard]] typename TReduction::return_type block_for_each(TIterator itBegin, TIterator itEnd, const TTLS& tls, TFunction&& rFunction)
{
    return BlockPartition<TIterator>(itBegin, itEnd).template for_each<TReduction>(tls, std::forward<TFunction>(std::forward<TFunction>(rFunction)));
}

/** @brief simplified version of the basic loop (without reduction) to enable template type deduction
 *  @tparam TContainerType A standard-conforming container type.
 *  @tparam TFunctionType Functor operating on @a TContainerType::value_type.
 *  @param v - containers to be looped upon
 *  @param func - must be a unary function accepting as input TContainerType::value_type&
 */
template <class TContainerType,
          class TFunctionType,
          std::enable_if_t<!std::is_same_v<
            std::iterator_traits<typename decltype(std::declval<std::remove_cv_t<TContainerType>>().begin())::value_type>,
            void
          >, bool> = true
         >
void block_for_each(TContainerType &&v, TFunctionType &&func)
{
    block_for_each(v.begin(), v.end(), std::forward<TFunctionType>(func));
}

/** @brief simplified version of the basic loop with reduction to enable template type deduction
 *  @tparam TReducer Reduction type to apply. See @ref SumReduction as an example.
 *  @tparam TContainerType A standard-conforming container type.
 *  @tparam TFunctionType Functor operating on @a TContainerType::value_type.
 *  @param v - containers to be looped upon
 *  @param func - must be a unary function accepting as input TContainerType::value_type&
 */
template <class TReducer,
          class TContainerType,
          class TFunctionType,
          std::enable_if_t<!std::is_same_v<
            std::iterator_traits<typename decltype(std::declval<std::remove_cv_t<TContainerType>>().begin())::value_type>,
            void
          >, bool> = true
         >
[[nodiscard]] typename TReducer::return_type block_for_each(TContainerType &&v, TFunctionType &&func)
{
    return block_for_each<TReducer>(v.begin(), v.end(), std::forward<TFunctionType>(func));
}

/** @brief simplified version of the basic loop with thread local storage (TLS) to enable template type deduction
 *  @tparam TContainerType A standard-conforming container type.
 *  @tparam TThreadLocalStorage Copy constructible thread-local type.
 *  @tparam TFunctionType Functor operating on @a TContainerType::value_type.
 *  @param v - containers to be looped upon
 *  @param tls - thread local storage
 *  @param func - must be a function accepting as input TContainerType::value_type& and the thread local storage
 */
template <class TContainerType,
          class TThreadLocalStorage,
          class TFunctionType,
          std::enable_if_t<!std::is_same_v<
            std::iterator_traits<typename decltype(std::declval<std::remove_cv_t<TContainerType>>().begin())::value_type>,
            void
          >, bool> = true
         >
void block_for_each(TContainerType &&v, const TThreadLocalStorage& tls, TFunctionType &&func)
{
    block_for_each(v.begin(), v.end(), tls, std::forward<TFunctionType>(func));
}

/** @brief simplified version of the basic loop with reduction and thread local storage (TLS) to enable template type deduction
 *  @tparam TReducer Reduction type to apply. See @ref SumReduction as an example.
 *  @tparam TContainerType A standard-conforming container type.
 *  @tparam TThreadLocalStorage Copy constructible thread-local type.
 *  @tparam TFunctionType Functor operating on @a TContainerType::value_type.
 *  @param v - containers to be looped upon
 *  @param tls - thread local storage
 *  @param func - must be a function accepting as input TContainerType::value_type& and the thread local storage
 */
template <class TReducer,
          class TContainerType,
          class TThreadLocalStorage,
          class TFunctionType,
          std::enable_if_t<!std::is_same_v<
            std::iterator_traits<typename decltype(std::declval<std::remove_cv_t<TContainerType>>().begin())::value_type>,
            void
          >, bool> = true
         >
[[nodiscard]] typename TReducer::return_type block_for_each(TContainerType &&v, const TThreadLocalStorage& tls, TFunctionType &&func)
{
    return block_for_each<TReducer>(v.begin(), v.end(), tls, std::forward<TFunctionType>(func));
}

//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
/** @brief This class is useful for index iteration over containers
 *  @param TIndexType type of index to be used in the loop
 *  @param TMaxThreads - maximum number of threads allowed in the partitioning.
 *                       must be known at compile time to avoid heap allocations in the partitioning
 */
template<class TIndexType=std::size_t, int TMaxThreads=Globals::MaxAllowedThreads>
class IndexPartition
{
public:

/** @brief constructor using the size of the partition to be used
 *  @param Size - the size of the partition
 *  @param Nchunks - number of threads to be used in the loop (must be lower than TMaxThreads)
 */
    IndexPartition(TIndexType Size,
                   int Nchunks = ParallelUtilities::GetNumThreads())
    {
        KRATOS_ERROR_IF(Nchunks < 1) << "Number of chunks must be > 0 (and not " << Nchunks << ")" << std::endl;

        if (Size == 0) {
            mNchunks = Nchunks;
        } else {
            // in case the container is smaller than the number of chunks
            mNchunks = std::min(static_cast<int>(Size), Nchunks);
        }

        const int block_partition_size = Size / mNchunks;
        mBlockPartition[0] = 0;
        mBlockPartition[mNchunks] = Size;
        for (int i=1; i<mNchunks; i++) {
            mBlockPartition[i] = mBlockPartition[i-1] + block_partition_size;
        }

    }

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
        KRATOS_PREPARE_CATCH_THREAD_EXCEPTION

        #pragma omp parallel for schedule(runtime)
        for (int i=0; i<mNchunks; ++i) {
            KRATOS_TRY
            for (auto k = mBlockPartition[i]; k < mBlockPartition[i+1]; ++k) {
                f(k); //note that we pass a reference to the value, not the iterator
            }
            KRATOS_CATCH_THREAD_EXCEPTION
        }
        KRATOS_CHECK_AND_THROW_THREAD_EXCEPTION
    }

    /** version with reduction to be called for each index in the partition
     * function f is expected to return the values to be reduced
     * @param TReducer - template parameter specifying the type of reducer to be applied
     * @param f - must be a unary function accepting as input IndexType
     */
    template <class TReducer, class TUnaryFunction>
    [[nodiscard]] inline typename TReducer::return_type for_each(TUnaryFunction &&f)
    {
        KRATOS_PREPARE_CATCH_THREAD_EXCEPTION

        TReducer global_reducer;
        #pragma omp parallel for schedule(runtime)
        for (int i=0; i<mNchunks; ++i) {
            KRATOS_TRY
            TReducer local_reducer;
            for (auto k = mBlockPartition[i]; k < mBlockPartition[i+1]; ++k) {
                local_reducer.LocalReduce(f(k));
            }
            global_reducer.ThreadSafeReduce(local_reducer);
            KRATOS_CATCH_THREAD_EXCEPTION
        }
        KRATOS_CHECK_AND_THROW_THREAD_EXCEPTION
        return global_reducer.GetValue();
    }


    /** @brief loop with thread local storage (TLS). f called on every entry in rData
     * @param TThreadLocalStorage template parameter specifying the thread local storage
     * @param f - must be a function accepting as input IndexType and the thread local storage
     */
    template <class TThreadLocalStorage, class TFunction>
    inline void for_each(const TThreadLocalStorage& rThreadLocalStoragePrototype, TFunction &&f)
    {
        static_assert(std::is_copy_constructible<TThreadLocalStorage>::value, "TThreadLocalStorage must be copy constructible!");

        KRATOS_PREPARE_CATCH_THREAD_EXCEPTION

        #pragma omp parallel
        {
            // copy the prototype to create the thread local storage
            TThreadLocalStorage thread_local_storage(rThreadLocalStoragePrototype);

            #pragma omp for
            for (int i=0; i<mNchunks; ++i) {
                KRATOS_TRY
                for (auto k = mBlockPartition[i]; k < mBlockPartition[i+1]; ++k) {
                    f(k, thread_local_storage); //note that we pass a reference to the value, not the iterator
                }
                KRATOS_CATCH_THREAD_EXCEPTION
            }
        }
        KRATOS_CHECK_AND_THROW_THREAD_EXCEPTION
    }

    /** version with reduction and thread local storage (TLS) to be called for each index in the partition
     * function f is expected to return the values to be reduced
     * @param TReducer - template parameter specifying the type of reducer to be applied
     * @param TThreadLocalStorage template parameter specifying the thread local storage
     * @param f - must be a function accepting as input IndexType and the thread local storage
     */
    template <class TReducer, class TThreadLocalStorage, class TFunction>
    [[nodiscard]] inline typename TReducer::return_type for_each(const TThreadLocalStorage& rThreadLocalStoragePrototype, TFunction &&f)
    {
        static_assert(std::is_copy_constructible<TThreadLocalStorage>::value, "TThreadLocalStorage must be copy constructible!");

        KRATOS_PREPARE_CATCH_THREAD_EXCEPTION

        TReducer global_reducer;

        #pragma omp parallel
        {
            // copy the prototype to create the thread local storage
            TThreadLocalStorage thread_local_storage(rThreadLocalStoragePrototype);

            #pragma omp for
            for (int i=0; i<mNchunks; ++i) {
                KRATOS_TRY
                TReducer local_reducer;
                for (auto k = mBlockPartition[i]; k < mBlockPartition[i+1]; ++k) {
                    local_reducer.LocalReduce(f(k, thread_local_storage));
                }
                global_reducer.ThreadSafeReduce(local_reducer);
                KRATOS_CATCH_THREAD_EXCEPTION
            }
        }
        KRATOS_CHECK_AND_THROW_THREAD_EXCEPTION

        return global_reducer.GetValue();
    }

private:
    int mNchunks;
    std::array<TIndexType, TMaxThreads> mBlockPartition;
};

} // namespace Kratos.

#undef KRATOS_PREPARE_CATCH_THREAD_EXCEPTION
#undef KRATOS_CATCH_THREAD_EXCEPTION
#undef KRATOS_CHECK_AND_THROW_THREAD_EXCEPTION
