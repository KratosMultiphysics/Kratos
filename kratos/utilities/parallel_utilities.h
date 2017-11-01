#pragma once

//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license:
//kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "utilities/openmp_utils.h"

namespace Kratos {
namespace Parallel {
///@addtogroup Kratos
///@{

/// type definitions
typedef omp_lock_t LockType;
typedef int ThreadIdType;  // do not expect this to be an int. It is not for
                           // C++11 parallelism

inline LockType ConstructLock() {
    LockType my_lock;
    omp_init_lock(&my_lock);
    return my_lock;
}

inline void SetLock(LockType &Lock) { omp_set_lock(&Lock); }

inline void UnSetLock(LockType &Lock) { omp_unset_lock(&Lock); }

inline int GetNumThreads() {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

inline ThreadIdType GetThreadHashId() {
#ifdef _OPENMP
    return omp_get_thread_num();
#else
    return 0;
#endif
}

///@name Kratos Globals
///@{
// this function implements a parallel for, in which f is a lambda function
// which captures its parameters once per very time it is called
template <class InputIt, class UnaryFunction>
inline void parallel_for(InputIt first, const int n, UnaryFunction &&f) {
#pragma omp parallel for firstprivate(first)
    for (int i = 0; i < n; ++i) {
        f(first + i);
    }
}

// parallel for executed by chunks. Here the function is called on every "block
// of unknowns", so that the lambda captures are done just once
// note that the function b must implement the loop from start to end
// work allocation is static
template <class TContainterType, class UnaryFunction>
inline void parallel_for(TContainterType &data, UnaryFunction &&f) {
    parallel_for(data.begin(), data.size(), f);
}

// the function f is executed in parallel. The parameter passed to f is the
// number of the thread executing it.
// one lambda capture is passed per every thread
template <class UnaryFunction>
inline void execute(UnaryFunction &&f) {
    const int NumThreads = OpenMPUtils::GetNumThreads();

#pragma omp parallel for firstprivate(NumThreads)
    for (int i = 0; i < NumThreads; ++i) {
        f(i);
    }
}

// parallel for executed by chunks. Here the function is called on every "block
// of unknowns", so that the lambda captures are done just once
// note that the function b must implement the loop from start to end
// work allocation is static
template <class InputIt, class BinaryFunction>
inline void block_parallel_for(InputIt first, const int n, BinaryFunction &&b) {
    const int NumThreads = OpenMPUtils::GetNumThreads();
    OpenMPUtils::PartitionVector Partitions;
    OpenMPUtils::DivideInPartitions(n, NumThreads, Partitions);

#pragma omp parallel for firstprivate(NumThreads)
    for (int i = 0; i < NumThreads; ++i) {
        b(first + Partitions[i], first + Partitions[i + 1]);
    }
}

// parallel for executed by chunks. Here the function is called on every "block
// of unknowns", so that the lambda captures are done just once
// note that the function b must implement the loop from start to end
// work allocation is dynamic
template <class InputIt, class BinaryFunction>
inline void block_parallel_for(
    InputIt first, const int n, const int NChunks, BinaryFunction &&b) {
    OpenMPUtils::PartitionVector Partitions;
    OpenMPUtils::DivideInPartitions(n, NChunks, Partitions);

#pragma omp parallel for firstprivate(first, NChunks) schedule(dynamic)
    for (int i = 0; i < NChunks; ++i) {
        b(first + Partitions[i], first + Partitions[i + 1]);
    }
}

template <class TContainterType, class BinaryFunction>
inline void block_parallel_for(TContainterType &data, BinaryFunction &&b) {
    Kratos::Parallel::block_parallel_for(data.begin(), data.size(), b);
}

template <class TContainterType, class BinaryFunction>
inline void block_parallel_for(
    TContainterType &data, const int NChunks, BinaryFunction &&b) {
    Kratos::Parallel::block_parallel_for(data.begin(), data.size(), NChunks, b);
}

inline void AtomicAdd(const double &value, double &target) {
#pragma omp atomic
    target += value;
}
inline void AtomicAdd(const int &value, int &target) {
#pragma omp atomic
    target += value;
}

template <class TVectorType1, class TVectorType2>
inline void AtomicAddVector(const TVectorType1 &value, TVectorType2 &target) {
    for (unsigned int i = 0; i < target.size(); i++)
        AtomicAdd(value[i], target[i]);
}

inline void AtomicSub(const double &value, double &target) {
#pragma omp atomic
    target -= value;
}
inline void AtomicSub(const int &value, int &target) {
#pragma omp atomic
    target -= value;
}

template <class TVectorType1, class TVectorType2>
inline void AtomicSubVector(const TVectorType1 &value, TVectorType2 &target) {
    for (unsigned int i = 0; i < target.size(); i++)
        AtomicSub(value[i], target[i]);
}

inline void AtomicAssign(const double &value, double &target) {
#pragma omp atomic write
    target = value;
}
inline void AtomicAssign(const int &value, int &target) {
#pragma omp atomic write
    target = value;
}

template <class TVectorType1, class TVectorType2>
inline void AtomicAssignVector(
    const TVectorType1 &value, TVectorType2 &target) {
    for (unsigned int i = 0; i < target.size(); i++)
        AtomicAssign(value[i], target[i]);
}
///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

}  // namespace Parallel

}  // namespace Kratos.
