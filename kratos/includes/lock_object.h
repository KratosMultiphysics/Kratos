//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Philipp Bucher (https://github.com/philbucher)
//

#pragma once

// System includes
#include <mutex>

// External includes
#ifdef KRATOS_SMP_OPENMP
#include <omp.h>
#endif

// Project includes
#include "includes/define.h"


namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/// This class defines and stores a lock and gives an interface to it.
/** The class makes a tiny wrapper over shared memory locking mechanisms
 * it is compliant with C++ Lockable
 * see https://en.cppreference.com/w/cpp/named_req/Lockable
 */
class LockObject
{
public:
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    LockObject() noexcept
    {
#ifdef KRATOS_SMP_OPENMP
        omp_init_lock(&mLock);
#endif
    }

    /// Copy constructor.
    LockObject(LockObject const& rOther) = delete;

    /// Destructor.
    ~LockObject() noexcept
    {
#ifdef KRATOS_SMP_OPENMP
        omp_destroy_lock(&mLock);
#endif
    }

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    LockObject& operator=(LockObject const& rOther) = delete;

    ///@}
    ///@name Access
    ///@{

    inline void lock() const
    {
#ifdef KRATOS_SMP_CXX11
        mLock.lock();
#elif KRATOS_SMP_OPENMP
        omp_set_lock(&mLock);
#endif
    }

    KRATOS_DEPRECATED_MESSAGE("Please use lock instead")
    inline void SetLock() const
    {
        this->lock();
    }

    inline void unlock() const
    {
#ifdef KRATOS_SMP_CXX11
        mLock.unlock();
#elif KRATOS_SMP_OPENMP
        omp_unset_lock(&mLock);
#endif
    }

    KRATOS_DEPRECATED_MESSAGE("Please use unlock instead")
    inline void UnSetLock() const
    {
        this->unlock();
    }

    inline bool try_lock() const
    {
#ifdef KRATOS_SMP_CXX11
        return mLock.try_lock();
#elif KRATOS_SMP_OPENMP
        return omp_test_lock(&mLock);
#endif
        return true;
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

#ifdef KRATOS_SMP_CXX11
        mutable std::mutex mLock;
#elif KRATOS_SMP_OPENMP
	    mutable omp_lock_t mLock;
#endif

    ///@}

}; // Class LockObject

///@}

///@} addtogroup block

}  // namespace Kratos.
