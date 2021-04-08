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
//
//

#if !defined(KRATOS_LOCK_OBJECT_H_INCLUDED)
#define KRATOS_LOCK_OBJECT_H_INCLUDED

// System includes

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
    LockObject() noexcept {
#ifdef KRATOS_SMP_OPENMP
			omp_init_lock(&mLock);
#endif
    }

    /// Copy constructor.
    LockObject(LockObject const& rOther) = delete;

    /// Move constructor.
    KRATOS_DEPRECATED_MESSAGE("The move constructor is deprecated and will be removed in the future!")
    LockObject(LockObject&& rOther) noexcept
#ifdef KRATOS_SMP_OPENMP
			: mLock(rOther.mLock)
#endif
    {
#ifdef KRATOS_SMP_OPENMP
        static_assert(std::is_move_constructible<omp_lock_t>::value, "omp_lock_t is not move constructible!");
			omp_init_lock(&mLock);
#endif
    }

    /// Destructor.
    virtual ~LockObject() noexcept
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
        //does nothing if openMP is not present
#ifdef KRATOS_SMP_OPENMP
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
        //does nothing if openMP is not present
#ifdef KRATOS_SMP_OPENMP
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
#ifdef KRATOS_SMP_OPENMP
        return omp_test_lock(&mLock);
#endif
        return true;
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

#ifdef KRATOS_SMP_OPENMP
	    mutable omp_lock_t mLock;
#endif

    ///@}

}; // Class LockObject

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_LOCK_OBJECT_H_INCLUDED defined
