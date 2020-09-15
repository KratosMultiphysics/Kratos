//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, https://github.com/philbucher
//                   Suneth Warnakulasuriya, https://github.com/sunethwarna
//


#if !defined(KRATOS_SCOPE_LOCK_H_INCLUDED)
#define KRATOS_SCOPE_LOCK_H_INCLUDED


// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/lock_object.h"


namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/// Class for locking a scope.
/** This class is used for locking a scope, i.e. to make this scope act as a critical section.
*/
class KRATOS_API(KRATOS_CORE) ScopeLock
{
public:
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ScopeLock();

    /// Destructor.
    ~ScopeLock();

    ScopeLock(const ScopeLock&) = delete;
    ScopeLock& operator=(const ScopeLock&) = delete;

    ///@}
    ///@name Operations
    ///@{

    ///@}

private:
    ///@name Member Variables
    ///@{

    static LockObject mLock;

    ///@}

}; // Class ScopeLock

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SCOPE_LOCK_H_INCLUDED defined
