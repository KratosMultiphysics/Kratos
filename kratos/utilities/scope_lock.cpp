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


// System includes

// External includes

// Project includes
#include "utilities/scope_lock.h"


namespace Kratos {

ScopeLock::ScopeLock()
{
    mLock.SetLock();
}

ScopeLock::~ScopeLock()
{
    mLock.UnSetLock();
}

LockObject ScopeLock::mLock;

}  // namespace Kratos.
