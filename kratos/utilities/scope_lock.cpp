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
#include "includes/define.h"
#include "utilities/scope_lock.h"


namespace Kratos {

ScopeLock::ScopeLock()
{
    KRATOS_ERROR_IF_NOT(ScopeLock::mIsInitialized) << "ScopeLock is not initialized!" << std::endl;
    mpLock->SetLock();
}

ScopeLock::~ScopeLock()
{
    mpLock->UnSetLock();
}

void ScopeLock::InitializeLock()
{
    KRATOS_ERROR_IF(ScopeLock::mIsInitialized) << "ScopeLock was already initialized!" << std::endl;

    static LockObject static_lock;
    mpLock = &static_lock;
    mIsInitialized = true;
}

LockObject* ScopeLock::mpLock = nullptr;
bool ScopeLock::mIsInitialized = false;

}  // namespace Kratos.
