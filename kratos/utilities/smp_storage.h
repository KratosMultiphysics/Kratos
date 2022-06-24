//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(KRATOS_SMP_STORAGE_H_INCLUDED )
#define  KRATOS_SMP_STORAGE_H_INCLUDED

// System includes
#ifdef KRATOS_SMP_OPENMP
#include <vector>
#include "utilities/openmp_utils.h"
#endif

// Project includes
#include "utilities/parallel_utilities.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

template<class TStorageType>
class SMPStorage
{
public:
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SMPStorage(){
        #ifdef KRATOS_SMP_OPENMP
            mThreadLocalStorageContainer.resize(ParallelUtilities::GetNumThreads());
        #endif
    }

    /// Destructor.
    ~SMPStorage() = default;

    /// Deleted copy constructor
    SMPStorage(SMPStorage const& rOther) = delete;

    /// Deleted assignment operator
    SMPStorage& operator=(SMPStorage const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    inline TStorageType& GetThreadLocalStorage()
    {
        #ifdef KRATOS_SMP_OPENMP
            return mThreadLocalStorageContainer[OpenMPUtils::ThisThread()];
        #elif defined(KRATOS_SMP_CXX11)
            thread_local TStorageType tls;
            return tls;
        #else
            return mThreadLocalStorageContainer;
        #endif
    }

    ///@}
private:
    ///@}
    ///@name Private Members
    ///@{

#ifdef KRATOS_SMP_OPENMP
    std::vector<TStorageType> mThreadLocalStorageContainer;
#elif defined(KRATOS_SMP_CXX11)
#else
    TStorageType mThreadLocalStorageContainer;
#endif

    ///@}

}; // Class SMPStorage

///@}

///@} addtogroup KratosCore

}  // namespace Kratos.

#endif // KRATOS_SMP_STORAGE_H_INCLUDED  defined
