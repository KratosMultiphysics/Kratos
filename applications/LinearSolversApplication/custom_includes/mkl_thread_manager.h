// KRATOS  _     _                       ____        _
//        | |   (_)_ __   ___  __ _ _ __/ ___|  ___ | |_   _____ _ __ ___
//        | |   | | '_ \ / _ \/ _` | '__\___ \ / _ \| \ \ / / _ \ '__/ __|
//        | |___| | | | |  __/ (_| | |   ___) | (_) | |\ V /  __/ |  \__ |
//        |_____|_|_| |_|\___|\__,_|_|  |____/ \___/|_| \_/ \___|_|  |___/ Application
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes
#ifdef USE_EIGEN_MKL
#include "mkl.h" // MKL_INT, dcsrilu0, dcsrilut
#endif

// Project includes
#include "includes/thread_manager.h"

namespace Kratos
{
///@name Kratos Globals
///@{

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

/**
 * @class MKLThreadManager
 * @brief Thread manager implementation using Intel MKL for parallel operations in Kratos.
 * @details The MKLThreadManager class derives from ThreadManager and provides
 * concrete implementations for managing threads using the Intel Math Kernel Library (MKL).
 * It allows querying and setting the number of threads used for parallel execution via MKL,
 * as well as retrieving the number of available processing units as seen by MKL.
 * This class is intended to be used when Intel MKL is the underlying threading model.
 * @see ThreadManager
 */
class KRATOS_API(KRATOS_CORE) MKLThreadManager
    : public ThreadManager
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MKLThreadManager
    KRATOS_CLASS_POINTER_DEFINITION(MKLThreadManager);

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Returns the number of threads currently being used by MKL.
     * @details This method overrides the base class implementation and returns the number of threads
     *          currently set for MKL parallel regions.
     * @return The number of MKL threads in use.
     */
    int GetNumThreads() const override
    {
    #ifdef USE_EIGEN_MKL
        const int nthreads = mkl_get_max_threads();
        KRATOS_DEBUG_ERROR_IF(nthreads <= 0) << "GetNumThreads would devolve nthreads = " << nthreads << " which is not possible" << std::endl;
        return nthreads;
    #else
        KRATOS_ERROR << "Calling MKLThreadManager::GetNumThreads when OpenMP is not enabled!" << std::endl;
        return 1;
    #endif
    }

    /**
     * @brief Returns the number of processing units available to MKL.
     * @details This method overrides the base class implementation and returns the number of logical
     *          processors (cores or hardware threads) available for MKL to use.
     * @return The number of available processing units.
     */
    int GetNumProcs() const override
    {
    #ifdef USE_EIGEN_MKL
        const int num_procs = mkl_get_max_threads();
        KRATOS_DEBUG_ERROR_IF(num_procs <= 0) << "GetNumProcs would devolve num_procs = " << num_procs << " which is not possible" << std::endl;
        return num_procs;
    #else
        KRATOS_ERROR << "Calling MKLThreadManager::GetNumProcs when OpenMP is not enabled!" << std::endl;
        return 1;
    #endif
    }

    /**
     * @brief Sets the number of threads to be used by MKL.
     * @details This method overrides the base class implementation and sets the number of threads
     *          that MKL will use for parallel regions.
     * @param NumThreads The number of threads to set for MKL.
     */
    void SetNumThreads(const int NumThreads) override
    {
    #ifdef USE_EIGEN_MKL
        mkl_set_num_threads(NumThreads);
    #else
        KRATOS_ERROR << "Calling MKLThreadManager::SetNumThreads when OpenMP is not enabled!" << std::endl;
    #endif
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "MKLThreadManager" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MKLThreadManager";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {

    }

    ///@}
}; /// class MKLThreadManager

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  MKLThreadManager& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MKLThreadManager& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}
} /// namespace Kratos
