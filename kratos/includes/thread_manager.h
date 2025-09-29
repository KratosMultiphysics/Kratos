//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/smart_pointers.h"

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
 * @class ThreadManager
 * @brief Base class for managing thread-related operations in Kratos.
 * @details The ThreadManager class provides an interface for controlling and configuring
 * the number of threads used in parallel operations. It is intended to be used
 * as a base class for platform-specific or parallelization-library-specific
 * thread managers. Derived classes should implement the thread management logic
 * appropriate for the underlying threading model (e.g. OMP, MKL).
 * @note The base implementation of SetNumThreads will throw an error if called.
 *       Derived classes must override this method to provide the actual
 *       functionality.
 */
class KRATOS_API(KRATOS_CORE) ThreadManager
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ThreadManager
    KRATOS_CLASS_POINTER_DEFINITION(ThreadManager);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ThreadManager() = default;

    /// Destructor.
    ~ThreadManager() = default;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Returns the number of threads currently being used.
     * @details This is a virtual method intended to be overridden by derived classes.
     *          The base implementation throws an error if called.
     * @return The number of threads (default implementation returns 1).
     */
    virtual int GetNumThreads() const;

    /**
     * @brief Returns the number of processing units available.
     * @details This is a virtual method intended to be overridden by derived classes.
     *          The base implementation throws an error if called.
     * @return The number of processing units (default implementation returns 1).
     */
    virtual int GetNumProcs() const;

    /**
     * @brief Sets the current number of threads to be used.
     * @details This is a virtual method intended to be overridden by derived classes.
     *          The base implementation throws an error if called.
     * @param NumThreads The number of threads to set.
     */
    virtual void SetNumThreads(const int NumThreads);

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "ThreadManager" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ThreadManager";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {

    }

    ///@}
}; /// class ThreadManager

/**
 * @class CXX11ThreadManager
 * @brief Thread manager implementation using C++11 threads for parallel operations in Kratos.
 * @details The CXX11ThreadManager class derives from ThreadManager and provides
 * concrete implementations for managing threads using the C++11 threading library.
 * It allows querying and setting the number of threads used for parallel execution,
 * as well as retrieving the number of available processing units.
 * This class is intended to be used when C++11 is the underlying threading model.
 * @see ThreadManager
 */
class KRATOS_API(KRATOS_CORE) CXX11ThreadManager
    : public ThreadManager
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CXX11ThreadManager
    KRATOS_CLASS_POINTER_DEFINITION(CXX11ThreadManager);

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Returns the number of threads currently being used by C++11.
     * @details This method overrides the base class implementation and returns the number of threads
     *          currently set for C++11 parallel regions.
     * @return The number of C++11 threads in use.
     */
    int GetNumThreads() const override;

    /**
     * @brief Returns the number of processing units available to C++11.
     * @details This method overrides the base class implementation and returns the number of logical
     *          processors (cores or hardware threads) available for C++11 to use.
     * @return The number of available processing units.
     */
    int GetNumProcs() const override;

    /**
     * @brief Sets the number of threads to be used by C++11.
     * @details This method overrides the base class implementation and sets the number of threads
     *          that C++11 will use for parallel regions.
     * @param NumThreads The number of threads to set for C++11.
     */
    void SetNumThreads(const int NumThreads) override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "CXX11ThreadManager" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "CXX11ThreadManager";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {

    }

    ///@}
}; /// class CXX11ThreadManager

/**
 * @class OMPThreadManager
 * @brief Thread manager implementation using OpenMP for parallel operations in Kratos.
 * @details The OMPThreadManager class derives from ThreadManager and provides
 * concrete implementations for managing threads using the OpenMP parallelization library.
 * It allows querying and setting the number of threads used for parallel execution,
 * as well as retrieving the number of available processing units.
 * This class is intended to be used when OpenMP is the underlying threading model.
 * @see ThreadManager
 */
class KRATOS_API(KRATOS_CORE) OMPThreadManager
    : public ThreadManager
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of OMPThreadManager
    KRATOS_CLASS_POINTER_DEFINITION(OMPThreadManager);

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Returns the number of threads currently being used by OpenMP.
     * @details This method overrides the base class implementation and returns the number of threads
     *          currently set for OpenMP parallel regions.
     * @return The number of OpenMP threads in use.
     */
    int GetNumThreads() const override;

    /**
     * @brief Returns the number of processing units available to OpenMP.
     * @details This method overrides the base class implementation and returns the number of logical
     *          processors (cores or hardware threads) available for OpenMP to use.
     * @return The number of available processing units.
     */
    int GetNumProcs() const override;

    /**
     * @brief Sets the number of threads to be used by OpenMP.
     * @details This method overrides the base class implementation and sets the number of threads
     *          that OpenMP will use for parallel regions.
     * @param NumThreads The number of threads to set for OpenMP.
     */
    void SetNumThreads(const int NumThreads) override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "OMPThreadManager" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "OMPThreadManager";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {

    }

    ///@}
}; /// class OMPThreadManager

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ThreadManager& rThis);

inline std::istream& operator >> (std::istream& rIStream,
                                  CXX11ThreadManager& rThis);

inline std::istream& operator >> (std::istream& rIStream,
                                  OMPThreadManager& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ThreadManager& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

inline std::ostream& operator << (std::ostream& rOStream,
                                  const CXX11ThreadManager& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

inline std::ostream& operator << (std::ostream& rOStream,
                                  const OMPThreadManager& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}
} /// namespace Kratos