//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Vicente Mataix Ferrandiz
//

#pragma once

// System includes
#include <memory>
#include <utility>
#include <vector>

// External includes

// Project includes

namespace Kratos
{
///@addtogroup Kratos Core
///@{

///@name Kratos Classes
///@{

/**
 * @class DataCommunicatorRequest
 * @ingroup KratosCore
 * @brief Handle to an asynchronous (non-blocking) communication started by a DataCommunicator.
 * @details This is the MPI-agnostic counterpart of an MPI_Request. It is returned by the non-blocking
 * methods of DataCommunicator (ISend, IRecv, IBarrier, IBroadcast, ISumAll, ...). The operation is
 * only guaranteed to be finished (and the involved buffers safe to reuse or read) once Wait() returned
 * or Test() returned true.
 * A default-constructed request represents an already completed operation (this is what the serial
 * DataCommunicator returns).
 * If a request is destroyed while its operation is still pending, the destructor waits for its completion,
 * so that the communication buffers are never left in use by an operation nobody can complete.
 * The class is move-only, since it owns the underlying operation.
 */
class DataCommunicatorRequest
{
public:
    ///@name Type Definitions
    ///@{

    /// Interface to be implemented by each communication backend (e.g. MPI).
    class Implementation
    {
    public:
        /// Destructor. Implementations must complete a still pending operation here.
        virtual ~Implementation() = default;

        /// Block until the operation completes (wrapper for MPI_Wait).
        virtual void Wait() = 0;

        /// Check if the operation completed, without blocking (wrapper for MPI_Test).
        virtual bool Test() = 0;
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor. Creates an already completed request.
    DataCommunicatorRequest() = default;

    /// Constructor taking ownership of a backend-specific pending operation.
    explicit DataCommunicatorRequest(std::unique_ptr<Implementation> pImplementation)
        : mpImplementation(std::move(pImplementation))
    {
    }

    /// Move constructor.
    DataCommunicatorRequest(DataCommunicatorRequest&& rOther) = default;

    /// Move assignment. A pending operation held by this request is completed first.
    DataCommunicatorRequest& operator=(DataCommunicatorRequest&& rOther) = default;

    /// Copy constructor (deleted, requests own their operation).
    DataCommunicatorRequest(const DataCommunicatorRequest& rOther) = delete;

    /// Copy assignment (deleted, requests own their operation).
    DataCommunicatorRequest& operator=(const DataCommunicatorRequest& rOther) = delete;

    /// Destructor. Waits for a still pending operation.
    ~DataCommunicatorRequest() = default;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Block until the operation completes.
     * @details This is a wrapper for MPI_Wait. Calling it on a completed request does nothing.
     */
    void Wait()
    {
        if (mpImplementation) {
            mpImplementation->Wait();
            mpImplementation.reset();
        }
    }

    /**
     * @brief Check whether the operation has completed, without blocking.
     * @details This is a wrapper for MPI_Test.
     * @return True if the operation is completed, false otherwise.
     */
    bool Test()
    {
        if (mpImplementation && mpImplementation->Test()) {
            mpImplementation.reset();
        }
        return IsCompleted();
    }

    /**
     * @brief Block until all the given operations complete.
     * @details Equivalent to MPI_Waitall. Since all operations were already started,
     * waiting for them one after the other cannot deadlock.
     * @param rRequests The requests to wait for.
     */
    static void WaitAll(std::vector<DataCommunicatorRequest>& rRequests)
    {
        for (auto& r_request : rRequests) {
            r_request.Wait();
        }
    }

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief Check whether the operation is known to be completed.
     * @details This does not query the backend: use Test() to make progress on a pending operation.
     * @return True if Wait() was called or Test() returned true (or the request was never pending).
     */
    bool IsCompleted() const
    {
        return !mpImplementation;
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    std::unique_ptr<Implementation> mpImplementation;

    ///@}

}; // Class DataCommunicatorRequest

///@}

///@} addtogroup block

} // namespace Kratos.
