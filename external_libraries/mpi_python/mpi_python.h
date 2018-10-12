//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//                   Michael Andre
//                   Philipp Bucher
//

#ifndef KRATOS_MPI_PYTHON_H
#define KRATOS_MPI_PYTHON_H

// System includes
#include <vector>
#include <numeric>
#include "mpi.h"

// External includes

// Project includes


namespace Kratos {

/// A wrapper arround the MPI_Comm MPI communicator class.
/** This class is used to provide basic MPI communication
 * capabilities from Python. At its core, it stores an
 * MPI Communicator object, which can be used by the caller
 * to perform basic actions such as synchronizing the different
 * processes using a barrier or asking about MPI rank and size.
 */
class PythonMPIComm
{
public:

    /// Default constructor, providing an interface to MPI_COMM_WORLD.
    PythonMPIComm() :
            mComm(MPI_COMM_WORLD )
    {
    }

    /// Constructor taking a custom MPI Communicator.
    /**
     * @param Comm MPI Communicator
     */
    PythonMPIComm(MPI_Comm Comm) :
            mComm(Comm)
    {
    }

    ~PythonMPIComm()
    {
    }

    /// Returns the process' rank.
    /** @return Identifier for the MPI process, larger or equal than 0 and smaller than the total number of processes.
     */
    int rank()
    {
        int rank;
        MPI_Comm_rank(mComm, &rank);
        return rank;
    }

    /// Returns the MPI size.
    /** @return Total number of MPI processes
     */
    int size()
    {
        int size;
        MPI_Comm_size(mComm, &size);
        return size;
    }

    /// Stops execution until all MPI processes reach the call to this function.
    /** Used to provide synchronization between the different MPI processes
     */
    void barrier()
    {
        MPI_Barrier(mComm);
    }

//    void abort(int ErrCode)
//    {
//        MPI_Abort(mComm,ErrCode);
//    }

private:

    friend class PythonMPI;

    /// Return a reference to the internal MPI_Comm object wrapped by this class.
    MPI_Comm& GetMPIComm()
    {
        return mComm;
    }

    /// The MPI communicator wrapped by this class.
    MPI_Comm mComm;

};

/// A Python wrapper to common MPI funcions.
/** This class reimplements a very limited subset of the functionality available
 * in the Python interface, which provides enough capabilities for
 * the tasks commonly performed in Kratos scripts.
 */
class PythonMPI
{
public:

    enum MPI_Operation {
        MAX,
        MIN,
        SUM
    };

	/// Default constructor.
	/** Initializes MPI if required and defines a wrapper for MPI_COMM_WORLD,
	 * which can be accessed by calling GetWorld().
	 */
	PythonMPI()
	{
		int MpiIsInitialized = 0;

        MPI_Initialized(&MpiIsInitialized);

        if (MpiIsInitialized == 0)
        {
            int argc = 0;
            char* a = new char[1];
            *a = '\0';
            char** empty_argv = &a;

#if MPI_VERSION < 2
            MPI_Init(&argc, &empty_argv);
#else
                        int provided;
                        MPI_Init_thread(&argc, &empty_argv, MPI_THREAD_MULTIPLE, &provided);

                        if(provided < MPI_THREAD_MULTIPLE)
                        {
                            int rank;
                            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                            if(rank==0)
                                std::cout<< "MPI_Init_thread returns : " << provided << std::endl;

                        }
#endif
            delete[] a;
        }

        mWorld = PythonMPIComm(MPI_COMM_WORLD );
    }

    /// Destructor, finalizes MPI if necessary.
    ~PythonMPI()
    {
        int MpiIsFinalized = 0;

        MPI_Finalized(&MpiIsFinalized);

        if (MpiIsFinalized == 0)
            MPI_Finalize();
    }

    /**
     * @return A PythonMPIComm object wrapping MPI_COMM_WORLD
     */
    PythonMPIComm& GetWorld()
    {
        return mWorld;
    }

    /// Return MPI rank as given by the provided communicator.
    int rank(PythonMPIComm& rComm)
    {
        int rank;
        MPI_Comm_rank(rComm.GetMPIComm(), &rank);
        return rank;
    }

    /// Return MPI rank as given by MPI_COMM_WORLD.
    int rank()
    {
        return this->rank(mWorld);
    }

    /// Return MPI size as given by the provided communicator.
    int size(PythonMPIComm& rComm)
    {
        int size;
        MPI_Comm_size(rComm.GetMPIComm(), &size);
        return size;
    }

    /// Return MPI size as given by MPI_COMM_WORLD.
    int size()
    {
        return this->size(mWorld);
    }

    /// Synchronize processes using an MPI_Barrier call using given communicator.
    void barrier(PythonMPIComm& rComm)
    {
        MPI_Barrier(rComm.GetMPIComm());
    }

    void barrier()
    {
        this->barrier(mWorld);
    }

	/// Perform a MPI_Bcast operation.
	/**
	 * Broadcasting a value to all ranks
	 * @param rComm A communicator object.
	 * @param LocalValue The local value to be sent in the gather.
	 * @param RankToBroacastFrom The MPI rank of the process where the values will be broadcasted from.
	 * @return The broadcasted value on all ranks
	 */
	template<class TValueType>
	TValueType broadcast(PythonMPIComm& rComm,
                         TValueType LocalValue,
                         const int RankToBroacastFrom)
	{
        // Determine data type
		const MPI_Datatype DataType = this->GetMPIDatatype(LocalValue);

        MPI_Bcast(&LocalValue, 1, DataType, RankToBroacastFrom, rComm.GetMPIComm());

        return LocalValue;
    }

	/// Perform a MPI_Reduce operation.
	/**
	 * Perform a reduction given an MPI_Op-Type
	 * @param rComm A communicator object.
	 * @param LocalValue The local value to be sent in the gather.
	 * @param RankToReduceOn The MPI rank of the process where the values will be gathered.
	 * @param MPI_Operation The MPI_Op to be used for the reduction
	 * @return The reduced value on the for the RankToReduceOn thread
	 */
	template<class TValueType>
    TValueType reduce(PythonMPIComm& rComm,
                      const TValueType LocalValue,
                      const int RankToReduceOn,
                      const MPI_Operation MpiOp)
	{
		// Determine data type
		const MPI_Datatype DataType = this->GetMPIDatatype(LocalValue);

        TValueType result_val;
        MPI_Reduce(&LocalValue, &result_val, 1, DataType,
                   GetMPIOpType(MpiOp), RankToReduceOn, rComm.GetMPIComm());

        return result_val;
    }

	/// Perform a MPI_Allreduce operation.
	/**
	 * Perform a allreduction given an MPI_Op-Type
	 * @param rComm A communicator object.
	 * @param LocalValue The local value to be sent in the gather.
	 * @param MPI_Operation The MPI_Op to be used for the reduction
	 * @return The reduced value on all ranks
	 */
	template<class TValueType>
    TValueType allreduce(PythonMPIComm& rComm,
                         const TValueType LocalValue,
                         const MPI_Operation MpiOp)
	{
		// Determine data type
		const MPI_Datatype DataType = this->GetMPIDatatype(LocalValue);

        TValueType result_val;
        MPI_Allreduce(&LocalValue, &result_val, 1, DataType,
                      GetMPIOpType(MpiOp), rComm.GetMPIComm());

        return result_val;
    }

	/// Perform a MPI_Scatter operation.
	/**
	 * Scattering a vector of values from the RankToScatterFrom process.
	 * @param rComm A communicator object.
	 * @param rLocalValues The local values to be scattered.
	 * @param RankToScatterFrom The MPI rank of the process where the values will be scattered from.
	 * @return The local value that was scattered
	 */
    template<class TValueType>
	TValueType scatter(PythonMPIComm& rComm,
                       const std::vector<TValueType>& rLocalValues,
                       const int RankToScatterFrom)
    {
        // Determine data type
        const MPI_Datatype DataType = this->GetMPIDatatype(TValueType());
        int rank, size;
        MPI_Comm_rank(rComm.GetMPIComm(), &rank);
        MPI_Comm_size(rComm.GetMPIComm(), &size);
        if (rank == RankToScatterFrom && rLocalValues.size() != size)
            throw std::runtime_error("Wrong number of values to Scatter!");
        TValueType receive_val;
        // Communicate
        MPI_Scatter(rLocalValues.data(), 1, DataType, &receive_val,
                    1, DataType, RankToScatterFrom, rComm.GetMPIComm());
        return receive_val;
    }

	/// Perform a MPI_Scatterv operation.
	/**
	 * Scattering a number of vector with values from the RankToScatterFrom process.
	 * @param rComm A communicator object.
	 * @param rLocalValues The local values to be scattered.
	 * @param RankToScatterFrom The MPI rank of the process where the values will be scattered from.
	 * @return The vector of values that was scattered
	 */
    template<class TValueType>
	std::vector<TValueType> scatterv(PythonMPIComm& rComm,
                                     const std::vector<std::vector<TValueType>>& rLocalValues,
                                     const int RankToScatterFrom)
    {
        // Determine data type
        const MPI_Datatype DataType = this->GetMPIDatatype(TValueType());
        int rank, size;
        MPI_Comm_rank(rComm.GetMPIComm(), &rank);
        MPI_Comm_size(rComm.GetMPIComm(), &size);
        if (rank == RankToScatterFrom && rLocalValues.size() != size)
            throw std::runtime_error("Wrong number of values to Scatter!");

        std::vector<TValueType> send_buffer;
        std::vector<int> send_sizes(size);
        std::vector<int> displs(size);

        if (rank == RankToScatterFrom)
        {
            for (int i=0; i<size; ++i) {
                send_sizes[i] = rLocalValues[i].size();
            }
            // no padding is applied in the buffer
            const int sum_values = std::accumulate(send_sizes.begin(), send_sizes.end(), 0);
            send_buffer.resize(sum_values);

            int counter = 0;
            for (int i=0; i<size; ++i) {
                displs[i] = counter;
                for (int j=0; j<send_sizes[i]; ++j) {
                    send_buffer[counter] = rLocalValues[i][j];
                    ++counter;
                }
            }
        }

        // Communicate how much data each rank will receive
        int recv_count;
        MPI_Scatter(send_sizes.data(), 1, MPI_INT, &recv_count,
                    1, MPI_INT, RankToScatterFrom, rComm.GetMPIComm());

        std::vector<TValueType> scattered_vals(recv_count);

        // Communicate
        MPI_Scatterv(send_buffer.data(), send_sizes.data(), displs.data(),
                     DataType, scattered_vals.data(), recv_count,
                     DataType, RankToScatterFrom, rComm.GetMPIComm());

        return scattered_vals;
    }

	/// Perform a MPI_Gather operation.
	/**
	 * Provide a std::vector containing all local values to the RankToGatherOn process.
	 * @param rComm A communicator object.
	 * @param LocalValue The local value to be sent in the gather.
	 * @param RankToGatherOn The MPI rank of the process where the values will be gathered.
	 * @return A std::vector containing the local values in all processes, sorted by rank,
     * for the RankToGatherOn thread, an empty std::vector for other processes
	 */
	template<class TValueType>
	std::vector<TValueType> gather(PythonMPIComm& rComm,
                                   const TValueType LocalValue,
                                   const int RankToGatherOn)
    {
        // Determine data type
        const MPI_Datatype DataType = this->GetMPIDatatype(LocalValue);

        int rank, size;
        MPI_Comm_rank(rComm.GetMPIComm(), &rank);
        MPI_Comm_size(rComm.GetMPIComm(), &size);

        // Create recieve buffer
        std::vector<TValueType> global_values;
        if (rank == RankToGatherOn)
            global_values.resize(size);

        // Communicate
        MPI_Gather(&LocalValue, 1, DataType, global_values.data(),
                   1, DataType, RankToGatherOn, rComm.GetMPIComm());

        return global_values;
    }

    /// Perform an MPI_Gatherv operation.
    /**
     * Provide a std::vector containing all local values to the RankToGatherOn process.
     * @param rComm A communicator object.
     * @param rLocalValues The std::vector of local values to be sent in the gather.
     * @param RankToGatherOn The MPI rank of the process where the values will be gathered.
     * @return A std::vector containing the local value lists in all processes, sorted by rank,
     * for the RankToGatherOn thread, an empty std::vector for other processes
     */
    template<class TValueType>
    std::vector<std::vector<TValueType>> gatherv(PythonMPIComm& rComm,
                                   const std::vector<TValueType>& rLocalValues,
                                   const int RankToGatherOn)
    {
        // Determine data type
        const MPI_Datatype DataType = this->GetMPIDatatype(TValueType());

        int recv_block_size = 1;
        int size, rank, send_size;

        MPI_Comm_size(rComm.GetMPIComm(), &size);
        MPI_Comm_rank(rComm.GetMPIComm(), &rank);
        send_size = rLocalValues.size();

        std::vector<int> recv_sizes(size);
        std::vector<int> displs(size);
        MPI_Gather(&send_size, 1, MPI_INT, recv_sizes.data(),
                   1, MPI_INT, RankToGatherOn, rComm.GetMPIComm());

        // calculate receive buffer block size for root
        if (rank == RankToGatherOn)
        {
            for (int i = 0; i < size; ++i)
                recv_block_size = (recv_block_size < recv_sizes[i]) ? recv_sizes[i] : recv_block_size;
            for (int i = 0; i < size; ++i)
                displs[i] = i * recv_block_size;
        }

        std::vector<TValueType> recv_buffer(recv_block_size * size);

        // gather local arrays at root
        MPI_Gatherv(rLocalValues.data(), send_size, DataType, recv_buffer.data(),
                    recv_sizes.data(), displs.data(), DataType, RankToGatherOn, rComm.GetMPIComm());

        std::vector<std::vector<TValueType>> condensed_vector;

        // now condensing the vector such that is has only the number of actual
        // elements that come from each rank (i.e. removing the buffering that is
        // needed for the mpi-call)
        if (rank == RankToGatherOn)
        {
            condensed_vector.resize(size);
            const auto buffer_begin = recv_buffer.begin();
            for (int i = 0; i < size; ++i)
            {
                const auto iter_start = buffer_begin + i * recv_block_size;
                const auto iter_end = buffer_begin + i * recv_block_size + recv_sizes[i];
                condensed_vector[i] = std::vector<TValueType>(iter_start, iter_end);
            }
        }

        return condensed_vector;
    }

    /// Perform an MPI_allgather operation.
    /**
     * Provide a std::vector containing all local values to all processes.
     * @param rComm A communicator object.
     * @param LocalValue The local value to be sent in the gather.
     * @return A std::vector containing the local values in all processes, sorted by rank.
     */
    template<class TValueType>
    std::vector<TValueType> allgather(PythonMPIComm& rComm,
                                      const TValueType LocalValue)
    {
        // Determine data type
        const MPI_Datatype DataType = this->GetMPIDatatype(LocalValue);

        int size;
        MPI_Comm_size(rComm.GetMPIComm(), &size);

        // Create recieve buffer
        std::vector<TValueType> global_values(size);

        // Communicate
        MPI_Allgather(&LocalValue, 1, DataType, global_values.data(),
                      1, DataType, rComm.GetMPIComm());

        return global_values;
    }

private:

    MPI_Comm& GetMPIComm(PythonMPIComm Comm)
    {
        return Comm.GetMPIComm();
    }

    /// An auxiliary function to determine the MPI_Datatype corresponding to a given C type
    template<class T>
    inline MPI_Datatype GetMPIDatatype(const T& Value);

    /// An auxiliary function to determine the MPI_Op corresponding to a given Enum type
    /// This is necessary bcs MPI_Op cannot be directly exposed to Python
    inline MPI_Op GetMPIOpType(const MPI_Operation MPIOpEnum)
    {
        switch(MPIOpEnum) {
	        case MAX: return MPI_MAX;
	        case MIN: return MPI_MIN;
	        case SUM: return MPI_SUM;
        	default: break;
        }
    }

	int mArgc;

    char** mArgv;

    PythonMPIComm mWorld;

};

template<>
inline MPI_Datatype PythonMPI::GetMPIDatatype<int>(const int& Value)
{
    return MPI_INT ;
}

template<>
inline MPI_Datatype PythonMPI::GetMPIDatatype<double>(const double& Value)
{
    return MPI_DOUBLE ;
}

PythonMPI& GetMPIInterface()
{
    static PythonMPI ThePythonMPI;
    return ThePythonMPI;
}

} // Namespace Kratos

#endif // KRATOS_MPI_PYTHON_H
