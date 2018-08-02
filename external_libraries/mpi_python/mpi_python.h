#ifndef KRATOS_MPI_PYTHON_H
#define KRATOS_MPI_PYTHON_H

#include <stdio.h>
#include <iostream>

#include "mpi.h"
#include <pybind11/pybind11.h>


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
 * in the Python interface for boost::mpi, which provides enough capabilities for
 * the tasks commonly performed in Kratos scripts.
 */
class PythonMPI
{
public:

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

	/// Perform a MPI_Gather operation.
	/**
	 * Provide a list containing all local values to the Root process.
	 * @param rComm A communicator object.
	 * @param LocalValue The local value to be sent in the gather.
	 * @param Root The MPI rank of the process where the valued will be gathered.
	 * @return A Python list containing the local values in all processes, sorted by rank, for the Root thread, an empty Python list for other processes
	 */
	template<class TValueType>
	pybind11::list gather(PythonMPIComm& rComm, TValueType LocalValue,
	        int Root)
	{
		// Determime data type
		MPI_Datatype DataType = this->GetMPIDatatype(LocalValue);

		int size;
		MPI_Comm_size(rComm.GetMPIComm(), &size);

		int rank;
		MPI_Comm_rank(rComm.GetMPIComm(), &rank);

		// Create recieve buffer
		TValueType* GlobalValues;
		if (rank == Root)
			GlobalValues = new TValueType[size];

		// Communicate
		MPI_Gather(&LocalValue, 1, DataType, GlobalValues, 1, DataType, Root,
		        rComm.GetMPIComm());

		// Copy output to a Python list
		pybind11::list Out;

		if (rank == Root)
		{
			for (int i = 0; i < size; i++)
			{
				auto val(GlobalValues[i]);
				Out.append(val);
			}
			delete[] GlobalValues;
		}

		return Out;
	}

	/// Perform an MPI_Gather operation.
	/**
	 * Provide a list containing all local values of type double to the Root process.
	 * @param rComm A communicator object.
	 * @param LocalValues The list of local values to be sent in the gather.
	 * @param Root The MPI rank of the process where the valued will be gathered.
	 * @return A Python list containing the local value lists in all processes, sorted by rank, for the Root thread, an empty Python list for other processes
	 */
	pybind11::list gather(PythonMPIComm& rComm,
	        pybind11::list LocalValues, int Root)
	{
		int RecvBlockSize = 1;
	        int Size, Rank, SendSize;
		int * RecvSize, * Displs;
		double * SendBuf, * RecvBuf;

		MPI_Comm_size(rComm.GetMPIComm(), &Size);
		MPI_Comm_rank(rComm.GetMPIComm(), &Rank);
		SendSize = len(LocalValues);

		RecvSize = new int[Size];
		Displs = new int[Size];
		MPI_Gather(&SendSize, 1, MPI_INT, RecvSize, 1, MPI_INT, Root,
		        rComm.GetMPIComm());

		// calculate receive buffer block size for root
		if (Rank == Root)
		{
			for (int i = 0; i < Size; i++)
				RecvBlockSize = (RecvBlockSize < RecvSize[i]) ? RecvSize[i] : RecvBlockSize;
			for (int i = 0; i < Size; i++)
			        Displs[i] = i * RecvBlockSize;
		}

		SendBuf = new double[SendSize];
		for (int i = 0; i < SendSize; i++)
		  SendBuf[i] = pybind11::cast<double>(LocalValues[i]);

		RecvBuf = new double[RecvBlockSize * Size];

		// gather local arrays at root
		MPI_Gatherv(SendBuf, SendSize, MPI_DOUBLE, RecvBuf, RecvSize, Displs, MPI_DOUBLE, Root, rComm.GetMPIComm());

		pybind11::list Out;

		// make python list of gathered arrays
		if (Rank == Root)
		{
			for (int i = 0; i < Size; i++)
			{
				int iblock = i * RecvBlockSize;
				pybind11::list Outi;
				for (int j = 0; j < RecvSize[i]; j++)
					Outi.append(RecvBuf[iblock + j]);
				Out.append(Outi);
			}
		}

		delete[] RecvSize;
		delete[] Displs;
		delete[] SendBuf;
		delete[] RecvBuf;

		return Out;
	}

	/// Perform an MPI_allgather operation.
	/**
	 * Provide a list containing all local values to all processes.
	 * @param rComm A communicator object.
	 * @param LocalValue The local value to be sent in the gather.
	 * @return A Python list containing the local values in all processes, sorted by rank.
	 */
	template<class TValueType>
	std::vector<TValueType> allgather(PythonMPIComm& rComm, TValueType LocalValue)
	{
		// Determime data type
		const MPI_Datatype DataType = this->GetMPIDatatype(LocalValue);

		int size;
		MPI_Comm_size(rComm.GetMPIComm(), &size);

		// Create recieve buffer
        std::vector<TValueType> global_values(size);

		// Communicate
		MPI_Allgather(&LocalValue, 1, DataType, global_values.data(), 1, DataType,
		        rComm.GetMPIComm());

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

}

#endif // KRATOS_MPI_PYTHON_H
