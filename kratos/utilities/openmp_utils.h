//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#ifndef KRATOS_OPENMP_UTILS_H
#define	KRATOS_OPENMP_UTILS_H

#include <stdio.h>
#include <vector>
#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#else
#include <ctime>
#endif

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/// Implements basic tasks for OpenMP parallelism and suitable scalar alternatives
/**
 This class defines utility functions that implement some basic OpenMP
 capabilities and an equivalent scalar alternative to use in compilations
 where OpenMP is not enabled. The idea is to allow Kratos developers to
 design their code in parallel, knowing that it will work in scalar runs
 as well.
 */
class OpenMPUtils
{
public:

    ///@name Type definitions
    ///@{

    /// Vector type for the output of DivideInPartitions method
    /**
     *  @see OpenMPUtils::DivideInPartitions
     */
    typedef std::vector<int> PartitionVector;

    ///@}
    ///@name Operations
    ///@{

    /// Wrapper for omp_get_max_threads().
    /**
     @return Maximum number of OpenMP threads that will be used in
     parallel regions.
     */
    static inline int GetNumThreads()
    {
#ifdef _OPENMP
        return omp_get_max_threads();
#else
        return 1;
#endif
    }

    /// Wrapper for omp_get_num_threads().
    /**
     @return Number of OpenMP threads in the current team.
     */
    static int GetCurrentNumberOfThreads()
    {
#ifdef _OPENMP
        return omp_get_num_threads();
#else
        return 1;
#endif
    }

    /// Wrapper for omp_get_num_procs().
    /**
     @return Number of processors available to the device.
     */
    static int GetNumberOfProcessors()
    {
#ifdef _OPENMP
        return omp_get_num_procs();
#else
        return 1;
#endif
    }

    /// Wrapper for omp_get_dynamic().
    /**
     @return Dynamic teams are enabled.
     */
    static int IsDynamic()
    {
#ifdef _OPENMP
        return omp_get_dynamic();
#else
        return 0;
#endif
    }

    /// Wrapper for omp_get_thread_num().
    /**
     @return The thread number for this thread, 0 if scalar run.
     */
    static inline int ThisThread()
    {
#ifdef _OPENMP
        return omp_get_thread_num();
#else
        return 0;
#endif
    }

    /// Wrapper for omp_in_parallel().
    /**
     @return Maximum number of OpenMP threads that will be used in
     parallel regions.
     */
    static inline int IsInParallel()
    {
#ifdef _OPENMP
        return omp_in_parallel();
#else
        return 0;
#endif
    }

    /// Timing routine.
    /**
     Determine the current time by calling an appropiate
     (scalar or parallel) timer class.
     @return Current time
     */
    static double GetCurrentTime()
    {
#ifndef _OPENMP
        return std::clock()/static_cast<double>(CLOCKS_PER_SEC);
#else
        return  omp_get_wtime();
#endif
    }

    /// Divide an array of length NumTerms between NumThreads threads.
    /**
     Creates a std::vector containing NumThreads + 1 terms, where term k
     is the first and position of the array that corresponds to thread k.
     The k+1 term is the end of the array, so that the vector can be used
     to iterate the array between 'k' and 'k+1' in each thread.
     @param NumTerms Number of objects to be divided between the threads.
     @param NumThreads The number of parallel threads that will be used.
     @param Partitions This object will contain the begin and end positions
     for each thread.
     */
    static inline void DivideInPartitions(
        const int NumTerms,
        const int NumThreads,
        PartitionVector& Partitions)
    {
        Partitions.resize(NumThreads + 1);
        int PartitionSize = NumTerms / NumThreads;
        Partitions[0] = 0;
        Partitions[NumThreads] = NumTerms;
        for(int i = 1; i < NumThreads; i++)
            Partitions[i] = Partitions[i-1] + PartitionSize ;
    }

    /// Generate a partition for an std::vector-like array, providing iterators to the begin and end positions for each thread.
    /** This function assumes that the vector class will have an iterator type and implement begin(), end() and size() methods.
      * @param rVector An arary containing the elements to be distributed between the threads.
      * @param rBegin Iterator pointing to the first element in rVector to be used in the current thread.
      * @param rEnd Iterator pointing to the end position for the current thread in rVector.
      */
    template< class TVector >
    static void PartitionedIterators(TVector& rVector,
                                     typename TVector::iterator& rBegin,
                                     typename TVector::iterator& rEnd)
    {
#ifdef _OPENMP
        int NumTerms = rVector.size();
        int ThreadNum = omp_get_thread_num();
        int NumThreads = omp_get_max_threads();
        int PartitionSize = NumTerms / NumThreads;
        // Set Partition start
        rBegin = rVector.begin() + ThreadNum * PartitionSize;
        // Partition ends after 'PartitionSize' terms, except if this is the last partition
        if ( (ThreadNum + 1) != NumThreads )
            rEnd = rBegin + PartitionSize;
        else
            rEnd = rVector.end();
#else
        rBegin = rVector.begin();
        rEnd = rVector.end();
#endif
    }

    /// A function to set the number of threads from Python.
    /**
     This is an auxiliary mainly intended for test purposes, to help with the
     detection of race conditions.
     @param NumThreads Number of threads to use in parallel regions. Note
     that values greater than the environment variable OMP_NUM_THREADS
     will be ignored.
     */
    static inline void SetNumThreads(int NumThreads = 1)
    {
#ifdef _OPENMP

      int procs    = omp_get_num_procs();
      if( procs < NumThreads ){
	std::cout<<" WARNING: Maximimun number of threads is EXCEEDED "<<std::endl;
	/* Set thread number */
	omp_set_num_threads(procs);
	std::cout<<" Number of Threads Set To : "<<procs<<std::endl;
      }
      else{
	/* Set thread number */
	omp_set_num_threads(NumThreads);
      }

#endif
    }

    /**
     A method to print the OMP information
     */
    static inline void PrintOMPInfo()
    {
#ifdef _OPENMP

      int nthreads,tid, procs, maxt, inpar, dynamic, nested;

      /* Start parallel region */

#pragma omp parallel private(nthreads, tid)
      {
	/* Obtain thread number */
	tid = omp_get_thread_num();

	/* Only master thread does this */
	if (tid == 0)
	  {
	    printf("  Thread %d getting environment info...\n", tid);

	    /* Get environment information */
	    procs    = omp_get_num_procs();
	    nthreads = omp_get_num_threads();
	    maxt     = omp_get_max_threads();
	    inpar    = omp_in_parallel();
	    //omp_set_dynamic(true);
	    dynamic  = omp_get_dynamic();
	    //omp_set_nested(true);
	    nested   = omp_get_nested();

	    /* Print environment information */
	    printf( "  | ------------ OMP IN USE --------- |\n");
	    printf( "  | Machine number of processors  = %d |\n", procs);
	    printf( "  | Number of threads set         = %d |\n", nthreads);
	    printf( "  | Max threads in use            = %d |\n", maxt);
	    printf( "  | In parallel?                  = %d |\n", inpar);
	    printf( "  | Dynamic threads enabled?      = %d |\n", dynamic);
	    printf( "  | Nested parallelism supported? = %d |\n", nested);
	    printf( "  | --------------------------------- |\n");


	    if( procs < nthreads )
	      std::cout<<" ( WARNING: Maximimun number of threads is EXCEEDED )"<<std::endl;

	  }

      }

#endif
    }

    template<class T>
    static inline void CreatePartition(unsigned int number_of_threads, const int number_of_rows, T& partitions)
    {
        partitions.resize(number_of_threads+1);
        int partition_size = number_of_rows / number_of_threads;
        partitions[0] = 0;
        partitions[number_of_threads] = number_of_rows;
        for(unsigned int i = 1; i<number_of_threads; i++)
            partitions[i] = partitions[i-1] + partition_size ;
    }

    ///@} //Operations
};

///@} //Kratos classes

///@} addtogroup block
}

#endif	/* KRATOS_OPENMP_UTILS_H */

