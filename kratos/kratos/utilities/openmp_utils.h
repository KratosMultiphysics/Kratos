/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

/*
 * File:   openmp_utils.h
 * Author: jcotela
 *
 * Created on 15 June 2010, 10:40
 */

#ifndef KRATOS_OPENMP_UTILS_H
#define	KRATOS_OPENMP_UTILS_H

#include <stdio.h>
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
            #ifdef _OPENMP
            Partitions.resize(NumThreads + 1);
            int PartitionSize = NumTerms / NumThreads;
            Partitions[0] = 0;
            Partitions[NumThreads] = NumTerms;
            for(int i = 1; i < NumThreads; i++)
                Partitions[i] = Partitions[i-1] + PartitionSize ;
            #else
            Partitions.resize(2);
            Partitions[0] = 0;
            Partitions[1] = NumTerms;
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
        static inline void SetNumThreads(int NumThreads)
        {
            #ifdef _OPENMP
            omp_set_num_threads(NumThreads);
            std::cout << "Maximum number of threads is now " << omp_get_max_threads() << std::endl;
            #endif
        }

       static inline void CreatePartition(unsigned int number_of_threads, const int number_of_rows, vector<unsigned int>& partitions)
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

