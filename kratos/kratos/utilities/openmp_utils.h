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

#ifdef _OPENMP
#include <omp.h>
#endif

namespace Kratos
{
    /* This class defines utility functions that implement some basic OpenMP
     * capabilities and an equivalent scalar alternative to use in compilations
     * where OpenMP is not enabled. The idea is to allow Kratos developers to
     * design their code in parallel, knowing that it will work in scalar runs
     * as well.
     */
    class OpenMPUtils
    {
    public:

        typedef std::vector<int> PartitionVector;

        /* Wrapper for omp_get_max_threads() */
        static inline int GetNumThreads()
        {
            #ifdef _OPENMP
            return omp_get_max_threads();
            #else
            return 1;
            #endif
        }

        /* Wrapper for omp_get_thread_num() */
        static inline int ThisThread()
        {
            #ifdef _OPENMP
            return omp_get_thread_num();
            #else
            return 0;
            #endif
        }

        /* Divide an array of length NumTerms between NumThreads threads.
         * Creates a std::vector containing NumThreads + 1 terms, where term k
         * is the first and position of the array that corresponds to thread k.
         * The k+1 term is the end of the array, so that the vector can be used
         * to iterate the array between 'k' and 'k+1' in each thread.
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
    };
}

#endif	/* KRATOS_OPENMP_UTILS_H */

