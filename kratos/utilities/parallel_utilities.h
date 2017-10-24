#pragma once

//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//
                



// System includes

// External includes 

// Project includes
#include "includes/define.h"
#include "utilities/openmp_utils.h"


namespace Kratos
{
    namespace Parallel
    {
    ///@addtogroup Kratos
    ///@{

    ///@name Kratos Globals
    ///@{
        //this function implements a parallel for, in which f is a lambda function which captures its parameters once per very time it is called
        template<class InputIt, class UnaryFunction>
        inline void parallel_for(InputIt first, const int n, UnaryFunction f)
        {
            
            #pragma omp parallel for firstprivate(first)
            for (int i = 0; i < n; ++i)
            {
                f(first+i);
            }
        }
        
        //parallel for executed by chunks. Here the function is called on every "block of unknowns", so that the lambda captures are done just once
        //note that the function b must implement the loop from start to end
        //work allocation is static
        template<class TContainterType, class BinaryFunction>
        inline void parallel_for(TContainterType& data, UnaryFunction f)
        {
            ::Parallel::parallel_for(data.begin(), data.size(), f)
        }
        
        
        
        
        
        
        
        
        
        //the function f is executed in parallel. The parameter passed to f is the number of the thread executing it. 
        //one lambda capture is passed per every thread
        template< class UnaryFunction>
        inline execute(UnaryFunction f)
        {
            const int NumThreads = OpenMPUtils::GetNumThreads();

            PartitionVector& Partitions;
            #pragma omp parallel for firstprivate(first) static
            for (int i = 0; i < NumThreads; ++i) 
            {
                f(i);
            }
        }
        

        
        
        
        
        
        
        
        
        //parallel for executed by chunks. Here the function is called on every "block of unknowns", so that the lambda captures are done just once
        //note that the function b must implement the loop from start to end
        //work allocation is static
        template<class InputIt, class BinaryFunction>
        inline void block_parallel_for(InputIt first, const int n, BinaryFunction b)
        {
            int NumThreads = OpenMPUtils::GetNumThreads();
            OpenMPUtils::PartitionVector Partitions();
            OpenMPUtils::DivideInPartitions(n,NumThreads,Partitions);

            PartitionVector& Partitions;
            #pragma omp parallel for firstprivate(first) static
            for (int i = 0; i < NumThreads; ++i) 
            {
                b(first+Partitions[i], first+Partitions[i+1]);
            }
        }

        //parallel for executed by chunks. Here the function is called on every "block of unknowns", so that the lambda captures are done just once
        //note that the function b must implement the loop from start to end
        //work allocation is dynamic
        template<class InputIt, class BinaryFunction>
        inline void block_parallel_for(InputIt first, const int n, cosnt int NChunks, BinaryFunction f)
        {
            OpenMPUtils::PartitionVector Partitions();
            OpenMPUtils::DivideInPartitions(n,NChunks,Partitions);

            PartitionVector& Partitions;
            #pragma omp parallel for firstprivate(first, NChunks) dynamic
            for (int i = 0; i < NChunks; ++i) 
            {
                f(first+Partitions[i], first+Partitions[i+1]);
            }
        }
        
        template<class TContainterType, class BinaryFunction>
        inline void block_parallel_for(TContainterType& data, BinaryFunction b)
        {
            ::Parallel::block_parallel_for(data.begin(), data.size(), b);
        }

        template<class TContainterType, class BinaryFunction>
        inline void block_parallel_for(TContainterType& data, const int NChunks,  BinaryFunction b)
        {
            ::Parallel::block_parallel_for(data.begin(), data.size(), NChunks, b);
        }
        
        
        
        
        
        
        
        
        
        
        
        inline void AtomicAdd(const double& value, double& target)
        {
            #pragma omp atomic
            target += value;
        }
        inline void AtomicAdd(const int& value, int& target)
        {
            #pragma omp atomic
            target += value;
        }
        
        inline void AtomicSub(const double& value, double& target)
        {
            #pragma omp atomic
            target -= value;
        }
        inline void AtomicSub(const int& value, int& target)
        {
            #pragma omp atomic
            target -= value;
        }

        inline void AtomicAssign(const double& value, double& target)
        {
            #pragma omp atomic write
            target = value;
        }
        inline void AtomicAssign(const int& value, int& target)
        {
            #pragma omp atomic write
            target = value;
        }
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
    
    ///@} 
    
    ///@name Type Definitions       
    ///@{ 
    
    
    ///@} 
    ///@name Input and output 
    ///@{ 
            
    ///@}

    ///@} addtogroup block
  
    } //namespace Parallel
    
}  // namespace Kratos.

#endif // KRATOS_FILENAME_H_INCLUDED  defined 


