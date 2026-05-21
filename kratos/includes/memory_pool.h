//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   
//


#if !defined(KRATOS_MEMORY_POOL_H_INCLUDED )
#define  KRATOS_MEMORY_POOL_H_INCLUDED

#include <memory>
#include <vector>
#include <sstream>

#include "includes/fixed_size_memory_pool.h"

namespace Kratos
{
  ///@addtogroup KratosCore
  ///@{

  ///@name Kratos Classes
  ///@{

  /// MemoryPool is the smallest building block of Kratos memory management.
  /** The memory management of Kratos is implemented based on the design
	  given in Modern C++ Design by A. Alexandrescu and MemoryPool is the third
	  layer of it "AKA SmallObjAllocator" holding a pool of fixed size memory pools.
  */
  class MemoryPool
    {
    public:    
	  ///@name Type Definitions
	  ///@{

		using PoolsContainerType = std::vector<FixedSizeMemoryPool*>;

	  ///@}
	  ///@name Life Cycle
      ///@{

	  /// Copy constructor is deleted.
	  MemoryPool(MemoryPool const& rOther) = delete;

      /// Destructor 
	  virtual ~MemoryPool() {
		   for (auto i_pool = GetInstance().mPools.begin(); i_pool != GetInstance().mPools.end(); i_pool++)
			   delete *i_pool;
	  }

      ///@}
      ///@name Operators
      ///@{

	  /// Assignment operator is deleted.
	  MemoryPool& operator=(MemoryPool const& rOther) = delete;

      ///@}
      ///@name Operations
      ///@{

	  ///// Adding a memory pool by name and get the pointer. The name is only for report
	  //template <typename TPoolType = FixedSizeMemoryPool>
	  //static FixedSizeMemoryPool* AddPool(std::string const& PoolName, TPoolType* pHeapAllocatedPool) {
		 // KRATOS_ERROR_IF(HasPool(PoolName)) << "A duplicated memory pool found! The pool \"" << PoolName << "\" is already added." << std::endl;
		 // GetInstance().mPools[PoolName] = pHeapAllocatedPool;
		 // return pHeapAllocatedPool;
	  //}

	  static void* Allocate(std::size_t ObjectSizeInBytes) {
		  return GetPoolWithBlockSize(ObjectSizeInBytes)->Allocate();
	  }

	  static void Deallocate(void* pPointrerToRelease, std::size_t ObjectSizeInBytes) {
		  GetPoolWithBlockSize(ObjectSizeInBytes)->Deallocate(pPointrerToRelease);
	  }


	  ///@}
      ///@name Access
      ///@{

	  static MemoryPool& GetInstance() {
		  static MemoryPool instance;
		  return instance;
	  }

	  static std::size_t GetNumberOfPools()  {
		  return GetInstance().mPools.size();
	  }

	  static FixedSizeMemoryPool* GetPoolWithBlockSize(std::size_t BlockSize) {
		  PoolsContainerType& r_pools = GetInstance().mPools;

		  // I would avoid it by defining a max block size, but before that I should profile to see if it really slows down. Pooyan.
		  if (r_pools.size() <= BlockSize) { // This check is extra but is to avoid critical each time
#pragma omp critical
		  {
			  if (r_pools.size() <= BlockSize) // checking again to be sure some other thread doesn't change it meanwhile
				  r_pools.resize(BlockSize + 1, nullptr);
		  }
		  }

		  if (r_pools[BlockSize] == nullptr) { // This check is extra but is to avoid critical each time
#pragma omp critical
		  {
			  if (r_pools[BlockSize] == nullptr) // checking again to be sure some other thread doesn't change it meanwhile
				  r_pools[BlockSize] = new FixedSizeMemoryPool(BlockSize);
		  }
		  }
		  return r_pools[BlockSize];
	  }

      ///@}
      ///@name Inquiry
      ///@{

	  //static bool HasPool(std::string const& PoolName) {
		 // return (GetInstance().mPools.find(PoolName) != GetInstance().mPools.end());
	  //}

	  static std::size_t MemoryUsed() {
		  std::size_t result = sizeof(MemoryPool);
		  for (auto i_pool = GetInstance().mPools.begin(); i_pool != GetInstance().mPools.end(); i_pool++)
			  if(*i_pool != nullptr)
				result += (*i_pool)->MemoryUsed();
		  return result;
	  }

	  static std::size_t MemoryOverhead() {
		  std::size_t result = sizeof(MemoryPool);
		  for (auto i_pool = GetInstance().mPools.begin(); i_pool != GetInstance().mPools.end(); i_pool++)
			  if (*i_pool != nullptr)
				  result += (*i_pool)->MemoryOverhead();
		  return result;
	  }

      ///@}
      ///@name Input and output
      ///@{

      /// Turn back information as a string.
	  static std::string Info()  {
		  std::stringstream buffer("MemoryPool");
		  std::size_t memory_used = MemoryUsed();
		  std::size_t memory_overhead = MemoryOverhead();
		  double overhead_percentage = memory_overhead;
		  if (memory_overhead < memory_used)
			  overhead_percentage = static_cast<double>(memory_overhead)/(memory_used - memory_overhead);
		  overhead_percentage *= 100.00;

		  buffer << "Total memory usage: " 
			  << SizeInBytesToString(MemoryUsed()) << " bytes and memory overhead " 
			  << SizeInBytesToString(MemoryOverhead()) << "(" << overhead_percentage << "%)" << std::endl;

		  return buffer.str();
	  }

 
      ///@}


    private:

		///@name Life Cycle
		///@{

		/// MemoryPool cannot be created from outside. To ensure that the one created by instance is the only one.
		MemoryPool(){}

		///@}
      ///@name Member Variables
      ///@{

		PoolsContainerType mPools;

      ///@}
      ///@name Operations
      ///@{

		static std::string SizeInBytesToString(std::size_t Bytes) {
			std::stringstream buffer;
			double result = Bytes;
			constexpr int units_size = 5;
			constexpr char units[units_size+1] = { ' ', 'k','M','G','T', 'T'};
			int i = 0;
			for (; i < units_size; i++)
				if (result > 1024)
				{
					result /= 1024;
				}
				else
					break;
			buffer << result << units[i];

			return buffer.str();
		}

     ///@}

    }; // Class MemoryPool

  ///@}
  ///@name Input and output
  ///@{

  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MEMORY_POOL_H_INCLUDED  defined
