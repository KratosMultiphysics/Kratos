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


#if !defined(KRATOS_FIXED_SIZE_MEMORY_POOL_H_INCLUDED )
#define  KRATOS_FIXED_SIZE_MEMORY_POOL_H_INCLUDED

#include <vector>
#include <list>

#include "includes/thread_fixed_size_memory_pool.h"

namespace Kratos
{
  ///@addtogroup KratosCore
  ///@{

  ///@name Kratos Classes
  ///@{

  /// FixedSizeMemoryPool is the multi-thread manager of Kratos memory management.
  /** The memory management of Kratos is implemented based on the design
	  given in Modern C++ Design by A. Alexandrescu. However a new layer is added
	  over the ThreadFixedSizeMemoryPool in order to control the multi-thread
	  creation and destruction of the chunks. This layer keeps ThreadFixedSizeMemoryPool
	  for each thread and is in charge of calling their methods in a thread safe way
	  and also feeds them with new chunk (or reusing the released ones)
  */
  class FixedSizeMemoryPool : public LockObject
    {
    public:
	  ///@name Type Definitions
	  ///@{

		using SizeType = Chunk::SizeType;

	  ///@}

		static constexpr SizeType DefaultChunkSize = 1024*1024; // 1M byte

	  ///@name Life Cycle
      ///@{

      /// Default constructor is deleted.
      FixedSizeMemoryPool() = delete;

	  /// Copy constructor is deleted.
	  FixedSizeMemoryPool(FixedSizeMemoryPool const& rOther) = delete;

	  /// The constructor to be called
	  FixedSizeMemoryPool(std::size_t BlockSizeInBytes, SizeType ChunkSize = DefaultChunkSize)
		  : LockObject()
		  , mChunkSize(ChunkSize)
	  {
		  for (int i_thread = 0; i_thread < OpenMPUtils::GetCurrentNumberOfThreads(); i_thread++)
			  mThreadsPool.emplace_back(BlockSizeInBytes, ChunkSize, i_thread);
	  }

      /// Destructor
	  virtual ~FixedSizeMemoryPool() {

	  }

      ///@}
      ///@name Operators
      ///@{

	  /// Assignment operator is deleted.
	  FixedSizeMemoryPool& operator=(FixedSizeMemoryPool const& rOther) = delete;

      ///@}
      ///@name Operations
      ///@{

	  /// This function does not throw and returns zero if cannot allocate
	  void* Allocate() {
		  void* p_result = mThreadsPool[OpenMPUtils::ThisThread()].Allocate();
		  return p_result;
	  }

	  void Deallocate(void* pPointrerToRelease) {

		  if (mThreadsPool[OpenMPUtils::ThisThread()].Deallocate(pPointrerToRelease))
		  {
			  return;
		  }

		  for (int i_thread = 0; i_thread < OpenMPUtils::GetCurrentNumberOfThreads(); i_thread++)
			  if (i_thread != OpenMPUtils::ThisThread())
				  if (mThreadsPool[i_thread].Deallocate(pPointrerToRelease)) {
					  return;
				  }

		  KRATOS_ERROR << "The Pointer with address " << pPointrerToRelease << " was not found in this pool" << std::endl;
	  }

	  std::size_t MemoryUsed() const {
		  std::size_t memory_used = sizeof(FixedSizeMemoryPool);
		  for (int i_thread = 0; i_thread < OpenMPUtils::GetCurrentNumberOfThreads(); i_thread++)
			  memory_used += mThreadsPool[i_thread].MemoryUsed();
		  return memory_used;
	  }

	  std::size_t MemoryOverhead() const {
		  std::size_t memory_overhead = sizeof(FixedSizeMemoryPool);
		  for (int i_thread = 0; i_thread < OpenMPUtils::GetCurrentNumberOfThreads(); i_thread++)
			  memory_overhead += mThreadsPool[i_thread].MemoryOverhead();
		  return memory_overhead;
	  }

	  std::size_t GetNumberOfAllocatedChunks() const {
		  std::size_t number_of_allocated_chunks = 0;
		  for (int i_thread = 0; i_thread < OpenMPUtils::GetCurrentNumberOfThreads(); i_thread++)
			  number_of_allocated_chunks += mThreadsPool[i_thread].GetNumberOfChunks() - mThreadsPool[i_thread].GetNumberOfReleasedChunks();
		  return number_of_allocated_chunks;
	  }

	  ///@}
      ///@name Access
      ///@{

	  SizeType ChunkSize() const {
		  return mChunkSize;
	  }

	  ///@}
      ///@name Inquiry
      ///@{

      ///@}
      ///@name Input and output
      ///@{

      /// Turn back information as a string.
	  std::string Info() const {
		  return "FixedSizeMemoryPool";
	  }

      /// Print information about this object.
      void PrintInfo(std::ostream& rOStream) const {
		  rOStream << Info();
	  }

      /// Print object's data.
	  void PrintData(std::ostream& rOStream) const {
		  std::size_t memory_used = MemoryUsed();
		  std::size_t memory_overhead = MemoryOverhead();
		  double overhead_percentage = memory_overhead;
		  if (memory_overhead < memory_used)
			  overhead_percentage = static_cast<double>(memory_overhead)/(memory_used - memory_overhead);
		  overhead_percentage *= 100.00;

		  rOStream << GetNumberOfAllocatedChunks() << " Chunks of "
			  << SizeInBytesToString(ChunkSize()) << " bytes each. Total memory usage: "
			  << SizeInBytesToString(MemoryUsed()) << " bytes and memory overhead "
			  << SizeInBytesToString(MemoryOverhead()) << "(" << overhead_percentage << "%)" << std::endl;
	  }

      ///@}


    private:
      ///@name Member Variables
      ///@{

		SizeType mChunkSize;
		std::vector<ThreadFixedSizeMemoryPool> mThreadsPool;

      ///@}
      ///@name Operations
      ///@{

		std::string SizeInBytesToString(std::size_t Bytes) const {
			std::stringstream buffer;
			double result = Bytes;
			constexpr int units_size = 5;
			constexpr char units[units_size] = { ' ', 'k','M','G','T' };
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

    }; // Class FixedSizeMemoryPool

  ///@}
  ///@name Input and output
  ///@{


  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const FixedSizeMemoryPool& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_FIXED_SIZE_MEMORY_POOL_H_INCLUDED  defined
