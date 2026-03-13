//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//
//

#pragma once

// System includes
#include <vector>
#include <list>
#include <algorithm>
#include <atomic>

// External includes
#include "concurrentqueue/concurrentqueue.h"

// Project includes
#include "includes/chunk.h"

namespace Kratos
{
	std::atomic_flag ThreadFixedSizeMemoryPoolLock = ATOMIC_FLAG_INIT ;

  ///@addtogroup KratosCore
  ///@{

  ///@name Kratos Classes
  ///@{

  /// ThreadFixedSizeMemoryPool holds chunks belong to a certain thread and operate over them.
  /** The memory management of Kratos is implemented based on the design
	  given in Modern C++ Design by A. Alexandrescu and ThreadFixedSizeMemoryPool is the second
	  layer of it "Called FixedAllocator" holding chunks.
	  This class is the owner of the chunks of this thread and the only one who can
	  allocate in them. The rest of the threads would only deallocate objects.
  */
  class ThreadFixedSizeMemoryPool
    {
    public:
	  ///@name Type Definitions
	  ///@{

		using SizeType = Chunk::SizeType;

		using ChunkList = std::list<Chunk*>;

	  ///@}

		static constexpr SizeType MaximumEmptyChunksToKeep = 4;
		static constexpr SizeType MaximumPointersToKeep = 1024;

	  ///@name Life Cycle
      ///@{

      /// Default constructor is deleted.
      ThreadFixedSizeMemoryPool() = delete;

	  /// Copy constructor is deleted.
	  ThreadFixedSizeMemoryPool(ThreadFixedSizeMemoryPool const& rOther) = delete;

	  /// Move constructor to be used in STL containers.
	  ThreadFixedSizeMemoryPool(ThreadFixedSizeMemoryPool&& rOther) = default;

	  /// The constructor to be called
	  ThreadFixedSizeMemoryPool(std::size_t BlockSizeInBytes, SizeType ChunkSize, std::size_t ThreadNumber)
		  : mBlockSizeInBytes(BlockSizeInBytes)
		  , mChunkSize(ChunkSize)
		  , mThreadNumber(ThreadNumber)
		  , mChunks()
		  , mAvailableChunks()
		  , mNumberOfReleasedChunks(0)
		  , mpCurrentChunk(nullptr)
	  {
	  }

      /// Destructor
	  virtual ~ThreadFixedSizeMemoryPool() {

	  }

      ///@}
      ///@name Operators
      ///@{

	  /// Assignment operator is deleted.
	  ThreadFixedSizeMemoryPool& operator=(ThreadFixedSizeMemoryPool const& rOther) = delete;

      ///@}
      ///@name Operations
      ///@{

	/// This function does not throw and returns zero if cannot allocate
	void* Allocate() {
		void* p_result = nullptr;
		if(mAvailablePointers.try_dequeue(p_result))
			return p_result;

		if(mpCurrentChunk == nullptr){
			AddChunk();
		}

		if (mpCurrentChunk->IsFull())
			if (!mAvailableChunks.try_dequeue(mpCurrentChunk))
				AddChunk();


		if (mpCurrentChunk->IsReleased())
			mpCurrentChunk->Initialize();
		KRATOS_CHECK_IS_FALSE(mpCurrentChunk->IsFull());
		p_result = mpCurrentChunk->Allocate();
		KRATOS_DEBUG_CHECK_NOT_EQUAL(p_result, nullptr);


		return p_result;
	}

	  bool Deallocate(void* pPointerToRelease) {
		  if(mAvailablePointers.size_approx() < MaximumPointersToKeep){
			  mAvailablePointers.enqueue(pPointerToRelease);
			  return true;
		  }


			  if (mpCurrentChunk->Has(pPointerToRelease)) {
				  DeallocateFromAvailableChunk(pPointerToRelease, mpCurrentChunk);
				  return true;
			  }

		  for (auto i_chunk = mChunks.begin(); i_chunk != mChunks.end(); i_chunk++)
		  {
			  if (i_chunk->Has(pPointerToRelease)) {
				  if (i_chunk->IsFull())
					DeallocateFromFullChunk(pPointerToRelease, &(*i_chunk));
				  else {
					  DeallocateFromAvailableChunk(pPointerToRelease, &(*i_chunk));
				  }
				  return true;
			  }
		  }


		  return false;
	  }

	  void Release() {
		  mChunks.clear();
		  mNumberOfReleasedChunks += mChunks.size();
		  Chunk* p_dummy;
		  while(mAvailableChunks.try_dequeue(p_dummy));
	  }

	  SizeType ChunkSize() const {
		  return mChunkSize;
	  }

	  std::size_t MemoryUsed() const {
		  std::size_t memory_used = sizeof(ThreadFixedSizeMemoryPool) + (mChunks.size() * 2 * sizeof(std::size_t));
		  for (auto i_chunk = mChunks.begin(); i_chunk != mChunks.end(); i_chunk++)
			  memory_used += i_chunk->MemoryUsed();
		  return memory_used;
	  }

	  std::size_t MemoryOverhead() const {
		  std::size_t memory_overhead = sizeof(ThreadFixedSizeMemoryPool) + (mAvailableChunks.size_approx() * sizeof(std::size_t));
		  for (auto i_chunk = mChunks.begin(); i_chunk != mChunks.end(); i_chunk++)
			  memory_overhead += i_chunk->MemoryOverhead();
		  return memory_overhead;
	  }

	  void lock(){
		while(std::atomic_flag_test_and_set_explicit(&ThreadFixedSizeMemoryPoolLock, std::memory_order_acquire))
             ; // spin until the lock is acquired
	  }

	  void unlock(){
		std::atomic_flag_clear_explicit(&ThreadFixedSizeMemoryPoolLock, std::memory_order_release);
	  }

	  ///@}
      ///@name Access
      ///@{

	  std::size_t GetNumberOfChunks() const {
		  return mChunks.size();
	  }

	  std::size_t GetNumberOfAvailableChunks() const {
		  return mAvailableChunks.size_approx();
	  }

	  std::size_t GetNumberOfReleasedChunks() const {
		  return mNumberOfReleasedChunks;
	  }

	  ///@}
      ///@name Inquiry
      ///@{

	  bool HasAvailableChunk() {
		  return !(mAvailableChunks.size_approx() == 0);
	  }


      ///@}
      ///@name Input and output
      ///@{

      /// Turn back information as a string.
	  std::string Info() const {
		  return "ThreadFixedSizeMemoryPool";
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

		  rOStream << GetNumberOfChunks() << " Chunks of "
			  << SizeInBytesToString(ChunkSize()) << " bytes each. Total memory usage: "
			  << SizeInBytesToString(MemoryUsed()) << " bytes and memory overhead "
			  << SizeInBytesToString(MemoryOverhead()) << "(" << overhead_percentage << "%)" << std::endl;
	  }

      ///@}


    private:
      ///@name Member Variables
      ///@{

		std::size_t mBlockSizeInBytes;
		SizeType mChunkSize;
		std::size_t mThreadNumber;
		std::list<Chunk> mChunks;
		moodycamel::ConcurrentQueue<Chunk*> mAvailableChunks;
		std::size_t mNumberOfReleasedChunks;
		moodycamel::ConcurrentQueue<void*> mAvailablePointers;
		Chunk* mpCurrentChunk= nullptr;


      ///@}
      ///@name Operations
      ///@{

		void AddChunk() {
			if (mThreadNumber != static_cast<std::size_t>(OpenMPUtils::ThisThread()))
                KRATOS_ERROR << "Trying to add chunk in thread " << mThreadNumber << " pool by thread " << static_cast<std::size_t>(OpenMPUtils::ThisThread());

			KRATOS_DEBUG_CHECK_EQUAL(mThreadNumber, static_cast<std::size_t>(OpenMPUtils::ThisThread()));

            mChunks.emplace_back(mBlockSizeInBytes, mChunkSize);
            mpCurrentChunk = &(mChunks.back());
            mpCurrentChunk->Initialize();
			// std::cout << "crating " << *p_available_chunk << std::endl;
		}

		void DeallocateFromAvailableChunk(void* pPointrerToRelease, Chunk* pChunk) {
			pChunk->Deallocate(pPointrerToRelease);
			if (pChunk->IsEmpty())
				if (mAvailableChunks.size_approx() - mNumberOfReleasedChunks > MaximumEmptyChunksToKeep)
					ReleaseChunk(pChunk);
		}

		void DeallocateFromFullChunk(void* pPointrerToRelease, Chunk* pChunk) {
			// It will be available after deallocating but is not in the list yet
			mAvailableChunks.enqueue(pChunk);
			pChunk->Deallocate(pPointrerToRelease);
			if (pChunk->IsEmpty()) // a rare case where a chunk has only one block! simptom of bad configuration
				if (mAvailableChunks.size_approx() - mNumberOfReleasedChunks > MaximumEmptyChunksToKeep)
					ReleaseChunk(pChunk);
		}

		void ReleaseChunk(Chunk* pChunk) {
			pChunk->Release();
			// mAvailableChunks.splice(mAvailableChunks.end(), mAvailableChunks, iChunk);
			mNumberOfReleasedChunks++;
		}

		// std::list<Chunk*> const& GetAvailableChunks() {
		// 	return mAvailableChunks;
		// }

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

    }; // Class ThreadFixedSizeMemoryPool

  ///@}
  ///@name Input and output
  ///@{


  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const ThreadFixedSizeMemoryPool& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.
