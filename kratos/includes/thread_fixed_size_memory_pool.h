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


#if !defined(KRATOS_THREAD_FIXED_SIZE_MEMORY_POOL_H_INCLUDED )
#define  KRATOS_THREAD_FIXED_SIZE_MEMORY_POOL_H_INCLUDED

#include <vector>
#include <list>
#include <algorithm>

#include "includes/chunk.h"

namespace Kratos
{
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
  class ThreadFixedSizeMemoryPool : public LockObject
    {
    public:
	  ///@name Type Definitions
	  ///@{

		using DataType = Chunk::DataType;

		using SizeType = DataType;

		using ChunkList = std::list<Chunk*>;

	  ///@}

		static constexpr SizeType MaximumEmptyChunksToKeep = 100000000;

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
		  : LockObject()
		  , mBlockSizeInBytes(BlockSizeInBytes)
		  , mChunkSize(ChunkSize)
		  , mThreadNumber(ThreadNumber)
		  , mChunks()
		  , mAvailableChunks()
		  , mNumberOfReleasedChunks(0)
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
		  if (mAvailableChunks.empty())
			  AddChunk();

		  void* p_result = nullptr;
		  Chunk& r_available_chunk = *(GetAvailableChunks().front());
		  if (r_available_chunk.IsReleased())
			  r_available_chunk.Initialize();
		  KRATOS_CHECK_IS_FALSE(r_available_chunk.IsFull());
		  p_result = r_available_chunk.Allocate();
		  KRATOS_DEBUG_CHECK_NOT_EQUAL(p_result, nullptr);
		  if (r_available_chunk.IsFull())
			  mAvailableChunks.pop_front();
		  return p_result;
	  }

	  bool Deallocate(void* pPointrerToRelease) {
		  if (!mAvailableChunks.empty()) {
			  auto i_chunk = mAvailableChunks.begin();
			  if ((*i_chunk)->Has(pPointrerToRelease)) {
				  DeallocateFromAvailableChunk(pPointrerToRelease, i_chunk);
				  return true;
			  }
		  }

		  for (auto i_chunk = mChunks.begin(); i_chunk != mChunks.end(); i_chunk++)
		  {
			  if (i_chunk->Has(pPointrerToRelease)) {
				  if (i_chunk->IsFull())
					DeallocateFromFullChunk(pPointrerToRelease, &(*i_chunk));
				  else {
					  auto i_available_chunk = std::find(mAvailableChunks.begin(), mAvailableChunks.end(), &(*i_chunk));
					  KRATOS_DEBUG_CHECK_NOT_EQUAL(i_available_chunk, mAvailableChunks.end()); // Un explicable!
					  DeallocateFromAvailableChunk(pPointrerToRelease, i_available_chunk);
				  }
				  return true;
			  }
		  }

		  return false;
	  }

	  void Release() {
		  mChunks.clear();
		  mNumberOfReleasedChunks += mChunks.size();
		  mAvailableChunks.clear();
		  for (auto i_chunk = mChunks.begin(); i_chunk != mChunks.end(); i_chunk++)
			  mAvailableChunks.push_front(&(*i_chunk));
	  }

	  SizeType ChunkSize() const {
		  return mChunkSize;
	  }

	  std::size_t MemoryUsed() const {
		  std::size_t memory_used = sizeof(ThreadFixedSizeMemoryPool) + (mAvailableChunks.size() * 2 * sizeof(std::size_t));
		  for (auto i_chunk = mChunks.begin(); i_chunk != mChunks.end(); i_chunk++)
			  memory_used += i_chunk->MemoryUsed();
		  return memory_used;
	  }

	  std::size_t MemoryOverhead() const {
		  std::size_t memory_overhead = sizeof(ThreadFixedSizeMemoryPool) + (mAvailableChunks.size() * sizeof(std::size_t));
		  for (auto i_chunk = mChunks.begin(); i_chunk != mChunks.end(); i_chunk++)
			  memory_overhead += i_chunk->MemoryOverhead();
		  return memory_overhead;
	  }

	  ///@}
      ///@name Access
      ///@{

	  std::size_t GetNumberOfChunks() const {
		  return mChunks.size();
	  }

	  std::size_t GetNumberOfAvailableChunks() const {
		  return mAvailableChunks.size();
	  }

	  std::size_t GetNumberOfReleasedChunks() const {
		  return mNumberOfReleasedChunks;
	  }

	  ///@}
      ///@name Inquiry
      ///@{

	  bool HasAvailableChunk() {
		  return !(mAvailableChunks.empty());
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
		std::list<Chunk*> mAvailableChunks;
		std::size_t mNumberOfReleasedChunks;



      ///@}
      ///@name Operations
      ///@{

		void AddChunk() {
			if (mThreadNumber != static_cast<std::size_t>(GetThreadNumber()))
                KRATOS_ERROR;

			KRATOS_DEBUG_CHECK_EQUAL(mThreadNumber, static_cast<std::size_t>(GetThreadNumber()));

            mChunks.emplace_back(mBlockSizeInBytes, mChunkSize);
            Chunk* p_available_chunk = &(mChunks.back());
            p_available_chunk->Initialize();
            mAvailableChunks.push_front(p_available_chunk);
		}

		void DeallocateFromAvailableChunk(void* pPointrerToRelease, ChunkList::iterator iChunk) {
			(*iChunk)->Deallocate(pPointrerToRelease);
			if ((*iChunk)->IsEmpty())
				if (mAvailableChunks.size() - mNumberOfReleasedChunks > MaximumEmptyChunksToKeep)
					ReleaseChunk(iChunk);
		}

		void DeallocateFromFullChunk(void* pPointrerToRelease, Chunk* pChunk) {
			// It will be available after deallocating but is not in the list yet
			mAvailableChunks.push_front(pChunk);
			auto i_chunk = mAvailableChunks.begin();
			(*i_chunk)->Deallocate(pPointrerToRelease);
			if ((*i_chunk)->IsEmpty()) // a rare case where a chunk has only one block! simptom of bad configuration
				if (mAvailableChunks.size() - mNumberOfReleasedChunks > MaximumEmptyChunksToKeep)
					ReleaseChunk(i_chunk);
		}

		void ReleaseChunk(ChunkList::iterator iChunk) {
			//(*iChunk)->Release();
			//mAvailableChunks.splice(mAvailableChunks.end(), mAvailableChunks, iChunk);
			//mNumberOfReleasedChunks++;
		}

		std::list<Chunk*> const& GetAvailableChunks() {
			return mAvailableChunks;
		}

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

#endif // KRATOS_THREAD_FIXED_SIZE_MEMORY_POOL_H_INCLUDED  defined
