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

#include "includes/chunk.h"

namespace Kratos
{
  ///@addtogroup KratosCore
  ///@{

  ///@name Kratos Classes
  ///@{

  /// FixedSizeMemoryPool is the smallest building block of Kratos memory management.
  /** The memory management of Kratos is implemented based on the design
	  given in Modern C++ Design by A. Alexandrescu and FixedSizeMemoryPool is the second
	  layer of it "AKA FixedAllocator" holding chunks.
  */
  class FixedSizeMemoryPool
    {
    public:    
	  ///@name Type Definitions
	  ///@{

		using DataType = Chunk::DataType;

		using SizeType = DataType;

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
		  : mBlockSizeInBytes(BlockSizeInBytes)
		  , mChunkSize(ChunkSize)
	  {
		  AddChunk();
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
		  void* p_result = nullptr;
#pragma omp critical
		  {
			  if (mAvailableChunksIndices.empty())
				  AddChunk();

			  Chunk& r_available_chunk = mChunks[mAvailableChunksIndices.back()];
			  if (r_available_chunk.IsReleased())
				  r_available_chunk.Initialize();
			  p_result = r_available_chunk.Allocate();
			  if (r_available_chunk.IsFull()) mAvailableChunksIndices.pop_back();
		  }
		  return p_result;
	  }

	  void Deallocate(void* pPointrerToRelease) {
		  if (!mAvailableChunksIndices.empty()) {
			  std::size_t chunk_index = mAvailableChunksIndices.back();
			  if (mChunks[chunk_index].Has(pPointrerToRelease)) {
				  DeallocateFromChunk(pPointrerToRelease, chunk_index);
				  return;
			  }
		  }

		  for (std::size_t i_chunk_index = 0; i_chunk_index < mChunks.size(); i_chunk_index++)
		  {
			  if (mChunks[i_chunk_index].Has(pPointrerToRelease)) {
				  DeallocateFromChunk(pPointrerToRelease, i_chunk_index);
				  return;
			  }
		  }

		  KRATOS_ERROR << "The Pointer with address " << pPointrerToRelease << " was not found in this pool" << std::endl;
	  }

	  void ReleaseChunk(std::size_t ChunkIndex) {
		  if(mAvailableChunksIndices.size() > 1)
			mChunks[ChunkIndex].Release();
	  }

	  void Release() {
		  mChunks.clear();
	  }

	  SizeType ChunkSize() const {
		  return mChunkSize;
	  }

	  std::size_t MemoryUsed() const {
		  std::size_t memory_used = sizeof(FixedSizeMemoryPool) + (mAvailableChunksIndices.size() * sizeof(std::size_t));
		  for (auto i_chunk = mChunks.begin(); i_chunk != mChunks.end(); i_chunk++)
			  memory_used += i_chunk->MemoryUsed();
		  return memory_used;
	  }

	  std::size_t MemoryOverhead() const {
		  std::size_t memory_overhead = sizeof(FixedSizeMemoryPool) + (mAvailableChunksIndices.size() * sizeof(std::size_t));
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
		std::vector<Chunk> mChunks;
		std::vector<std::size_t> mAvailableChunksIndices;

      ///@}
      ///@name Operations
      ///@{

		void AddChunk() {
			mChunks.push_back(Chunk(mBlockSizeInBytes, mChunkSize));
			mAvailableChunksIndices.push_back(mChunks.size() - 1);
		}

		void DeallocateFromChunk(void* pPointrerToRelease, std::size_t ChunkIndex) {
			Chunk& r_chunk = mChunks[ChunkIndex];
			if (r_chunk.IsFull()) // It will be available after deallocating but is not in the list yet
				mAvailableChunksIndices.push_back(ChunkIndex);
			r_chunk.Deallocate(pPointrerToRelease);
			if (r_chunk.IsEmpty())
				ReleaseChunk(ChunkIndex);
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
