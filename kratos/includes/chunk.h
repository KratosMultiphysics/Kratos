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


#if !defined(KRATOS_CHUNK_H_INCLUDED )
#define  KRATOS_CHUNK_H_INCLUDED

#include <atomic>

#include "includes/lock_object.h"
#include "utilities/openmp_utils.h"

namespace Kratos
{
  ///@addtogroup KratosCore
  ///@{

  ///@name Kratos Classes
  ///@{

  /// Chunk is the smallest building block of Kratos memory management.
  /** The memory management of Kratos is implemented based on the design
	  given in Modern C++ Design by A. Alexandrescu and chunk is the lower
	  layer of it holding a chunk of NumberOfBlocks objects of BlockSize.
	  This implemnetation is designed for larg chunk size (i.e. 1M)
	  imposes more overhead than the reference for stroing a header
	  containing the chunk size and block size.
  */
  class Chunk 
    {
    public:

	  ///@name Type Definitions
	  ///@{

		// The reference algorithm uses unsigned char as data type
		// Here I use std::int64_t to ensure better alignment considering
		// that objects in Kratos are "always" larger than a double
		using BlockType = std::int64_t;

		using SizeType = std::size_t;


	  ///@}
	  ///@name Life Cycle
      ///@{

      /// Default constructor is deleted.
      Chunk() = delete;

	  /// Copy constructor is deleted.
	  Chunk(Chunk const& rOther) = delete;

	  Chunk(Chunk&& rOther) = delete;

	  /// The constructor to be called
	  Chunk(std::size_t BlockSizeInBytes, SizeType SizeInBytes) noexcept
		  : mpData(nullptr)
		  , mpEnd(nullptr)
		  , mpUninitializedMemory(nullptr)
		  , mSize(SizeInBytes)
		  , mBlockSizeInBytes(BlockSizeInBytes)
		  , mNumberOfAvailableBlocks(0) // to be initialized in initialize
		  , mFirstAvailableBlock(mpData){}

      /// Destructor is not virtual. This class can not be drived.
	  ~Chunk() noexcept {
			  delete[] mpData;
	  }

      ///@}
      ///@name Operators
      ///@{

	  /// Assignment operator is deleted.
	  Chunk& operator=(Chunk const& rOther) = delete;

      ///@}
      ///@name Operations
      ///@{
	  void Initialize() {
		  const std::size_t block_size_after_alignment = GetBlockSize(mBlockSizeInBytes);
		  const std::size_t data_size = DataSize();
		  mpData = new BlockType[data_size];
		  mpEnd = mpData + data_size;
  		  mFirstAvailableBlock = mpData;
		  mNumberOfAvailableBlocks = AllocatableDataSize() / block_size_after_alignment;
		  mpUninitializedMemory = mpData;
		  

		  *mpData = 0; // The first entry of the link list to the next one
	  }

	  /// This function does not throw and returns zero if cannot allocate
	  void* Allocate() {
		  KRATOS_DEBUG_CHECK_NOT_EQUAL(mpData, nullptr);

		  if (mNumberOfAvailableBlocks == 0)
			  return nullptr;

		    const std::size_t block_size_after_alignment = GetBlockSize(mBlockSizeInBytes);
			lock();  
			BlockType * p_result = mFirstAvailableBlock;
			if(p_result >= mpUninitializedMemory) {
				mFirstAvailableBlock += block_size_after_alignment;
				if(mpUninitializedMemory < mpEnd) {
					mpUninitializedMemory += block_size_after_alignment;
				}
			}
				
			KRATOS_DEBUG_CHECK(Has(p_result));
			unlock();

			mNumberOfAvailableBlocks--;
	
		  return p_result;
	  }

	  void Deallocate(void* pPointrerToRelease) {
		  if (pPointrerToRelease == nullptr)
			  return;

		  KRATOS_DEBUG_CHECK_NOT_EQUAL(mpData, nullptr);
		  // Range check at least in lower bound.
		  KRATOS_DEBUG_CHECK_GREATER_EQUAL(pPointrerToRelease, mpData);
		  BlockType* p_to_release = static_cast<BlockType*>(pPointrerToRelease);

		  // Alignment check
		  KRATOS_DEBUG_CHECK_EQUAL((p_to_release - mpData) % GetBlockSize(mBlockSizeInBytes), 0);
		  
		  lock();
		  *p_to_release = (BlockType)mFirstAvailableBlock.exchange(p_to_release);
		  unlock();

		  mNumberOfAvailableBlocks++;
		  pPointrerToRelease = nullptr;

	  }

	  void Release() {
		  if (mpData == nullptr)
			  return;
		  delete[] mpData;
		  mpData = nullptr;

	  }

	  std::size_t Size() const {
		  return mSize; // GetBlockSize(BlockSizeInBytes) * NumberOfBlocks * sizeof(BlockType);
	  }

	  std::size_t MemoryUsed() const {
		  if (mpData == nullptr)
			  return sizeof(Chunk);
		  return Size() + sizeof(Chunk);
	  }

	  std::size_t MemoryOverhead() const {
		  if (mpData == nullptr)
			  return sizeof(Chunk);
		  return MemoryUsed() - (mBlockSizeInBytes*(GetNumberOfBlocks() - mNumberOfAvailableBlocks));
	  }

	  ///@}
      ///@name Access
      ///@{

	  SizeType GetNumberOfBlocks() const {
		  return AllocatableDataSize() / GetBlockSize(mBlockSizeInBytes);
	  }

	  SizeType GetNumberOfAvailableBlocks() const {
		  return mNumberOfAvailableBlocks;
	  }

	  const BlockType* pGetData() const {
		  return mpData;
	  }

	  const BlockType* pDataBegin() const {
		  return mpData;
	  }

	  const BlockType* pDataEnd() const {
		  return mpData + DataSize();
	  }

	  ///@}
      ///@name Inquiry
      ///@{

	  bool HasAvailableBlock() const {
		  if (GetNumberOfAvailableBlocks() != 0)
			  return true;

		  return false;
	  }

	  bool Has(const void* pThePointer) const {
		  return ((pThePointer >= mpData) && (pThePointer <= (mpData + DataSize() - GetBlockSize(mBlockSizeInBytes))));
	  }

	  bool IsEmpty() const {
		  return (mNumberOfAvailableBlocks == GetNumberOfBlocks());
	  }

	  bool IsFull() {
		  return (mNumberOfAvailableBlocks == 0);
	  }

	  bool IsInitialized() const {
		  return mpData != nullptr;
	  }

	  bool IsReleased() const {
		  return mpData == nullptr;
	  }

	  void lock(){
		while(std::atomic_flag_test_and_set_explicit(&mLocked, std::memory_order_acquire))
             ; // spin until the lock is acquired
	  }

	  void unlock(){
		std::atomic_flag_clear_explicit(&mLocked, std::memory_order_release);
	  }

	  ///@}
      ///@name Input and output
      ///@{

      /// Turn back information as a string.
	  std::string Info() const {
		  return "Chunk";
	  }

      /// Print information about this object.
      void PrintInfo(std::ostream& rOStream) const {
		  rOStream << Info();
	  }

      /// Print object's data.
	  void PrintData(std::ostream& rOStream) const {
		  rOStream << " from " << pDataBegin() << " to " << pDataEnd() << " with size " << Size() << " bytes allocated as " << DataSize() << " size_t with " << static_cast<SizeType>(GetNumberOfAvailableBlocks()) << " available blocks" << std::endl;
	  }

      ///@}


    private:
      ///@name Member Variables
      ///@{

		BlockType * mpData;
		BlockType* mpEnd;
		BlockType* mpUninitializedMemory;
		SizeType mSize;
		SizeType mBlockSizeInBytes;
		std::atomic_uint mNumberOfAvailableBlocks;
		std::atomic<BlockType*> mFirstAvailableBlock;
		std::atomic_flag mLocked = ATOMIC_FLAG_INIT;

      ///@}
      ///@name Operations
      ///@{

		SizeType GetHeaderSize() const {
			return 0; // mFirstAvailableBlockIndex
		}

		SizeType DataSize() const {
			return mSize / sizeof(BlockType);
		}

		SizeType AllocatableDataSize() const {
			std::size_t header_size = GetHeaderSize();
			std::size_t data_size = DataSize();
			if(data_size > header_size)
				return data_size - header_size;
			return 0;
		}

		static std::size_t GetBlockSize(std::size_t BlockSizeInBytes) {
			return (BlockSizeInBytes + sizeof(BlockType) - 1) / sizeof(BlockType);
		}

		///@}

    }; // Class Chunk

  ///@}
  ///@name Operators
  ///@{

	inline bool operator << (Chunk const& rFirst, Chunk const& rSedond) {
		return rFirst.pGetData() < rSedond.pGetData();
	}

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const Chunk& rThis)
    {
      rThis.PrintInfo(rOStream);
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CHUNK_H_INCLUDED  defined
