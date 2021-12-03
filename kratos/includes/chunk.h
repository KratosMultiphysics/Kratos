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
  class Chunk : public LockObject
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
		  : LockObject()
		  , mpData(nullptr)
		  , mSize(SizeInBytes)
		  , mBlockSizeInBytes(BlockSizeInBytes)
		  , mNumberOfAvailableBlocks(0) // to be initialized in initialize
		  , mFirstAvailableBlockIndex(0){}

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
		  std::size_t block_size_after_alignment = GetBlockSize(mBlockSizeInBytes);
		  mpData = new BlockType[DataSize()];
  		  mFirstAvailableBlockIndex = 0;
		  SetNumberOfAvailableBlocks(AllocatableDataSize() / block_size_after_alignment);

		  *mpData = -1;
	  }

	  /// This function does not throw and returns zero if cannot allocate
	  void* Allocate() {
		  KRATOS_DEBUG_CHECK_NOT_EQUAL(mpData, nullptr);

		  if (GetNumberOfAvailableBlocks() == 0)
			  return nullptr;

		    lock();
			  
			BlockType * p_result = GetData() + (GetFirstAvailableBlockIndex() * GetBlockSize(mBlockSizeInBytes));
			if(*p_result < 0) {
				*p_result = -(*p_result);
				if(GetNumberOfAvailableBlocks() > 1) {
					*(p_result + GetBlockSize(mBlockSizeInBytes)) = -(*p_result + 1);
				}
			}
				
			KRATOS_DEBUG_CHECK(Has(p_result));
			mFirstAvailableBlockIndex = *p_result;

			

			SetNumberOfAvailableBlocks(GetNumberOfAvailableBlocks()-1);
	
			unlock();

		  return p_result;
	  }

	  void Deallocate(void* pPointrerToRelease) {
		  if (pPointrerToRelease == nullptr)
			  return;

		  KRATOS_DEBUG_CHECK_NOT_EQUAL(GetData(), nullptr);
		  // Range check at least in lower bound.
		  KRATOS_DEBUG_CHECK_GREATER_EQUAL(pPointrerToRelease, GetData());
		  BlockType* p_to_release = static_cast<BlockType*>(pPointrerToRelease);

		  // Alignment check
		  KRATOS_DEBUG_CHECK_EQUAL((p_to_release - GetData()) % GetBlockSize(mBlockSizeInBytes), 0);
		  lock();
		  *p_to_release = GetFirstAvailableBlockIndex();
		  mFirstAvailableBlockIndex = static_cast<SizeType>((p_to_release - GetData()) / GetBlockSize(mBlockSizeInBytes));

		  // Check if there is no truncation error
		  KRATOS_DEBUG_CHECK_EQUAL(GetFirstAvailableBlockIndex(), double(p_to_release - GetData()) / GetBlockSize(mBlockSizeInBytes));
		  SetNumberOfAvailableBlocks(GetNumberOfAvailableBlocks() + 1);
		  unlock();
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
		  return MemoryUsed() - (mBlockSizeInBytes*(GetNumberOfBlocks() - GetNumberOfAvailableBlocks()));
	  }

	  ///@}
      ///@name Access
      ///@{

	  SizeType GetNumberOfBlocks() const {
		  return AllocatableDataSize() / GetBlockSize(mBlockSizeInBytes);
	  }

	  SizeType GetNumberOfAvailableBlocks() const {
		  KRATOS_DEBUG_CHECK_NOT_EQUAL(mpData, nullptr);
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
		  return ((pThePointer >= GetData()) && (pThePointer <= (mpData + DataSize() - GetBlockSize(mBlockSizeInBytes))));
	  }

	  bool IsEmpty() const {
		  return (GetNumberOfAvailableBlocks() == GetNumberOfBlocks());
	  }

	  bool IsFull() {
		  return (GetNumberOfAvailableBlocks() == 0);
	  }

	  bool IsInitialized() const {
		  return mpData != nullptr;
	  }

	  bool IsReleased() const {
		  return mpData == nullptr;
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
		SizeType mSize;
		SizeType mBlockSizeInBytes;
		SizeType mNumberOfAvailableBlocks;
		SizeType mFirstAvailableBlockIndex;

      ///@}
      ///@name Operations
      ///@{

		SizeType GetHeaderSize() const {
			return 0; // mFirstAvailableBlockIndex
		}

		const BlockType* GetData() const {
			if (mpData == nullptr)
				return nullptr;
			return mpData + GetHeaderSize();
		}

		BlockType* GetData() {
			if (mpData == nullptr)
				return nullptr;
			return mpData + GetHeaderSize();
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

		SizeType GetFirstAvailableBlockIndex() const {
			KRATOS_DEBUG_CHECK_NOT_EQUAL(mpData, nullptr);
			return mFirstAvailableBlockIndex;
		}

		void SetNumberOfAvailableBlocks(SizeType NumberOfAvailableBlocks) {
			KRATOS_DEBUG_CHECK_NOT_EQUAL(mpData, nullptr);
			mNumberOfAvailableBlocks = NumberOfAvailableBlocks;
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
