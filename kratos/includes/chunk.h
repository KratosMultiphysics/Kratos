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
	  containing the chunk size and block size and thread private links.
  */
  class Chunk : public LockObject
    {
    public:

	  ///@name Type Definitions
	  ///@{

		// The reference algorithm uses unsigned char as data type
		// Here I use std::size_t to ensure better alignment considering
		// that objects in Kratos are "always" larger than a double
		using DataType = std::size_t;

		using SizeType = DataType;


	  ///@}
	  ///@name Life Cycle
      ///@{

      /// Default constructor is deleted.
      Chunk() = delete;

	  /// Copy constructor is deleted.
	  Chunk(Chunk const& rOther) = delete;

	  Chunk(Chunk&& rOther) noexcept
		  : LockObject(std::move(rOther))
		  , mpData(rOther.mpData)
		  , mSize(rOther.mSize)
		  , mBlockSizeInBytes(rOther.mBlockSizeInBytes)
		  , mOwnerThread(rOther.mOwnerThread)
	  {
		  rOther.mpData = nullptr;
	  }

	  /// The constructor to be called
	  Chunk(std::size_t BlockSizeInBytes, SizeType SizeInBytes) noexcept
		  : LockObject()
		  , mpData(nullptr)
		  , mSize(SizeInBytes)
		  , mBlockSizeInBytes(BlockSizeInBytes){
		  mOwnerThread = OpenMPUtils::ThisThread();
	  }

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
		  mOwnerThread = OpenMPUtils::ThisThread();  // initialization can change the owner thread
		  std::size_t block_size_after_alignment = GetBlockSize(mBlockSizeInBytes);
		  mpData = new DataType[DataSize()];
  		  SetFirstAvailableBlockIndex(0);
		  for (auto i_thread = 0; i_thread < OpenMPUtils::GetCurrentNumberOfThreads(); i_thread++)
			  SetNumberOfAvailableBlocks(0, i_thread);
		  SetNumberOfAvailableBlocks(AllocatableDataSize() / block_size_after_alignment);

		  SizeType i_block = 0;
		  for (auto p = GetData(); i_block < GetNumberOfBlocks(); p += block_size_after_alignment)
			  *p = ++i_block;

		  KRATOS_DEBUG_CHECK_EQUAL(i_block, GetNumberOfAvailableBlocks());
	  }

	  /// This function does not throw and returns zero if cannot allocate
	  void* Allocate() {
		  KRATOS_DEBUG_CHECK_EQUAL(mOwnerThread, OpenMPUtils::ThisThread()); // Allocate should be called only by owner thread
		  KRATOS_DEBUG_CHECK_NOT_EQUAL(mpData, nullptr);

		  if (GetTotalNumberOfAvailableBlocks() == 0)
			  return nullptr;

		  if (GetNumberOfAvailableBlocks() == 0) // Time to get blocks from other threads
		  {
			  lock();
			  for (auto i_thread = 0; i_thread < OpenMPUtils::GetCurrentNumberOfThreads(); i_thread++)
				  if (GetNumberOfAvailableBlocks(i_thread) != 0) {
					  SetFirstAvailableBlockIndex(GetFirstAvailableBlockIndex(i_thread));
					  SetNumberOfAvailableBlocks(GetNumberOfAvailableBlocks(i_thread));
					  SetNumberOfAvailableBlocks(0, i_thread);
					  break;
				  }
			  unlock();
		  }
		  DataType * p_result = GetData() + (GetFirstAvailableBlockIndex() * GetBlockSize(mBlockSizeInBytes));
		  KRATOS_DEBUG_CHECK(Has(p_result));
		  SetFirstAvailableBlockIndex(*p_result);
		  SetNumberOfAvailableBlocks(GetNumberOfAvailableBlocks()-1);


		  return p_result;
	  }

	  void Deallocate(void* pPointrerToRelease) {
		  if (pPointrerToRelease == nullptr)
			  return;

		  KRATOS_DEBUG_CHECK_NOT_EQUAL(GetData(), nullptr);
		  // Range check at least in lower bound.
		  KRATOS_DEBUG_CHECK_GREATER_EQUAL(pPointrerToRelease, GetData());
		  DataType* p_to_release = static_cast<DataType*>(pPointrerToRelease);

		  // Alignment check
		  KRATOS_DEBUG_CHECK_EQUAL((p_to_release - GetData()) % GetBlockSize(mBlockSizeInBytes), 0);
		  lock();
		  *p_to_release = GetFirstAvailableBlockIndex();
		  SetFirstAvailableBlockIndex(static_cast<SizeType>((p_to_release - GetData()) / GetBlockSize(mBlockSizeInBytes)));

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
		  return mSize; // GetBlockSize(BlockSizeInBytes) * NumberOfBlocks * sizeof(DataType);
	  }

	  std::size_t MemoryUsed() const {
		  if (mpData == nullptr)
			  return sizeof(Chunk);
		  return Size() + sizeof(Chunk);
	  }

	  std::size_t MemoryOverhead() const {
		  if (mpData == nullptr)
			  return sizeof(Chunk);
		  return MemoryUsed() - (mBlockSizeInBytes*(GetNumberOfBlocks() - GetTotalNumberOfAvailableBlocks()));
	  }

	  ///@}
      ///@name Access
      ///@{

	  SizeType GetNumberOfBlocks() const {
		  return AllocatableDataSize() / GetBlockSize(mBlockSizeInBytes);
	  }

	  SizeType GetNumberOfAvailableBlocks() const {
		  KRATOS_DEBUG_CHECK_NOT_EQUAL(mpData, nullptr);
		  // remember that the first n blocks are the FirstAvailableBlockIndex for each n threads
		  return (mpData+ OpenMPUtils::GetCurrentNumberOfThreads())[OpenMPUtils::ThisThread()];
	  }

	  const DataType* pGetData() const {
		  return mpData;
	  }


	  ///@}
      ///@name Inquiry
      ///@{

	  bool HasAvailableBlock() const {
		  if (GetNumberOfAvailableBlocks() != 0)
			  return true;
#ifdef _OPENMP
		  if (mOwnerThread == OpenMPUtils::ThisThread()) // The thread which chunk belongs to
		  {
			  lock();
			  for (auto i_thread = 0; i_thread < OpenMPUtils::GetCurrentNumberOfThreads(); i_thread++)
				  if (GetNumberOfAvailableBlocks(i_thread) != 0)
				  {
					  unlock();
					  return true;
				  }
			  unlock();
		  }
#endif
		  return false;
	  }

	  bool Has(void* pThePointer) const {
		  return ((pThePointer >= GetData()) && (pThePointer <= (mpData + Size() - mBlockSizeInBytes)));
	  }

	  bool IsEmpty() const {
		  return (GetTotalNumberOfAvailableBlocks() == GetNumberOfBlocks());
	  }

	  bool IsFull() {
		  return (GetTotalNumberOfAvailableBlocks() == 0);
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
		  rOStream << static_cast<SizeType>(GetTotalNumberOfAvailableBlocks()) << " blocks are available" << std::endl;
	  }

      ///@}


    private:
      ///@name Member Variables
      ///@{

		DataType * mpData;
		SizeType mSize;
		SizeType mBlockSizeInBytes;
		int mOwnerThread;

      ///@}
      ///@name Operations
      ///@{

		SizeType GetHeaderSize() const {
			return OpenMPUtils::GetCurrentNumberOfThreads() * 2; // mFirstAvailableBlockIndex and mNumberOfAvailableBlocks for each thread
		}

		const DataType* GetData() const {
			if (mpData == nullptr)
				return nullptr;
			return mpData + GetHeaderSize();
		}

		DataType* GetData() {
			if (mpData == nullptr)
				return nullptr;
			return mpData + GetHeaderSize();
		}

		SizeType DataSize() const {
			return mSize / sizeof(DataType);
		}

		SizeType AllocatableDataSize() const {
			std::size_t header_size = GetHeaderSize();
			std::size_t data_size = DataSize();
			if(data_size > header_size)
				return data_size - header_size;
			return 0;
		}

		static std::size_t GetBlockSize(std::size_t BlockSizeInBytes) {
			return (BlockSizeInBytes + sizeof(DataType) - 1) / sizeof(DataType);
		}

		void SetFirstAvailableBlockIndex(SizeType NewValue) {
			KRATOS_DEBUG_CHECK_NOT_EQUAL(mpData, nullptr);
			// remember that the first n blocks are the FirstAvailableBlockIndex for each n threads
			mpData[OpenMPUtils::ThisThread()] = NewValue;
		}

		SizeType GetFirstAvailableBlockIndex() const {
			KRATOS_DEBUG_CHECK_NOT_EQUAL(mpData, nullptr);
			// remember that the first n blocks are the FirstAvailableBlockIndex for each n threads
			return mpData[OpenMPUtils::ThisThread()];
		}

		SizeType GetFirstAvailableBlockIndex(SizeType TheThreadNumber) const {
			KRATOS_DEBUG_CHECK_NOT_EQUAL(mpData, nullptr);
			// remember that the first n blocks are the FirstAvailableBlockIndex for each n threads
			return mpData[TheThreadNumber];
		}

		void SetNumberOfAvailableBlocks(SizeType NumberOfAvailableBlocks, SizeType TheThreadNumber) {
			KRATOS_DEBUG_CHECK_NOT_EQUAL(mpData, nullptr);
			// remember that the first n blocks are the FirstAvailableBlockIndex for each n threads
			(mpData + OpenMPUtils::GetCurrentNumberOfThreads())[TheThreadNumber] = NumberOfAvailableBlocks;
		}

		void SetNumberOfAvailableBlocks(SizeType NumberOfAvailableBlocks) {
			KRATOS_DEBUG_CHECK_NOT_EQUAL(mpData, nullptr);
			// remember that the first n blocks are the FirstAvailableBlockIndex for each n threads
			(mpData + OpenMPUtils::GetCurrentNumberOfThreads())[OpenMPUtils::ThisThread()] = NumberOfAvailableBlocks;
		}


		SizeType GetNumberOfAvailableBlocks(SizeType TheThreadNumber) const {
			KRATOS_DEBUG_CHECK_NOT_EQUAL(mpData, nullptr);
			// remember that the first n blocks are the FirstAvailableBlockIndex for each n threads
			return (mpData + OpenMPUtils::GetCurrentNumberOfThreads())[TheThreadNumber];
		}

		SizeType GetTotalNumberOfAvailableBlocks() const {
			KRATOS_DEBUG_CHECK_NOT_EQUAL(mpData, nullptr);
			// remember that the first n blocks are the FirstAvailableBlockIndex for each n threads
			SizeType total_number_of_available_blocks = 0;
			lock();
			for (auto i_thread = 0; i_thread < OpenMPUtils::GetCurrentNumberOfThreads(); i_thread++)
				total_number_of_available_blocks += (mpData + OpenMPUtils::GetCurrentNumberOfThreads())[i_thread];
			unlock();
			return total_number_of_available_blocks;
		}

		///@}

    }; // Class Chunk

  ///@}
  ///@name Operators
  ///@{

	inline bool operator < (Chunk const& rFirst, Chunk const& rSedond) {
		return rFirst.pGetData() < rSedond.pGetData();
	}

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const Chunk& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CHUNK_H_INCLUDED  defined
