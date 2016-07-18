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

#ifdef _OPENMP
#include <omp.h>
#endif


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
  class Chunk
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

	  Chunk(Chunk&& rOther) 
		  : mpData(rOther.mpData)
		  , mSize(rOther.mSize)
		  , mBlockSizeInBytes(rOther.mBlockSizeInBytes)
		  , mThreadNumber(rOther.mThreadNumber)
#ifdef _OPENMP
		  , mLock(rOther.mLock)
#endif
	  {
		  rOther.mpData = nullptr;
	  }

	  /// The constructor to be called
	  Chunk(std::size_t BlockSizeInBytes, SizeType SizeInBytes)
		  : mSize(SizeInBytes)
		  , mBlockSizeInBytes(BlockSizeInBytes){
#ifdef _OPENMP
		  omp_init_lock(&mLock);
#endif
		  mThreadNumber = GetThreadNumber();
		  Initialize();
	  }

      /// Destructor is not virtual. This class can not be drived.
	  ~Chunk() {
			  delete[] mpData;
#ifdef _OPENMP
			  omp_destroy_lock(&mLock);
#endif
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
		  KRATOS_DEBUG_CHECK_EQUAL(GetThreadNumber(), mThreadNumber);// Initialize should be called only by owner thread!
		  std::size_t block_size_after_alignment = GetBlockSize(mBlockSizeInBytes);
		  mpData = new DataType[DataSize()];
  		  SetFirstAvailableBlockIndex(0); 
		  for (auto i_thread = 0; i_thread < GetMaxThreads(); i_thread++)
			  SetNumberOfAvailableBlocks(0, i_thread);
		  SetNumberOfAvailableBlocks(AllocatableDataSize() / block_size_after_alignment);

		  SizeType i_block = 0;
		  for (auto p = GetData(); i_block < GetNumberOfBlocks(); p += block_size_after_alignment)
			  *p = ++i_block;
	  }

	  /// This function does not throw and returns zero if cannot allocate
	  void* Allocate() {
		  KRATOS_DEBUG_CHECK_EQUAL(mThreadNumber, GetThreadNumber()); // Allocate should be called only by owner thread!
		  if (GetTotalNumberOfAvailableBlocks() == 0)
			  return nullptr;

		  if (GetNumberOfAvailableBlocks() == 0) // Time to get blocks from other threads
		  {
			  SetLock();
			  for (auto i_thread = 0; i_thread < GetMaxThreads(); i_thread++)
				  if (GetNumberOfAvailableBlocks(i_thread) != 0) {
					  SetFirstAvailableBlockIndex(GetFirstAvailableBlockIndex(i_thread));
					  SetNumberOfAvailableBlocks(GetNumberOfAvailableBlocks(i_thread));
					  SetNumberOfAvailableBlocks(0, i_thread);
					  break;
				  }
			  UnSetLock();
		  }
		  DataType* p_data = GetData();
		  auto first_index = GetFirstAvailableBlockIndex();
		  auto block_size = GetBlockSize(mBlockSizeInBytes);
		  auto thread_number = GetThreadNumber();
		  DataType * p_result = GetData() + (GetFirstAvailableBlockIndex() * GetBlockSize(mBlockSizeInBytes));
		  KRATOS_DEBUG_CHECK(Has(p_result));
		  SetFirstAvailableBlockIndex(*p_result);
		  SetNumberOfAvailableBlocks(GetNumberOfAvailableBlocks()-1);


		  return p_result;
	  }

	  void Deallocate(void* pPointrerToRelease) {

		  KRATOS_DEBUG_CHECK_NOT_EQUAL(GetData(), nullptr);
		  // Range check at least in lower bound.
		  KRATOS_DEBUG_CHECK_GREATER_EQUAL(pPointrerToRelease, GetData());
		  DataType* p_to_release = static_cast<DataType*>(pPointrerToRelease);

		  // Alignment check
		  KRATOS_DEBUG_CHECK_EQUAL((p_to_release - GetData()) % GetBlockSize(mBlockSizeInBytes), 0);

		  *p_to_release = GetFirstAvailableBlockIndex();
		  SetFirstAvailableBlockIndex(static_cast<SizeType>((p_to_release - GetData()) / GetBlockSize(mBlockSizeInBytes)));

		  // Check if there is no truncation error
		  KRATOS_DEBUG_CHECK_EQUAL(GetFirstAvailableBlockIndex(), double(p_to_release - GetData()) / GetBlockSize(mBlockSizeInBytes));
		  SetNumberOfAvailableBlocks(GetNumberOfAvailableBlocks() + 1);
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
		  return (mpData+ GetMaxThreads())[GetThreadNumber()];
	  }


#ifdef _OPENMP
	  omp_lock_t& GetLock() const
	  {
		  return mLock;
	  }
#endif

	  inline void SetLock() const
	  {
		  //does nothing if openMP is not present
#ifdef _OPENMP
		  omp_set_lock(&mLock);
#endif
	  }

	  inline void UnSetLock() const
	  {
		  //does nothing if openMP is not present
#ifdef _OPENMP
		  omp_unset_lock(&mLock);
#endif
	  }

	  ///@}
      ///@name Inquiry
      ///@{

	  bool HasAvailableBlock() const {
		  if (GetNumberOfAvailableBlocks() != 0)
			  return true;
#ifdef _OPENMP
		  if (mThreadNumber == GetThreadNumber()) // The thread which chunk belongs to
		  {
			  SetLock();
			  for (auto i_thread = 0; i_thread < GetMaxThreads(); i_thread++)
				  if (GetNumberOfAvailableBlocks(i_thread) != 0)
				  {
					  UnSetLock();
					  return true;
				  }
			  UnSetLock();
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
		SizeType mThreadNumber;
#ifdef _OPENMP
		mutable omp_lock_t mLock;
#endif

      ///@}
      ///@name Operations
      ///@{

	  /**
	  @return Maximum number of OpenMP threads that will be used in
	  parallel regions.
	  */
		static inline int GetMaxThreads()
		{
#ifdef _OPENMP
			return omp_get_max_threads();
#else
			return 1;
#endif
		}


		/**
		@return the number of OpenMP threads that will be used in
		parallel regions.
		*/
		static inline int GetThreadNumber()
		{
#ifdef _OPENMP
			return omp_get_thread_num();
#else
			return 0;
#endif
		}
		SizeType GetHeaderSize() const {
			return GetMaxThreads() * 2; // mFirstAvailableBlockIndex and mNumberOfAvailableBlocks for each thread
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
			mpData[GetThreadNumber()] = NewValue;
		}

		SizeType GetFirstAvailableBlockIndex() const {
			KRATOS_DEBUG_CHECK_NOT_EQUAL(mpData, nullptr);
			// remember that the first n blocks are the FirstAvailableBlockIndex for each n threads
			return mpData[GetThreadNumber()];
		}

		SizeType GetFirstAvailableBlockIndex(SizeType TheThreadNumber) const {
			KRATOS_DEBUG_CHECK_NOT_EQUAL(mpData, nullptr);
			// remember that the first n blocks are the FirstAvailableBlockIndex for each n threads
			return mpData[TheThreadNumber];
		}

		void SetNumberOfAvailableBlocks(SizeType NumberOfAvailableBlocks, SizeType TheThreadNumber) {
			KRATOS_DEBUG_CHECK_NOT_EQUAL(mpData, nullptr);
			// remember that the first n blocks are the FirstAvailableBlockIndex for each n threads
			(mpData + GetMaxThreads())[TheThreadNumber] = NumberOfAvailableBlocks;
		}

		void SetNumberOfAvailableBlocks(SizeType NumberOfAvailableBlocks) {
			KRATOS_DEBUG_CHECK_NOT_EQUAL(mpData, nullptr);
			// remember that the first n blocks are the FirstAvailableBlockIndex for each n threads
			(mpData + GetMaxThreads())[GetThreadNumber()] = NumberOfAvailableBlocks;
		}


		SizeType GetNumberOfAvailableBlocks(SizeType TheThreadNumber) const {
			KRATOS_DEBUG_CHECK_NOT_EQUAL(mpData, nullptr);
			// remember that the first n blocks are the FirstAvailableBlockIndex for each n threads
			return (mpData + GetMaxThreads())[TheThreadNumber];
		}

		SizeType GetTotalNumberOfAvailableBlocks() const {
			KRATOS_DEBUG_CHECK_NOT_EQUAL(mpData, nullptr);
			// remember that the first n blocks are the FirstAvailableBlockIndex for each n threads
			SizeType total_number_of_available_blocks = 0;
			SetLock();
			for (auto i_thread = 0; i_thread < GetMaxThreads(); i_thread++)
				total_number_of_available_blocks += (mpData + GetMaxThreads())[i_thread];
			UnSetLock();
			return total_number_of_available_blocks;
		}

		///@}

    }; // Class Chunk

  ///@}
  ///@name Input and output
  ///@{


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
