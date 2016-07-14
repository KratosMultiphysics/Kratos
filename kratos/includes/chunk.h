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
	  This implemnetation imposes more overhead than the reference for 
	  storing next chunk pointer in order to improve the performance.
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
		  , mFirstAvailableBlockIndex(rOther.mFirstAvailableBlockIndex)
		  , mNumberOfAvailableBlocks(rOther.mNumberOfAvailableBlocks)
	  {
		  rOther.mpData = nullptr;
	  }

	  /// The constructor to be called
	  Chunk(std::size_t BlockSizeInBytes, SizeType NumberOfBlocks) {
		  Initialize(BlockSizeInBytes, NumberOfBlocks);
	  }

      /// Destructor is not virtual. This class can not be drived.
	  ~Chunk() {
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
	  void Initialize(std::size_t BlockSizeInBytes, SizeType NumberOfBlocks) {
		  std::size_t block_size_after_alignment = GetBlockSize(BlockSizeInBytes);
		  mpData = new DataType[Size(BlockSizeInBytes, NumberOfBlocks)];
		  mFirstAvailableBlockIndex = 0;
		  mNumberOfAvailableBlocks = NumberOfBlocks;

		  SizeType i_block = 0;
		  for (auto p = mpData; i_block < NumberOfBlocks; p += block_size_after_alignment)
			  *p = ++i_block;
	  }

	  /// This function does not throw and returns zero if cannot allocate
	  void* Allocate(std::size_t BlockSizeInBytes) {
		  if (mNumberOfAvailableBlocks == 0)
			  return nullptr;

		  DataType * p_result = mpData + (mFirstAvailableBlockIndex * GetBlockSize(BlockSizeInBytes));
		  mFirstAvailableBlockIndex = *p_result;
		  mNumberOfAvailableBlocks--;

		  return p_result;
	  }

	  void Deallocate(void* pPointrerToRelease, std::size_t BlockSizeInBytes) {
		  KRATOS_DEBUG_CHECK_NOT_EQUAL(mpData, nullptr);
		  // Range check at least in lower bound.
		  KRATOS_DEBUG_CHECK_GREATER_EQUAL(pPointrerToRelease, mpData);
		  DataType* p_to_release = static_cast<DataType*>(pPointrerToRelease);

		  // Alignment check
		  KRATOS_DEBUG_CHECK_EQUAL((p_to_release - mpData) % GetBlockSize(BlockSizeInBytes), 0);

		  *p_to_release = mFirstAvailableBlockIndex;
		  mFirstAvailableBlockIndex = static_cast<SizeType>((p_to_release - mpData) / GetBlockSize(BlockSizeInBytes));

		  // Check if there is no truncation error
		  KRATOS_DEBUG_CHECK_EQUAL(mFirstAvailableBlockIndex, double(p_to_release - mpData) / GetBlockSize(BlockSizeInBytes));
		  mNumberOfAvailableBlocks++;
		  pPointrerToRelease = nullptr;
	  }

	  void Release() {
		  if (mpData == nullptr)
			  return;
		  delete[] mpData;
		  mpData = nullptr;
		  mNumberOfAvailableBlocks = 0;
		  mFirstAvailableBlockIndex = 0;
	  }

	  static std::size_t Size(std::size_t BlockSizeInBytes, SizeType NumberOfBlocks) {
		  return GetBlockSize(BlockSizeInBytes) * NumberOfBlocks * sizeof(DataType);
	  }

	  std::size_t MemoryUsed(std::size_t BlockSizeInBytes, SizeType NumberOfBlocks) const {
		  if (mpData == nullptr)
			  return sizeof(Chunk);
		  return Size(BlockSizeInBytes, NumberOfBlocks) + sizeof(Chunk);
	  }

	  std::size_t MemoryOverhead(std::size_t BlockSizeInBytes, SizeType NumberOfBlocks) const {
		  if (mpData == nullptr)
			  return sizeof(Chunk);
		  return MemoryUsed(BlockSizeInBytes, NumberOfBlocks) - (BlockSizeInBytes*(NumberOfBlocks - mNumberOfAvailableBlocks));
	  }

	  ///@}
      ///@name Access
      ///@{

	  SizeType GetNumberOfAvailableBlocks() const {
		  return mNumberOfAvailableBlocks;
	  }

      ///@}
      ///@name Inquiry
      ///@{

	  bool HasAvailableBlock() const {
		  return (mNumberOfAvailableBlocks != 0);
	  }

	  bool Has(void* pThePointer, std::size_t BlockSizeInBytes, SizeType NumberOfBlocks) const {
		  return ((pThePointer >= mpData) && (pThePointer <= (mpData + Size(BlockSizeInBytes, NumberOfBlocks) - BlockSizeInBytes)));
	  }

	  bool IsEmpty(SizeType NumberOfBlocks) const {
		  return (mNumberOfAvailableBlocks == NumberOfBlocks);
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
		  rOStream << static_cast<SizeType>(mNumberOfAvailableBlocks) << " blocks are available" << std::endl;
	  }

      ///@}


    private:
      ///@name Member Variables
      ///@{

		DataType * mpData;
		SizeType mFirstAvailableBlockIndex;
		SizeType mNumberOfAvailableBlocks;

      ///@}
      ///@name Operations
      ///@{

		static SizeType DataSize(std::size_t BlockSizeInBytes, SizeType NumberOfBlocks) {
			return GetBlockSize(BlockSizeInBytes) * NumberOfBlocks;
		}

		static std::size_t GetBlockSize(std::size_t BlockSizeInBytes) {
			return (BlockSizeInBytes + sizeof(DataType) - 1) / sizeof(DataType);
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
