//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos A. Roig
//

#if !defined(KRATOS_GLOBAL_POINTER_H_INCLUDED )
#define  KRATOS_GLOBAL_POINTER_H_INCLUDED

#include <iostream>

#include "includes/define.h"

#ifdef KRATOS_USING_MPI
#include "mpi.h"
#endif

namespace Kratos {

template<class TDataType>
class GlobalPointer {
private:

  TDataType GetPointer() { return mDataPointer; }

  /// Pointer to the data
  TDataType * mDataPointer;

#ifdef KRATOS_USING_MPI
  /// Rank is enabled only under mpi
  int mRank;
#endif

public:

 /** Default constructor
	* Default constructor
	* This should never be called as we need a local pointer to exists
	*/
	GlobalPointer() = delete;

  /** Constructor by Data
   * Constructor by Data
   * @param Data the data.
   */
  GlobalPointer(TDataType Data) = delete;

  /** Constructor by Data Pointer
   * Constructor by Data Pointer
   * @param DataPointer Pointer to the data.
   */
  GlobalPointer(TDataType * DataPointer)
    : mDataPointer(DataPointer)
#ifdef KRATOS_USING_MPI
    , mRank(GetLocalRank())
#endif
    {
  }

  /** Constructor by boost::shared_ptr
   * Constructor by boost::shared_ptr
   * @param DataPointer Boost Shared Pointer to the Data.
   */
  GlobalPointer(boost::shared_ptr<TDataType> DataPointer)
    : mDataPointer(DataPointer.get())
#ifdef KRATOS_USING_MPI
    , mRank(GetLocalRank())
#endif
    {
  }

  /** Constructor by boost::weak_ptr
   * Constructor by boost::weak_ptr
   * @param DataPointer Boost Weak Pointer to the Data.
   */
  GlobalPointer(boost::weak_ptr<TDataType> DataPointer)
    : mDataPointer(DataPointer.lock().get())
  #ifdef KRATOS_USING_MPI
    , mRank(GetLocalRank())
  #endif
    {
  }

  /** Constructor by boost::unique_ptr
   * Constructor by boost::unique_ptr
   * @param DataPointer Std Unique Pointer to the Data.
   */
  /// Note: Not currently even imported on kratos.
  // GlobalPointer(boost::unique_ptr<TDataType> DataPointer) = delete;

  /** Constructor by std::shared_ptr
   * Constructor by std::shared_ptr
   * @param DataPointer Std Shared Pointer to the Data.
   */
  GlobalPointer(std::shared_ptr<TDataType> DataPointer)
    : mDataPointer(DataPointer.get())
#ifdef KRATOS_USING_MPI
    , mRank(GetLocalRank())
#endif
    {
  }

  /** Constructor by std::weak_ptr
   * Constructor by std::weak_ptr
   * @param DataPointer Std Weak Pointer to the Data.
   */
  GlobalPointer(std::weak_ptr<TDataType> DataPointer)
    : mDataPointer(DataPointer.lock().get())
  #ifdef KRATOS_USING_MPI
    , mRank(GetLocalRank())
  #endif
      {
    }

  /** Constructor by std::unique_ptr
   * Constructor by std::unique_ptr
   * @param DataPointer Std Unique Pointer to the Data.
   */
  GlobalPointer(std::unique_ptr<TDataType> DataPointer) = delete;

  /** Copy constructor
   * Copy constructor
   * @ rOther: Copied GlobalPointer
   */
  GlobalPointer(const GlobalPointer & rOther)
    : mDataPointer(rOther.mDataPointer)
#ifdef KRATOS_USING_MPI
    , mRank(rOther.mRank)
#endif
    {
  }


  /** Default Destructor
   * Default Destructor.
   */
  ~GlobalPointer() {
  }

  /** Pointer Operator
  * Pointer Operator
  */
  TDataType & operator*() {
	  return *mDataPointer;
  }

  /** Pointer Operator
  * Pointer Operator
  */
  TDataType const& operator*() const {
	  return *mDataPointer;
  }

  /** Arrow Operator
   * Arrow Operator
   */
  TDataType * operator->() {
    return mDataPointer;
  }

  /** Assignment Operator
   * Assignment Operator
   */
  GlobalPointer & operator=(const GlobalPointer & rOther) {
    mDataPointer = rOther.mDataPointer;
#ifdef KRATOS_USING_MPI
    mRank = rOther.mRank;
#endif
    return *this;
  }

  /** Returns the rank of the data owner
   * Returns the rank of the data owner
   * @return rank of the data owner
   */
  int GetRank() const {
#ifdef KRATOS_USING_MPI
    return mRank;
#else
    return 0;
#endif
  }

  /** Fills buffer with the GlobalPoiter data
   * Fills buffer with the GlobalPoiter data
   * @param buffer Object data buffer
   */
  void Save(char * buffer) {
    memcpy(buffer, this, sizeof(GlobalPointer));
  }

  /** Restores the GlobalPoiter with the data from the buffer
   * Restores the GlobalPoiter with the data from the buffer
   * @param buffer Object data buffer
   */
  void Load(char * buffer) {
    memcpy(this, buffer, sizeof(GlobalPointer));
  }

private:
  /** Returns the rank of the current process
   * Returns the rank of the current process
   * @return Rank of the current process
   */
  int GetLocalRank() const {
#ifdef KRATOS_USING_MPI
    int mpi_rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    return mpi_rank;
#else
    return 0;
#endif
  }

};

} // namespace Kratos

#endif // KRATOS_GLOBAL_POINTER_H_INCLUDED
