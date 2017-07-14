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

namespace Kratos {

template<class TDataType>
class GlobalPointer {
private:

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
  GlobalPointer(TDataType * DataPointer, int Rank = 0)
    : mDataPointer(DataPointer)
#ifdef KRATOS_USING_MPI
    , mRank(Rank)
#endif
    {
  }

  /** Constructor by boost::shared_ptr
   * Constructor by boost::shared_ptr
   * @param DataPointer Boost Shared Pointer to the Data.
   */
  GlobalPointer(boost::shared_ptr<TDataType> DataPointer, int Rank = 0)
    : mDataPointer(DataPointer.get())
#ifdef KRATOS_USING_MPI
    , mRank(Rank)
#endif
    {
  }

  /** Constructor by boost::weak_ptr
   * Constructor by boost::weak_ptr
   * @param DataPointer Boost Weak Pointer to the Data.
   */
  GlobalPointer(boost::weak_ptr<TDataType> DataPointer, int Rank = 0)
    : mDataPointer(DataPointer.lock().get())
  #ifdef KRATOS_USING_MPI
    , mRank(Rank)
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
  GlobalPointer(std::shared_ptr<TDataType> DataPointer, int Rank = 0)
    : mDataPointer(DataPointer.get())
#ifdef KRATOS_USING_MPI
    , mRank(Rank)
#endif
    {
  }

  /** Constructor by std::weak_ptr
   * Constructor by std::weak_ptr
   * @param DataPointer Std Weak Pointer to the Data.
   */
  GlobalPointer(std::weak_ptr<TDataType> DataPointer, int Rank = 0)
    : mDataPointer(DataPointer.lock().get())
  #ifdef KRATOS_USING_MPI
    , mRank(Rank)
  #endif
      {
    }

  /** Constructor by std::unique_ptr
   * Constructor by std::unique_ptr
   * @param DataPointer Std Unique Pointer to the Data.
   */
  GlobalPointer(std::unique_ptr<TDataType> DataPointer, int Rank = 0) = delete;

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

  /** Move constructor
   * Move constructor
   * @ rOther: Copied GlobalPointer
   */
  GlobalPointer(const GlobalPointer && rOther)
    : mDataPointer(std::move(rOther.mDataPointer))
#ifdef KRATOS_USING_MPI
    , mRank(std::move(rOther.mRank))
#endif
    {
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

  /** Const Pointer Operator
  * Const Pointer Operator
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

  /** Const Arrow Operator
   * Const Arrow Operator
   */
  TDataType const* operator->() const {
    return mDataPointer;
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

  /** Returns the rank of the global pointer
   * Returns the rank of the global pointer data is located or 0 if no mpi
   * @return rank of the global pointer data or 0
   */
  int GetRank() {
#ifdef KRATOS_USING_MPI
    return this->rank;
#else
    return 0;
#endif
  }
};


} // namespace Kratos

#endif // KRATOS_GLOBAL_POINTER_H_INCLUDED
