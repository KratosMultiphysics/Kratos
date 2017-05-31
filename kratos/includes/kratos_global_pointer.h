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

  TDataType GetPointer() { return mBaseDataPtr; }

  /// Local pointer and its rank. This must never be seen outside this class.
  /// Do not make it protected if you derive from this class.
  TDataType * mBaseDataPtr;
  int mRank;

public:

	/** Default constructor
	* Default constructor
	* This should never be called as we need a local pointer to exists
	*/
	GlobalPointer() = delete;

  /** Constructor by a local BaseData
   * Constructor by a local pointer
   * @param BaseData BaseData of the local variable.
   */
  GlobalPointer(TDataType * BaseDataPtr)
    : mBaseDataPtr(BaseDataPtr)
    , mRank(GetLocalRank()) {
  }

  /** Constructor by boost::shared_ptr
   * Constructor by boost::shared_ptr
   * @param BaseData BaseData of the local variable.
   */
  GlobalPointer(boost::shared_ptr<TDataType> BaseSharedPtr)
    : mBaseDataPtr(&*BaseSharedPtr)
    , mRank(GetLocalRank()) {
  }

  /** Default Destructor
   * Default Destructor.
   */
  ~GlobalPointer() {
  }

  /**
   * Pointer Operator
   */
  TDataType & operator*() {
    return *mBaseDataPtr;
  }

  /**
   * Arrow Operator
   */
  TDataType * operator->() {
    return mBaseDataPtr;
  }

  /** Returns the rank of the BasePointer owner
   * Returns the rank of the BasePointer owner
   * @return Rank of the BasePointer owner
   */
  int GetRank() const {
    return mRank;
  }

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
