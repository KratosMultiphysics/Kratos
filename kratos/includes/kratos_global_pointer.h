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

template<class T>
class GlobalPointer {
private:

  /** Default constructor
   * Default constructor
   * This should never be called as we need a local pointer to exists
   */
  GlobalPointer() {}

  T GetPointer() { return mBasePointer; }

  /// Local pointer and its rank. This must never be seen outside this class.
  /// Do not make it protected if you derive from this class.
  T mBasePointer;
  int mRank;

public:

  /** Constructor by a local BasePointer
   * Constructor by a local pointer
   * @param BasePointer BasePointer of the local variable.
   */
  GlobalPointer(const T & BasePointer) {
    mBasePointer = BasePointer;
    mRank = GetLocalRank();
  }

  /** Default Destructor
   * Default Destructor.
   * The descturctor DOES NOT free the local pointer.
   */
  ~GlobalPointer() {
    mBasePointer = nullptr;
    mRank = -1;
  }

  /**
   * Pointer Operator
   */
  auto operator*() -> decltype(*mBasePointer) {
    if(GetLocalRank() == mRank) {
      return *mBasePointer;
    }

    KRATOS_ERROR << "Trying to access invalid memory space" << std::endl;
  }

  /**
   * Arrow Operator
   */
  auto operator->() -> decltype(mBasePointer) {
    if(GetLocalRank() == mRank) {
      return mBasePointer;
    }

    KRATOS_ERROR << "Trying to access invalid memory space" << std::endl;
  }

  /** Returns the rank of the BasePointer owner
   * Returns the rank of the BasePointer owner
   * @return Rank of the BasePointer owner
   */
  int GetRank() {
    return mRank;
  }

  /** Returns the number of ranks
   * Returns the number of ranks
   * @return Number of ranks
   */
  int GetSize() {
#ifdef KRATOS_USING_MPI
    int mpi_size;

    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    return mpi_size;
#else
    return 0;
#endif
  }

  /** Returns the rank of the current process
   * Returns the rank of the current process
   * @return Rank of the current process
   */
  int GetLocalRank() {
#ifdef KRATOS_USING_MPI
    int mpi_rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    return mpi_size;
#else
    return 0;
#endif
  }


};

} // namespace Kratos

#endif // KRATOS_GLOBAL_POINTER_H_INCLUDED
