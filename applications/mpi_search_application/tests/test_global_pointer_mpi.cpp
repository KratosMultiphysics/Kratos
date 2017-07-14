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

#include "boost/smart_ptr.hpp"

#include "testing/testing.h"
#include "includes/global_pointer.h"

#include "mpi.h"

namespace Kratos {
namespace Testing {

/// Parallel tests

KRATOS_TEST_CASE_IN_SUITE(GlobalPointerGatherRaw, KratosCoreFastSuit)
{
  int mpi_size;
  int mpi_rank;

  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  int sample_var = 1337 + mpi_rank;

	auto from_raw_origin = GlobalPointer<int>(&sample_var, mpi_rank);

  std::size_t gp_size = sizeof(GlobalPointer<int>);

  char * raw_gather_send_buffer = new char[gp_size];
  char * raw_gather_recv_buffer = new char[gp_size * mpi_size];

  KRATOS_CHECK_NOT_EQUAL(raw_gather_send_buffer, nullptr);
  KRATOS_CHECK_NOT_EQUAL(raw_gather_recv_buffer, nullptr);

  from_raw_origin.Save(raw_gather_send_buffer);

  MPI_Allgather(
    raw_gather_send_buffer, gp_size, MPI_CHAR,
    raw_gather_recv_buffer, gp_size, MPI_CHAR,
    MPI_COMM_WORLD
  );

  auto from_raw_remote = GlobalPointer<int>(nullptr);

  for(int i = 0; i < mpi_size; i += 1) {
    from_raw_remote.Load(&raw_gather_recv_buffer[i * gp_size]);
    KRATOS_CHECK_EQUAL(from_raw_remote.GetRank(), i);
    if(mpi_rank == i) {
      KRATOS_CHECK_EQUAL(*from_raw_remote, sample_var);
    }
  }

  delete[] raw_gather_send_buffer;
  delete[] raw_gather_recv_buffer;
}

} // namespace Testing
} // namespace Kratos
