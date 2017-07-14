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

namespace Kratos {
namespace Testing {

class TestClass {
  private:
    TestClass() {}

  public:
    TestClass(int MagicNbr) { this->mMagicNbr = MagicNbr; }
    TestClass(const TestClass & rOther) { this->mMagicNbr = rOther.mMagicNbr; }

    ~TestClass() {}

    int getVar() const { return this->mMagicNbr; }

    void setVar(int MagicNbr) { this->mMagicNbr = MagicNbr; }

    int mMagicNbr;
};

// Basic Type
KRATOS_TEST_CASE_IN_SUITE(GlobalPointerCreateRaw, KratosCoreFastSuit)
{
  int sample_var = 1337;

	auto from_raw = GlobalPointer<int>(&sample_var);

  KRATOS_CHECK_EQUAL(*from_raw, sample_var);
}

KRATOS_TEST_CASE_IN_SUITE(GlobalPointerCreateConstRaw, KratosCoreFastSuit)
{
  const int sample_var = 1337;

	auto from_raw = GlobalPointer<const int>(&sample_var);

  KRATOS_CHECK_EQUAL(*from_raw, sample_var);
}

KRATOS_TEST_CASE_IN_SUITE(GlobalPointerModifyRaw, KratosCoreFastSuit)
{
  int sample_var = 1337;
  int new_value = 42;

	auto from_raw = GlobalPointer<int>(&sample_var);
  *from_raw = new_value;

  KRATOS_CHECK_EQUAL(*from_raw, new_value);
}

// Class
KRATOS_TEST_CASE_IN_SUITE(GlobalPointerCreateClass, KratosCoreFastSuit)
{
  TestClass sample_var(1337);

	auto from_raw = GlobalPointer<TestClass>(&sample_var);

  KRATOS_CHECK_EQUAL(from_raw->getVar(), sample_var.getVar());
  KRATOS_CHECK_EQUAL((*from_raw).getVar(), sample_var.getVar());
}

KRATOS_TEST_CASE_IN_SUITE(GlobalPointerCreateConstClass, KratosCoreFastSuit)
{
  const TestClass sample_var(1337);

	auto from_raw = GlobalPointer<const TestClass>(&sample_var);

  KRATOS_CHECK_EQUAL(from_raw->getVar(), sample_var.getVar());
  KRATOS_CHECK_EQUAL((*from_raw).getVar(), sample_var.getVar());
}

KRATOS_TEST_CASE_IN_SUITE(GlobalPointerModifyClass, KratosCoreFastSuit)
{
  TestClass sample_var(1337);

	auto from_raw = GlobalPointer<TestClass>(&sample_var);

  from_raw->setVar(42);
  sample_var.setVar(42);

  KRATOS_CHECK_EQUAL(from_raw->getVar(), sample_var.getVar());
  KRATOS_CHECK_EQUAL((*from_raw).getVar(), sample_var.getVar());
}

// Boost::shared_ptr
KRATOS_TEST_CASE_IN_SUITE(GlobalPointerCreateBoostSharedPtr, KratosCoreFastSuit)
{
  typedef boost::shared_ptr<TestClass> BoostPtrType;

  auto sample_var = BoostPtrType(new TestClass(1337));
	auto from_boost = GlobalPointer<TestClass>(sample_var);

  KRATOS_CHECK_EQUAL(from_boost->getVar(), sample_var->getVar());
  KRATOS_CHECK_EQUAL((*from_boost).getVar(), sample_var->getVar());
}

KRATOS_TEST_CASE_IN_SUITE(GlobalPointerCreateConstBoostSharedPtr, KratosCoreFastSuit)
{
  typedef boost::shared_ptr<TestClass> BoostPtrType;

  const auto sample_var = BoostPtrType(new TestClass(1337));
	auto from_boost = GlobalPointer<TestClass>(sample_var);

  KRATOS_CHECK_EQUAL(from_boost->getVar(), sample_var->getVar());
  KRATOS_CHECK_EQUAL((*from_boost).getVar(), sample_var->getVar());
}

KRATOS_TEST_CASE_IN_SUITE(GlobalPointerModifyBoostSharedPtr, KratosCoreFastSuit)
{
  typedef boost::shared_ptr<TestClass> BoostPtrType;

  auto sample_var = BoostPtrType(new TestClass(1337));
	auto from_boost = GlobalPointer<TestClass>(sample_var);

  from_boost->setVar(42);
  sample_var->setVar(42);

  KRATOS_CHECK_EQUAL(from_boost->getVar(), sample_var->getVar());
  KRATOS_CHECK_EQUAL((*from_boost).getVar(), sample_var->getVar());
}

// Boost::weak_ptr
KRATOS_TEST_CASE_IN_SUITE(GlobalPointerCreateBoostWeakPtr, KratosCoreFastSuit)
{
  typedef boost::shared_ptr<TestClass> BoostPtrType;
  typedef boost::weak_ptr<TestClass> BoostWeakPtrType;

  auto sample_var = BoostPtrType(new TestClass(1337));
  BoostWeakPtrType weak_var = sample_var;

	auto from_boost = GlobalPointer<TestClass>(sample_var);

  if(weak_var.lock()) {
    KRATOS_CHECK_EQUAL(from_boost->getVar(), weak_var.lock()->getVar());
    KRATOS_CHECK_EQUAL((*from_boost).getVar(), weak_var.lock()->getVar());
  } else {
    KRATOS_CHECK_EQUAL("Error", "Unable to lock boost::weakptr");
  }
}

KRATOS_TEST_CASE_IN_SUITE(GlobalPointerCreateConstBoostWeakPtr, KratosCoreFastSuit)
{
  typedef boost::shared_ptr<TestClass> BoostPtrType;
  typedef boost::weak_ptr<TestClass> BoostWeakPtrType;

  auto sample_var = BoostPtrType(new TestClass(1337));
  const BoostWeakPtrType weak_var = sample_var;

	auto from_boost = GlobalPointer<TestClass>(weak_var);

  if(weak_var.lock()) {
    KRATOS_CHECK_EQUAL(from_boost->getVar(), weak_var.lock()->getVar());
    KRATOS_CHECK_EQUAL((*from_boost).getVar(), weak_var.lock()->getVar());
  } else {
    KRATOS_CHECK_EQUAL("Error", "Unable to lock boost::weakptr");
  }
}

KRATOS_TEST_CASE_IN_SUITE(GlobalPointerModifyBoostWeakPtr, KratosCoreFastSuit)
{
  typedef boost::shared_ptr<TestClass> BoostPtrType;
  typedef boost::weak_ptr<TestClass> BoostWeakPtrType;

  auto sample_var = BoostPtrType(new TestClass(1337));
  BoostWeakPtrType weak_var = sample_var;

	auto from_boost = GlobalPointer<TestClass>(sample_var);

  if(weak_var.lock()) {
    from_boost->setVar(42);
    weak_var.lock()->setVar(42);

    KRATOS_CHECK_EQUAL(from_boost->getVar(), weak_var.lock()->getVar());
    KRATOS_CHECK_EQUAL((*from_boost).getVar(), weak_var.lock()->getVar());
  } else {
    KRATOS_CHECK_EQUAL("Error", "Unable to lock boost::weakptr");
  }
}

/// Parallel tests

#ifdef KRATOS_USING_MPI
KRATOS_TEST_CASE_IN_SUITE(GlobalPointerGatherRaw, KratosCoreFastSuit)
{
  int mpi_size;
  int mpi_rank;

  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  int sample_var = 1337 + mpi_rank;

	auto from_raw_origin = GlobalPointer<int>(&sample_var);

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
#endif

// Test shared_ptr<global_ptr>

} // namespace Testing
} // namespace Kratos
