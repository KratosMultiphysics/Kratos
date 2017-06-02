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
    TestClass(int magicNumber) { this->mMagicNumber = magicNumber; }
    TestClass(const TestClass & rOther) { this->mMagicNumber = rOther.mMagicNumber; }

    ~TestClass() {}

    int getVar() const { return this->mMagicNumber; }

    void setVar(int magicNumber) { this->mMagicNumber = magicNumber; }

    int mMagicNumber;
};

KRATOS_TEST_CASE_IN_SUITE(GlobalPointerCreateRaw, KratosCoreFastSuit)
{
  int sampleVar = 1337;

	auto fromRaw = GlobalPointer<int>(&sampleVar);

  KRATOS_CHECK_EQUAL(*fromRaw, sampleVar);
}

KRATOS_TEST_CASE_IN_SUITE(GlobalPointerCreateConstRaw, KratosCoreFastSuit)
{
  const int sampleVar = 1337;

	auto fromRaw = GlobalPointer<const int>(&sampleVar);

  KRATOS_CHECK_EQUAL(*fromRaw, sampleVar);
}

KRATOS_TEST_CASE_IN_SUITE(GlobalPointerModifyRaw, KratosCoreFastSuit)
{
  int sampleVar = 1337;
  int newVal = 42;

	auto fromRaw = GlobalPointer<int>(&sampleVar);
  *fromRaw = newVal;

  KRATOS_CHECK_EQUAL(*fromRaw, newVal);
}

KRATOS_TEST_CASE_IN_SUITE(GlobalPointerCreateClass, KratosCoreFastSuit)
{
  TestClass sampleVar(1337);

	auto fromRaw = GlobalPointer<TestClass>(&sampleVar);

  KRATOS_CHECK_EQUAL(fromRaw->getVar(), sampleVar.getVar());
  KRATOS_CHECK_EQUAL((*fromRaw).getVar(), sampleVar.getVar());
}

KRATOS_TEST_CASE_IN_SUITE(GlobalPointerCreateConstClass, KratosCoreFastSuit)
{
  const TestClass sampleVar(1337);

	auto fromRaw = GlobalPointer<const TestClass>(&sampleVar);

  KRATOS_CHECK_EQUAL(fromRaw->getVar(), sampleVar.getVar());
  KRATOS_CHECK_EQUAL((*fromRaw).getVar(), sampleVar.getVar());
}

KRATOS_TEST_CASE_IN_SUITE(GlobalPointerModifyClass, KratosCoreFastSuit)
{
  TestClass sampleVar(1337);

	auto fromRaw = GlobalPointer<TestClass>(&sampleVar);

  fromRaw->setVar(42);
  sampleVar.setVar(42);

  KRATOS_CHECK_EQUAL(fromRaw->getVar(), sampleVar.getVar());
  KRATOS_CHECK_EQUAL((*fromRaw).getVar(), sampleVar.getVar());
}

// Test global_ptr<shared_ptr>
KRATOS_TEST_CASE_IN_SUITE(GlobalPointerCreateBoostSharedPtr, KratosCoreFastSuit)
{
  typedef boost::shared_ptr<TestClass> BoostPtrType;

  auto sampleVar = BoostPtrType(new TestClass(1337));
	auto fromBoost = GlobalPointer<TestClass>(sampleVar);

  KRATOS_CHECK_EQUAL(fromBoost->getVar(), sampleVar->getVar());
  KRATOS_CHECK_EQUAL((*fromBoost).getVar(), sampleVar->getVar());
}

KRATOS_TEST_CASE_IN_SUITE(GlobalPointerModifyBoostSharedPtr, KratosCoreFastSuit)
{
  typedef boost::shared_ptr<TestClass> BoostPtrType;

  auto sampleVar = BoostPtrType(new TestClass(1337));
	auto fromBoost = GlobalPointer<TestClass>(sampleVar);

  fromBoost->setVar(42);
  sampleVar->setVar(42);

  KRATOS_CHECK_EQUAL(fromBoost->getVar(), sampleVar->getVar());
  KRATOS_CHECK_EQUAL((*fromBoost).getVar(), sampleVar->getVar());
}

/// Parallel tests

#ifdef KRATOS_USING_MPI
KRATOS_TEST_CASE_IN_SUITE(GlobalPointerGatherRaw, KratosCoreFastSuit)
{
  int mpi_size;
  int mpi_rank;

  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  int sampleVar = 1337 + mpi_rank;

	auto fromRawOrigin = GlobalPointer<int>(&sampleVar);

  std::size_t gpSize = sizeof(GlobalPointer<int>);

  char * rawGather = (char *)malloc(gpSize * mpi_size);
  KRATOS_CHECK_NOT_EQUAL(rawGather, nullptr);

  MPI_Allgather(
    fromRawOrigin.ToRaw(), gpSize, MPI_CHAR,
    rawGather,             gpSize, MPI_CHAR,
    MPI_COMM_WORLD
  );

  auto fromRawRemote = GlobalPointer<int>(nullptr);

  for(int i = 0; i < mpi_size; i += 1) {
    fromRawRemote.FromRaw(&rawGather[i * gpSize]);
    KRATOS_CHECK_EQUAL(fromRawRemote.GetRank(), i);
    if(mpi_rank == i) {
      KRATOS_CHECK_EQUAL(*fromRawRemote, sampleVar);
    }
  }

  free(rawGather);
}
#endif

// Test shared_ptr<global_ptr>

} // namespace Testing
} // namespace Kratos
