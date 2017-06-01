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
#include "tests/geometries/test_geometry.h"
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

    int getVar() { return this->mMagicNumber; }
    void setVar(int magicNumber) { this->mMagicNumber = magicNumber; }

    int mMagicNumber;
};

KRATOS_TEST_CASE_IN_SUITE(GlobalPointerCreateRaw, KratosCoreFastSuit)
{
  int sampleVar = 1337;

	auto fromRaw = GlobalPointer<int>(&sampleVar);

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

// Not being tested. Only relevant in DEBUG
/**
 * This checks the access of invalid ranks while running in serial.
 * DO NOT USE this in a production code.
 */
// KRATOS_TEST_CASE_IN_SUITE(GlobalPointerInvalidRankRaw, KratosCoreFastSuit)
// {
//   int sampleVar = 1337;
//   int newRank = 2;
//   int * sampleRawPointer = &sampleVar;
//
// 	auto fromRaw = GlobalPointer<int>(sampleRawPointer);
//
//   // Inject a new value in the rank private member
//   int * rankIndirection = (int *)((long)(&fromRaw) + sizeof(void *) * 1);
//   *rankIndirection = newRank;
//
//   KRATOS_CHECK_EQUAL(fromRaw.GetRank(), newRank);
//   KRATOS_CHECK_EXCEPTION_RAISED(*fromRaw, Exception);
// }

KRATOS_TEST_CASE_IN_SUITE(GlobalPointerCreateClass, KratosCoreFastSuit)
{
  TestClass sampleVar(1337);

	auto fromRaw = GlobalPointer<TestClass>(&sampleVar);

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

#ifdef KRATOS_USING_MPI
// Parallel tests:
KRATOS_TEST_CASE_IN_SUITE(GlobalPointerGatherRaw, KratosCoreFastSuit)
{
  int sampleVar = 1337;

	auto fromRaw = GlobalPointer<int>(&sampleVar);

  KRATOS_CHECK_EQUAL(*fromRaw, sampleVar);
}
#endif

// Test shared_ptr<global_ptr>

} // namespace Testing
} // namespace Kratos
