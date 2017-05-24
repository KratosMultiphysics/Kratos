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
#include "includes/kratos_global_pointer.h"

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
  int * sampleRawPointer = &sampleVar;

	auto fromRaw = GlobalPointer<int*>(sampleRawPointer);

  KRATOS_CHECK_EQUAL(*fromRaw, sampleVar);
}

KRATOS_TEST_CASE_IN_SUITE(GlobalPointerModifyRaw, KratosCoreFastSuit)
{
  int sampleVar = 1337;
  int newVal = 42;

  int * sampleRawPointer = &sampleVar;

	auto fromRaw = GlobalPointer<int*>(sampleRawPointer);
  *fromRaw = newVal;

  KRATOS_CHECK_EQUAL(*fromRaw, newVal);
}

KRATOS_TEST_CASE_IN_SUITE(GlobalPointerCreateClass, KratosCoreFastSuit)
{
  TestClass sampleVar(1337);
  TestClass * sampleRawPointer = &sampleVar;

	auto fromRaw = GlobalPointer<TestClass*>(sampleRawPointer);

  KRATOS_CHECK_EQUAL(fromRaw->getVar(), sampleVar.getVar());
  KRATOS_CHECK_EQUAL((*fromRaw).getVar(), sampleVar.getVar());
}

KRATOS_TEST_CASE_IN_SUITE(GlobalPointerModifyClass, KratosCoreFastSuit)
{
  TestClass sampleVar(1337);
  TestClass * sampleRawPointer = &sampleVar;

	auto fromRaw = GlobalPointer<TestClass*>(sampleRawPointer);

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
	auto fromBoost = GlobalPointer<BoostPtrType>(sampleVar);

  KRATOS_CHECK_EQUAL(fromBoost->getVar(), sampleVar->getVar());
  KRATOS_CHECK_EQUAL((*fromBoost).getVar(), sampleVar->getVar());
}

KRATOS_TEST_CASE_IN_SUITE(GlobalPointerModifyBoostSharedPtr, KratosCoreFastSuit)
{
  typedef boost::shared_ptr<TestClass> BoostPtrType;

  auto sampleVar = BoostPtrType(new TestClass(1337));
	auto fromBoost = GlobalPointer<BoostPtrType>(sampleVar);

  fromBoost->setVar(42);
  sampleVar->setVar(42);

  KRATOS_CHECK_EQUAL(fromBoost->getVar(), sampleVar->getVar());
  KRATOS_CHECK_EQUAL((*fromBoost).getVar(), sampleVar->getVar());
}

// Make sure if no one references the pointer, the global pointer does NOT
// interfeer with the share_ptr de-allocation mechanism
KRATOS_TEST_CASE_IN_SUITE(GlobalPointerExpiredPointer, KratosCoreFastSuit)
{
  typedef boost::shared_ptr<TestClass> BoostPtrType;

  BoostPtrType * sampleVar = new BoostPtrType(new TestClass(1337));
	auto fromBoost = GlobalPointer<BoostPtrType>(*sampleVar);

  auto realPtr = &*fromBoost;

  KRATOS_CHECK_EQUAL(fromBoost->getVar(), (*sampleVar)->getVar());

  delete sampleVar;
  // KRATOS_CHECK_EXCEPTION_RAISED((*sampleVar)->getVar(), Exception);
  KRATOS_CHECK_EXCEPTION_RAISED(fromBoost->getVar(), Exception);

  // This is not a tests but a warning
  // Beware that, while possible, this line is incorrect:
  KRATOS_CHECK_EQUAL(realPtr->getVar(), (*sampleVar)->getVar());
}

// Test shared_ptr<global_ptr>

} // namespace Testing
} // namespace Kratos
