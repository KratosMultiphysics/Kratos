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
KRATOS_TEST_CASE_IN_SUITE(GlobalPointerCreateRaw, KratosCoreFastSuite)
{
    int sample_var = 1337;

    auto from_raw = GlobalPointer<int>(&sample_var);

    KRATOS_CHECK_EQUAL(*from_raw, sample_var);
}

KRATOS_TEST_CASE_IN_SUITE(GlobalPointerCreateConstRaw, KratosCoreFastSuite)
{
    const int sample_var = 1337;

    auto from_raw = GlobalPointer<const int>(&sample_var);

    KRATOS_CHECK_EQUAL(*from_raw, sample_var);
}

KRATOS_TEST_CASE_IN_SUITE(GlobalPointerModifyRaw, KratosCoreFastSuite)
{
    int sample_var = 1337;
    int new_value = 42;

    auto from_raw = GlobalPointer<int>(&sample_var);
    *from_raw = new_value;

    KRATOS_CHECK_EQUAL(*from_raw, new_value);
}

// Class
KRATOS_TEST_CASE_IN_SUITE(GlobalPointerCreateClass, KratosCoreFastSuite)
{
    TestClass sample_var(1337);

    auto from_raw = GlobalPointer<TestClass>(&sample_var);

    KRATOS_CHECK_EQUAL(from_raw->getVar(), sample_var.getVar());
    KRATOS_CHECK_EQUAL((*from_raw).getVar(), sample_var.getVar());
}

KRATOS_TEST_CASE_IN_SUITE(GlobalPointerCreateConstClass, KratosCoreFastSuite)
{
    const TestClass sample_var(1337);

	auto from_raw = GlobalPointer<const TestClass>(&sample_var);

    KRATOS_CHECK_EQUAL(from_raw->getVar(), sample_var.getVar());
    KRATOS_CHECK_EQUAL((*from_raw).getVar(), sample_var.getVar());
}

KRATOS_TEST_CASE_IN_SUITE(GlobalPointerModifyClass, KratosCoreFastSuite)
{
    TestClass sample_var(1337);

    auto from_raw = GlobalPointer<TestClass>(&sample_var);

    from_raw->setVar(42);
    sample_var.setVar(42);

    KRATOS_CHECK_EQUAL(from_raw->getVar(), sample_var.getVar());
    KRATOS_CHECK_EQUAL((*from_raw).getVar(), sample_var.getVar());
}

// Kratos::shared_ptr
KRATOS_TEST_CASE_IN_SUITE(GlobalPointerCreateSharedPtr, KratosCoreFastSuite)
{
    typedef Kratos::shared_ptr<TestClass> PtrType;

    auto sample_var = PtrType(new TestClass(1337));
    auto from_shared_ptr = GlobalPointer<TestClass>(sample_var);

    KRATOS_CHECK_EQUAL(from_shared_ptr->getVar(), sample_var->getVar());
    KRATOS_CHECK_EQUAL((*from_shared_ptr).getVar(), sample_var->getVar());
}

KRATOS_TEST_CASE_IN_SUITE(GlobalPointerCreateConstSharedPtr, KratosCoreFastSuite)
{
    typedef Kratos::shared_ptr<TestClass> PtrType;

    const auto sample_var = PtrType(new TestClass(1337));
    auto from_shared_ptr = GlobalPointer<TestClass>(sample_var);

    KRATOS_CHECK_EQUAL(from_shared_ptr->getVar(), sample_var->getVar());
    KRATOS_CHECK_EQUAL((*from_shared_ptr).getVar(), sample_var->getVar());
}

KRATOS_TEST_CASE_IN_SUITE(GlobalPointerModifySharedPtr, KratosCoreFastSuite)
{
    typedef Kratos::shared_ptr<TestClass> PtrType;

    auto sample_var = PtrType(new TestClass(1337));
    auto from_shared_ptr = GlobalPointer<TestClass>(sample_var);

    from_shared_ptr->setVar(42);
    sample_var->setVar(42);

    KRATOS_CHECK_EQUAL(from_shared_ptr->getVar(), sample_var->getVar());
    KRATOS_CHECK_EQUAL((*from_shared_ptr).getVar(), sample_var->getVar());
}

// Kratos::weak_ptr
KRATOS_TEST_CASE_IN_SUITE(GlobalPointerCreateWeakPtr, KratosCoreFastSuite)
{
    typedef Kratos::shared_ptr<TestClass> PtrType;
    typedef Kratos::weak_ptr<TestClass> WeakPtrType;

    auto sample_var = PtrType(new TestClass(1337));
    WeakPtrType weak_var = sample_var;

    auto from_shared_ptr = GlobalPointer<TestClass>(sample_var);

    if(weak_var.lock()) {
        KRATOS_CHECK_EQUAL(from_shared_ptr->getVar(), weak_var.lock()->getVar());
        KRATOS_CHECK_EQUAL((*from_shared_ptr).getVar(), weak_var.lock()->getVar());
    } else {
        KRATOS_CHECK_EQUAL(strcmp("Error", "Unable to lock Kratos::weakptr"), 0);
    }
}

KRATOS_TEST_CASE_IN_SUITE(GlobalPointerCreateConstWeakPtr, KratosCoreFastSuite)
{
    typedef Kratos::shared_ptr<TestClass> PtrType;
    typedef Kratos::weak_ptr<TestClass> WeakPtrType;

    auto sample_var = PtrType(new TestClass(1337));
    const WeakPtrType weak_var = sample_var;

    auto from_shared_ptr = GlobalPointer<TestClass>(weak_var);

    if(weak_var.lock()) {
        KRATOS_CHECK_EQUAL(from_shared_ptr->getVar(), weak_var.lock()->getVar());
        KRATOS_CHECK_EQUAL((*from_shared_ptr).getVar(), weak_var.lock()->getVar());
    } else {
        KRATOS_CHECK_EQUAL(strcmp("Error", "Unable to lock Kratos::weakptr"), 0);
    }
}

KRATOS_TEST_CASE_IN_SUITE(GlobalPointerModifyWeakPtr, KratosCoreFastSuite)
{
    typedef Kratos::shared_ptr<TestClass> PtrType;
    typedef Kratos::weak_ptr<TestClass> WeakPtrType;

    auto sample_var = PtrType(new TestClass(1337));
    WeakPtrType weak_var = sample_var;

    auto from_shared_ptr = GlobalPointer<TestClass>(sample_var);

    if(weak_var.lock()) {
        from_shared_ptr->setVar(42);
        weak_var.lock()->setVar(42);

        KRATOS_CHECK_EQUAL(from_shared_ptr->getVar(), weak_var.lock()->getVar());
        KRATOS_CHECK_EQUAL((*from_shared_ptr).getVar(), weak_var.lock()->getVar());
    } else {
        KRATOS_CHECK_EQUAL(strcmp("Error", "Unable to lock Kratos::weakptr"), 0);
    }
}

} // namespace Testing
} // namespace Kratos
