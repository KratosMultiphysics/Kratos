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
    explicit TestClass(int MagicNbr) { this->mMagicNbr = MagicNbr; }
    TestClass(const TestClass & rOther) { this->mMagicNbr = rOther.mMagicNbr; }

    ~TestClass() {}

    int getVar() const { return this->mMagicNbr; }

    void setVar(int MagicNbr) { this->mMagicNbr = MagicNbr; }

    int mMagicNbr;
};

// Basic Type
TEST(GlobalPointerCreateRaw, KratosCoreFastSuite)
{
    int sample_var = 1337;

    auto from_raw = GlobalPointer<int>(&sample_var);

    KRATOS_EXPECT_EQ(*from_raw, sample_var);
}

TEST(GlobalPointerCreateConstRaw, KratosCoreFastSuite)
{
    const int sample_var = 1337;

    auto from_raw = GlobalPointer<const int>(&sample_var);

    KRATOS_EXPECT_EQ(*from_raw, sample_var);
}

TEST(GlobalPointerModifyRaw, KratosCoreFastSuite)
{
    int sample_var = 1337;
    int new_value = 42;

    auto from_raw = GlobalPointer<int>(&sample_var);
    *from_raw = new_value;

    KRATOS_EXPECT_EQ(*from_raw, new_value);
}

// Class
TEST(GlobalPointerCreateClass, KratosCoreFastSuite)
{
    TestClass sample_var(1337);

    auto from_raw = GlobalPointer<TestClass>(&sample_var);

    KRATOS_EXPECT_EQ(from_raw->getVar(), sample_var.getVar());
    KRATOS_EXPECT_EQ((*from_raw).getVar(), sample_var.getVar());
}

TEST(GlobalPointerCreateConstClass, KratosCoreFastSuite)
{
    const TestClass sample_var(1337);

	auto from_raw = GlobalPointer<const TestClass>(&sample_var);

    KRATOS_EXPECT_EQ(from_raw->getVar(), sample_var.getVar());
    KRATOS_EXPECT_EQ((*from_raw).getVar(), sample_var.getVar());
}

TEST(GlobalPointerModifyClass, KratosCoreFastSuite)
{
    TestClass sample_var(1337);

    auto from_raw = GlobalPointer<TestClass>(&sample_var);

    from_raw->setVar(42);
    sample_var.setVar(42);

    KRATOS_EXPECT_EQ(from_raw->getVar(), sample_var.getVar());
    KRATOS_EXPECT_EQ((*from_raw).getVar(), sample_var.getVar());
}

// Kratos::shared_ptr
TEST(GlobalPointerCreateSharedPtr, KratosCoreFastSuite)
{
    typedef Kratos::shared_ptr<TestClass> PtrType;

    auto sample_var = PtrType(new TestClass(1337));
    auto from_shared_ptr = GlobalPointer<TestClass>(sample_var);

    KRATOS_EXPECT_EQ(from_shared_ptr->getVar(), sample_var->getVar());
    KRATOS_EXPECT_EQ((*from_shared_ptr).getVar(), sample_var->getVar());
}

TEST(GlobalPointerCreateConstSharedPtr, KratosCoreFastSuite)
{
    typedef Kratos::shared_ptr<TestClass> PtrType;

    const auto sample_var = PtrType(new TestClass(1337));
    auto from_shared_ptr = GlobalPointer<TestClass>(sample_var);

    KRATOS_EXPECT_EQ(from_shared_ptr->getVar(), sample_var->getVar());
    KRATOS_EXPECT_EQ((*from_shared_ptr).getVar(), sample_var->getVar());
}

TEST(GlobalPointerModifySharedPtr, KratosCoreFastSuite)
{
    typedef Kratos::shared_ptr<TestClass> PtrType;

    auto sample_var = PtrType(new TestClass(1337));
    auto from_shared_ptr = GlobalPointer<TestClass>(sample_var);

    from_shared_ptr->setVar(42);
    sample_var->setVar(42);

    KRATOS_EXPECT_EQ(from_shared_ptr->getVar(), sample_var->getVar());
    KRATOS_EXPECT_EQ((*from_shared_ptr).getVar(), sample_var->getVar());
}

// Kratos::weak_ptr
TEST(GlobalPointerCreateWeakPtr, KratosCoreFastSuite)
{
    typedef Kratos::shared_ptr<TestClass> PtrType;
    typedef Kratos::weak_ptr<TestClass> WeakPtrType;

    auto sample_var = PtrType(new TestClass(1337));
    WeakPtrType weak_var = sample_var;

    auto from_shared_ptr = GlobalPointer<TestClass>(sample_var);

    if(weak_var.lock()) {
        KRATOS_EXPECT_EQ(from_shared_ptr->getVar(), weak_var.lock()->getVar());
        KRATOS_EXPECT_EQ((*from_shared_ptr).getVar(), weak_var.lock()->getVar());
    } else {
        KRATOS_EXPECT_EQ(strcmp("Error", "Unable to lock Kratos::weakptr"), 0);
    }
}

TEST(GlobalPointerCreateConstWeakPtr, KratosCoreFastSuite)
{
    typedef Kratos::shared_ptr<TestClass> PtrType;
    typedef Kratos::weak_ptr<TestClass> WeakPtrType;

    auto sample_var = PtrType(new TestClass(1337));
    const WeakPtrType weak_var = sample_var;

    auto from_shared_ptr = GlobalPointer<TestClass>(weak_var);

    if(weak_var.lock()) {
        KRATOS_EXPECT_EQ(from_shared_ptr->getVar(), weak_var.lock()->getVar());
        KRATOS_EXPECT_EQ((*from_shared_ptr).getVar(), weak_var.lock()->getVar());
    } else {
        KRATOS_EXPECT_EQ(strcmp("Error", "Unable to lock Kratos::weakptr"), 0);
    }
}

TEST(GlobalPointerModifyWeakPtr, KratosCoreFastSuite)
{
    typedef Kratos::shared_ptr<TestClass> PtrType;
    typedef Kratos::weak_ptr<TestClass> WeakPtrType;

    auto sample_var = PtrType(new TestClass(1337));
    WeakPtrType weak_var = sample_var;

    auto from_shared_ptr = GlobalPointer<TestClass>(sample_var);

    if(weak_var.lock()) {
        from_shared_ptr->setVar(42);
        weak_var.lock()->setVar(42);

        KRATOS_EXPECT_EQ(from_shared_ptr->getVar(), weak_var.lock()->getVar());
        KRATOS_EXPECT_EQ((*from_shared_ptr).getVar(), weak_var.lock()->getVar());
    } else {
        KRATOS_EXPECT_EQ(strcmp("Error", "Unable to lock Kratos::weakptr"), 0);
    }
}

TEST(GlobalPointerCompare, KratosCoreFastSuite)
{
    std::array<int,3> values={21, 1, 8};

    auto global_ptr_1 = GlobalPointer<int>(&values[0]);
    auto global_ptr_2 = GlobalPointer<int>(&values[1]);
    GlobalPointerCompare<int> compare;
    KRATOS_EXPECT_TRUE(compare(global_ptr_1, global_ptr_2));
    KRATOS_EXPECT_FALSE(compare(global_ptr_2, global_ptr_1));
    KRATOS_EXPECT_FALSE(compare(global_ptr_1, global_ptr_1));
    KRATOS_EXPECT_FALSE(compare(global_ptr_2, global_ptr_2));
}


} // namespace Testing
} // namespace Kratos
