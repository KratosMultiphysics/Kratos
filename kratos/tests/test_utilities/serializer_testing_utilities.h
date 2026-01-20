#include "includes/serializer.h"

namespace Kratos::Testing
{

class KRATOS_API(KRATOS_CORE) TestClass 
{
public:
    int mValue;
    TestClass() : mValue(0) {}
    TestClass(int v) : mValue(v) {}
    virtual ~TestClass() = default;
private:
    friend class Kratos::Serializer;
    void save(Serializer& rSerializer) const {
        rSerializer.save("Value", mValue);
    }
    void load(Serializer& rSerializer) {
        rSerializer.load("Value", mValue);
    }
};

// Two dummy classes for testing (de)serialization of a derived class through
// a pointer to an abstract base class
class KRATOS_API(KRATOS_CORE) AbstractTestClass
{
public:
    virtual ~AbstractTestClass() = default;
    [[nodiscard]] virtual int foo() const = 0;

private:
    friend class Kratos::Serializer;
    virtual void save(Serializer&) const = 0;
    virtual void load(Serializer&) = 0;

    // The following members are required to use this class with intrusive_ptr
    friend void intrusive_ptr_add_ref(const AbstractTestClass* pInstance)
    {
        if (pInstance) ++(pInstance->mRefCount);
    }

    friend void intrusive_ptr_release(const AbstractTestClass* pInstance)
    {
        if (pInstance) {
            --(pInstance->mRefCount);
            if (pInstance->mRefCount == 0) delete pInstance;
        }
    }

    mutable std::size_t mRefCount = 0; // Must be mutable, since previous two members receive a pointer-to-const
};

class KRATOS_API(KRATOS_CORE) DerivedTestClass : public AbstractTestClass
{
public:
    explicit DerivedTestClass(int FooNumber = 0) : mFooNumber(FooNumber) {}
    ~DerivedTestClass() override = default;
    [[nodiscard]] int foo() const override { return mFooNumber; }

private:
    friend class Kratos::Serializer;

    void save(Serializer& rSerializer) const override
    {
        rSerializer.save("mFooNumber", mFooNumber);
    }

    void load(Serializer& rSerializer) override
    {
        rSerializer.load("mFooNumber", mFooNumber);
    }

    int mFooNumber;
};

class KRATOS_API(KRATOS_CORE) ScopedTestClassRegistration
{
public:
    ScopedTestClassRegistration()
    {
        Serializer::Register("TestClass", TestClass{});
        Serializer::Register("DerivedTestClass", DerivedTestClass{});
    }

    ~ScopedTestClassRegistration()
    {
        Serializer::Deregister("TestClass");
        Serializer::Deregister("DerivedTestClass");
    }

    // Since the destructor has been defined, use the Rule of Five
    ScopedTestClassRegistration(const ScopedTestClassRegistration&) = delete;
    ScopedTestClassRegistration& operator=(const ScopedTestClassRegistration&) = delete;

    ScopedTestClassRegistration(ScopedTestClassRegistration&&) noexcept = default;
    ScopedTestClassRegistration& operator=(ScopedTestClassRegistration&&) noexcept = default;
};

}