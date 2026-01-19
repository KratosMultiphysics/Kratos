
#include "testing/testing.h"
#include "includes/serializer.h"
#include "includes/stream_serializer.h"
#include "includes/define.h"

namespace Kratos::Testing 
{

class SerializerUniquePtrAliasingTestClass {
public:
    int mValue;
    SerializerUniquePtrAliasingTestClass() : mValue(0) {}
    SerializerUniquePtrAliasingTestClass(int v) : mValue(v) {}
    virtual ~SerializerUniquePtrAliasingTestClass() = default;
private:
    friend class Kratos::Serializer;
    void save(Serializer& rSerializer) const {
        rSerializer.save("Value", mValue);
    }
    void load(Serializer& rSerializer) {
        rSerializer.load("Value", mValue);
    }
};

class ScopedTestClassRegistration
{
public:
    ScopedTestClassRegistration()
    {
        Serializer::Register("SerializerUniquePtrAliasingTestClass", SerializerUniquePtrAliasingTestClass{});
    }

    ~ScopedTestClassRegistration()
    {
        Serializer::Deregister("SerializerUniquePtrAliasingTestClass");
    }

    ScopedTestClassRegistration(const ScopedTestClassRegistration&) = delete;
    ScopedTestClassRegistration& operator=(const ScopedTestClassRegistration&) = delete;
    ScopedTestClassRegistration(ScopedTestClassRegistration&&) noexcept = default;
    ScopedTestClassRegistration& operator=(ScopedTestClassRegistration&&) noexcept = default;
};

KRATOS_TEST_CASE_IN_SUITE(SerializerUniquePtrAliasing, KratosCoreFastSuite)
{
    StreamSerializer serializer;
    ScopedTestClassRegistration scoped_registration;
    const std::string tag_string("TestString");

    // Save: RawPtr first, then UniquePtr
    {
        auto p_unique = Kratos::make_unique<SerializerUniquePtrAliasingTestClass>(42);
        SerializerUniquePtrAliasingTestClass* p_raw = p_unique.get();
        
        serializer.save("raw", p_raw);
        serializer.save("unique", p_unique);
    }

    // Load
    {
        Kratos::unique_ptr<SerializerUniquePtrAliasingTestClass> p_unique_loaded;
        SerializerUniquePtrAliasingTestClass* p_raw_loaded = nullptr;

        serializer.SetLoadState(); // Reset for loading
        
        serializer.load("raw", p_raw_loaded);
        serializer.load("unique", p_unique_loaded);

        KRATOS_EXPECT_NE(p_unique_loaded, nullptr);
        if (p_unique_loaded) {
            KRATOS_EXPECT_EQ(p_unique_loaded->mValue, 42);
        }

        // Check Aliasing - this should now work with the fix
        KRATOS_EXPECT_EQ(p_raw_loaded, p_unique_loaded.get());
        
        if (p_raw_loaded) {
            KRATOS_EXPECT_EQ(p_raw_loaded->mValue, 42);
        }
    }
}
}
