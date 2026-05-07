//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos A. Roig
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/serializer.h"
#include "includes/stream_serializer.h"
#include "tests/test_utilities/serializer_testing_utilities.h"

namespace Kratos::Testing 
{

KRATOS_TEST_CASE_IN_SUITE(SerializerUniquePtrAliasing, KratosCoreFastSuite)
{
    StreamSerializer serializer;
    ScopedTestClassRegistration scoped_registration;
    const std::string tag_string("TestString");

    // Save: RawPtr first, then UniquePtr
    {
        auto p_unique = Kratos::make_unique<TestClass>(42);
        TestClass* p_raw = p_unique.get();
        
        serializer.save("raw", p_raw);
        serializer.save("unique", p_unique);
    }

    // Load
    {
        Kratos::unique_ptr<TestClass> p_unique_loaded;
        TestClass* p_raw_loaded = nullptr;

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
