//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//
//

// System includes

// External includes

// Project includes
#include "tests/test_utilities/serializer_test_utilities.h"

namespace Kratos::Testing 
{
    // TestClass Implementation
    TestClass::TestClass() : mValue(0) {}
    TestClass::TestClass(int v) : mValue(v) {}

    void TestClass::save(Serializer& rSerializer) const {
        rSerializer.save("Value", mValue);
    }
    void TestClass::load(Serializer& rSerializer) {
        rSerializer.load("Value", mValue);
    }

    // DerivedTestClass Implementation
    DerivedTestClass::DerivedTestClass(int FooNumber) : mFooNumber(FooNumber) {}
    
    int DerivedTestClass::foo() const { 
        return mFooNumber; 
    }

    void DerivedTestClass::save(Serializer& rSerializer) const {
        rSerializer.save("mFooNumber", mFooNumber);
    }

    void DerivedTestClass::load(Serializer& rSerializer) {
        rSerializer.load("mFooNumber", mFooNumber);
    }

    // ScopedTestClassRegistration Implementation
    ScopedTestClassRegistration::ScopedTestClassRegistration() {
        Serializer::Register("TestClass", TestClass{});
        Serializer::Register("DerivedTestClass", DerivedTestClass{});
    }

    ScopedTestClassRegistration::~ScopedTestClassRegistration() {
        Serializer::Deregister("TestClass");
        Serializer::Deregister("DerivedTestClass");
    }

}  // namespace Kratos::Testing