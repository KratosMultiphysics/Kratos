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

#pragma once

// System includes

// External includes

// Project includes
#include "includes/serializer.h"

namespace Kratos::Testing
{
    class KRATOS_API(KRATOS_TEST_UTILS) TestClass 
    {
    public:
        int mValue;
        TestClass();
        TestClass(int v);
        virtual ~TestClass() = default;
    private:
        friend class Kratos::Serializer;
        void save(Serializer& rSerializer) const;
        void load(Serializer& rSerializer);
    };

    class KRATOS_API(KRATOS_TEST_UTILS) AbstractTestClass
    {
    public:
        virtual ~AbstractTestClass() = default;
        [[nodiscard]] virtual int foo() const = 0;

    private:
        friend class Kratos::Serializer;
        virtual void save(Serializer&) const = 0;
        virtual void load(Serializer&) = 0;

        friend void intrusive_ptr_add_ref(const AbstractTestClass* pInstance) {
            if (pInstance) ++(pInstance->mRefCount);
        }

        friend void intrusive_ptr_release(const AbstractTestClass* pInstance) {
            if (pInstance) {
                --(pInstance->mRefCount);
                if (pInstance->mRefCount == 0) delete pInstance;
            }
        }
        mutable std::size_t mRefCount = 0;
    };

    class KRATOS_API(KRATOS_TEST_UTILS) DerivedTestClass : public AbstractTestClass
    {
    public:
        explicit DerivedTestClass(int FooNumber = 0);
        ~DerivedTestClass() override = default;
        [[nodiscard]] int foo() const override;

    private:
        friend class Kratos::Serializer;
        void save(Serializer& rSerializer) const override;
        void load(Serializer& rSerializer) override;
        int mFooNumber;
    };

    class KRATOS_API(KRATOS_TEST_UTILS) ScopedTestClassRegistration
    {
    public:
        ScopedTestClassRegistration();
        ~ScopedTestClassRegistration();

        ScopedTestClassRegistration(const ScopedTestClassRegistration&) = delete;
        ScopedTestClassRegistration& operator=(const ScopedTestClassRegistration&) = delete;
        ScopedTestClassRegistration(ScopedTestClassRegistration&&) noexcept = default;
        ScopedTestClassRegistration& operator=(ScopedTestClassRegistration&&) noexcept = default;
    };

} // namespace Kratos::Testing