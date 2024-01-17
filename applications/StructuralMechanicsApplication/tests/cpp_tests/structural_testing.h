// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Gennady Markelov
//
#pragma once

#include "testing/testing.h"
#include "structural_mechanics_application.h"

namespace Kratos::Testing
{
    class KratosStructuralMechanicsSuite : public KratosCoreFastSuite {

    protected:
        // Per-test-suite set-up.
        // Called before the first test in this test suite.
        // Can be omitted if not needed.
        void SetUp() override {
            structural_mechanics_app = new KratosStructuralMechanicsApplication();
            structural_mechanics_app->Register();

            // If `shared_resource_` is **not deleted** in `TearDownTestSuite()`,
            // reallocation should be prevented because `SetUpTestSuite()` may be called
            // in subclasses of FooTest and lead to memory leak.
            //
            // if (shared_resource_ == nullptr) {
            //   shared_resource_ = new ...;
            // }
        }

        // Per-test-suite tear-down.
        // Called after the last test in this test suite.
        // Can be omitted if not needed.
        void TearDown() override {
            delete structural_mechanics_app;
            structural_mechanics_app = nullptr;
        }

        KratosStructuralMechanicsApplication * structural_mechanics_app;

    };

    class KratosStructuralMechanicsFastSuite : public KratosStructuralMechanicsSuite {};
}