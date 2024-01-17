// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Gennady Markelov
//
#pragma once

#include "testing/testing.h"
#include "geo_mechanics_application.h"

namespace Kratos::Testing
{
    class KratosGeoMechanicsSuite : public KratosCoreFastSuite {

    protected:
        // Per-test-suite set-up.
        // Called before the first test in this test suite.
        // Can be omitted if not needed.
        void SetUp() override {
            geo_mechanics_app = new KratosGeoMechanicsApplication();
            geo_mechanics_app->Register();

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
            delete geo_mechanics_app;
            geo_mechanics_app = nullptr;
        }

        KratosGeoMechanicsApplication * geo_mechanics_app;

    };

    class KratosGeoMechanicsFastSuite : public KratosGeoMechanicsSuite {};
    class KratosGeoMechanicsIntegrationSuite : public KratosGeoMechanicsSuite {};
}