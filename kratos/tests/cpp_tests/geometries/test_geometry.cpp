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
//                   Vicente Mataix Ferrandiz
//
//

// System includes
#include <set>

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/geometry.h"
#include "tests/cpp_tests/geometries/test_geometry.h"

namespace Kratos {
namespace Testing {

    /// Auxiliar check functions (from geometry_tester.h)
    /// - All this functions should probably me moved somewhere else.

    /// Test self assigned geometry Id
    TEST(GeometryIdSelfAssigned, KratosCoreGeometriesFastSuite) {
        auto this_geometry = Geometry<Point>();

        KRATOS_EXPECT_FALSE(this_geometry.IsIdGeneratedFromString());
        KRATOS_EXPECT_TRUE(this_geometry.IsIdSelfAssigned());

        this_geometry.SetId(2);
        KRATOS_EXPECT_FALSE(this_geometry.IsIdGeneratedFromString());
        KRATOS_EXPECT_FALSE(this_geometry.IsIdSelfAssigned());

        this_geometry.SetId("ThisGeometry");
        KRATOS_EXPECT_TRUE(this_geometry.IsIdGeneratedFromString());
        KRATOS_EXPECT_FALSE(this_geometry.IsIdSelfAssigned());
    }

    /// Test geometry Id with name
    TEST(GeometryName, KratosCoreGeometriesFastSuite) {
        auto this_geometry = Geometry<Point>("Geometry1");

        KRATOS_EXPECT_TRUE(this_geometry.IsIdGeneratedFromString());
        KRATOS_EXPECT_FALSE(this_geometry.IsIdSelfAssigned());
        KRATOS_EXPECT_EQ(this_geometry.Id(), Geometry<Point>::GenerateId("Geometry1"));
    }

    /// Test geometry Id
    TEST(GeometryId, KratosCoreGeometriesFastSuite) {
        auto this_geometry = Geometry<Point>(1);

        KRATOS_EXPECT_FALSE(this_geometry.IsIdGeneratedFromString());
        KRATOS_EXPECT_FALSE(this_geometry.IsIdSelfAssigned());
        KRATOS_EXPECT_EQ(this_geometry.Id(), 1);

        // Check for higher Id.
        auto this_geometry_2 = Geometry<Point>(717);
        KRATOS_EXPECT_FALSE(this_geometry_2.IsIdGeneratedFromString());
        KRATOS_EXPECT_FALSE(this_geometry_2.IsIdSelfAssigned());
        KRATOS_EXPECT_EQ(this_geometry_2.Id(), 717);
    }
} // namespace Testing.
} // namespace Kratos.
