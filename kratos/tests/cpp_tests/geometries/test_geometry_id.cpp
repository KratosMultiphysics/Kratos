//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Teschemacher
//                   Pooyan Dadvand
//

// System includes
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/quadrilateral_2d_4.h"

#include "tests/cpp_tests/geometries/test_geometry.h"

#include "containers/pointer_vector.h"

namespace Kratos {
namespace Testing {

    Line3D2<Point>::Pointer GenerateGeometry() {
        return Kratos::make_shared<Line3D2<Point>>(
            Kratos::make_shared<Point>(0.0, 0.0, 0.0),
            Kratos::make_shared<Point>(1.0, 1.0, 1.0)
            );
    }

    ///// Test Geometry Id with string
    KRATOS_TEST_CASE_IN_SUITE(TestNameHash, KratosCoreGeometriesFastSuite) {
        auto test_id =Geometry<Point>::GetNameHash("GeometryTestName");

        KRATOS_WATCH(test_id)
    }

    ///// Test Geometry Id
    KRATOS_TEST_CASE_IN_SUITE(TestGeometryId, KratosCoreGeometriesFastSuite) {
        auto p_quad = GenerateGeometry();

        p_quad->SetId(1);

        KRATOS_CHECK_EQUAL(p_quad->Id(), 1);
    }

    ///// Test Geometry Id with name
    KRATOS_TEST_CASE_IN_SUITE(TestGeometryName, KratosCoreGeometriesFastSuite) {
        auto p_quad = GenerateGeometry();

        p_quad->SetId("Geometry1");

        KRATOS_CHECK_EXCEPTION_IS_THROWN(p_quad->Id(), "Assigned geometry id is name of type string.")

        KRATOS_CHECK_EQUAL(p_quad->GetHashId(), Geometry<Point>::GetNameHash("Geometry1"));
    }
} // namespace Testing.
} // namespace Kratos.
