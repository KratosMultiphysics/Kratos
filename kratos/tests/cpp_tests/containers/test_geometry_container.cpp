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
#include "geometries/line_3d_2.h"

#include "containers/geometry_container.h"

namespace Kratos {
namespace Testing {

    Line3D2<Point>::Pointer GenerateLineGeometry() {
        return Kratos::make_shared<Line3D2<Point>>(
            Kratos::make_shared<Point>(0.0, 0.0, 0.0),
            Kratos::make_shared<Point>(1.0, 1.0, 1.0)
            );
    }

    ///// Test Geometry Container
    KRATOS_TEST_CASE_IN_SUITE(TestgeometryContainer, KratosCoreGeometryContainerFastSuite) {
        auto geometry_container = GeometryContainer<Geometry<Point>>();

        auto p_line_1 = GenerateLineGeometry();
        p_line_1->SetId(1);

        geometry_container.AddGeometry(p_line_1);

        auto p_line_2 = GenerateLineGeometry();
        p_line_2->SetId(1);

        // check correct error if multiple geometries with sam id are added
        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            geometry_container.AddGeometry(p_line_2),
            "Geometry with Id: 1 exists already.");

        p_line_2->SetId(2);
        geometry_container.AddGeometry(p_line_2);

        // check correct number of geometries
        KRATOS_CHECK_EQUAL(geometry_container.NumberOfGeometries(), 2);

        // check adding with string
        auto p_line_3 = GenerateLineGeometry();
        p_line_3->SetId("GeometryLine1");
        geometry_container.AddGeometry(p_line_3);

        KRATOS_CHECK_EQUAL(geometry_container.NumberOfGeometries(), 3);

        // check if correct element is returned
        KRATOS_CHECK_EQUAL(geometry_container.GetGeometry(1).Id(), 1);
        KRATOS_CHECK_EQUAL(geometry_container.pGetGeometry(1)->Id(), 1);

        // check remove functions
        geometry_container.RemoveGeometry("GeometryLine1");
        geometry_container.RemoveGeometry(1);
        KRATOS_CHECK_EQUAL(geometry_container.NumberOfGeometries(), 1);
        KRATOS_CHECK_IS_FALSE(geometry_container.HasGeometry("GeometryLine1"));
    }
} // namespace Testing.
} // namespace Kratos.
