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

#include "containers/model.h"

namespace Kratos {
namespace Testing {

    Line3D2<Node<3>>::Pointer GenerateLineModelPartGeometryContainer() {
        PointerVector<Node<3>> points;

        points.push_back(Node<3>::Pointer(new Node<3>(1, 0, 5, 0)));
        points.push_back(Node<3>::Pointer(new Node<3>(2, 5, 5, 0)));

        return Kratos::make_shared<Line3D2<Node<3>>>(
            points
            );
    }

    ///// Test Geometry Container
    KRATOS_TEST_CASE_IN_SUITE(TestModelPartGeometryContainer, KratosCoreGeometryContainerFastSuite) {
        Model model;

        auto& model_part = model.CreateModelPart("Main");

        auto& model_part_lines = model_part.CreateSubModelPart("Lines");

        auto& model_part_no_lines = model_part.CreateSubModelPart("NoLines");

        auto p_line_1 = GenerateLineModelPartGeometryContainer();
        p_line_1->SetId(1);

        model_part_lines.AddGeometry(p_line_1);

        auto p_line_2 = GenerateLineModelPartGeometryContainer();
        p_line_2->SetId(1);

        // check correct error if multiple geometries with sam id are added
        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            model_part_lines.AddGeometry(p_line_2),
            "Geometry with Id: 1 exists already.");

        p_line_2->SetId(2);
        model_part_lines.AddGeometry(p_line_2);

        // check correct number of geometries
        KRATOS_CHECK_EQUAL(model_part_lines.NumberOfGeometries(), 2);
        KRATOS_CHECK_EQUAL(model_part.NumberOfGeometries(), 2);
        KRATOS_CHECK_EQUAL(model_part_no_lines.NumberOfGeometries(), 0);

        // check adding with string
        auto p_line_3 = GenerateLineModelPartGeometryContainer();
        p_line_3->SetId("GeometryLine1");
        model_part_lines.AddGeometry(p_line_3);

        KRATOS_CHECK_EQUAL(model_part_lines.NumberOfGeometries(), 3);

        // check if correct element is returned
        KRATOS_CHECK_EQUAL(model_part_lines.GetGeometry(1).Id(), 1);
        KRATOS_CHECK_EQUAL(model_part_lines.pGetGeometry(1)->Id(), 1);

        // check if correct element is returned from root
        KRATOS_CHECK_EQUAL(model_part.GetGeometry(1).Id(), 1);
        KRATOS_CHECK_EQUAL(model_part.pGetGeometry(1)->Id(), 1);

        // check remove functions
        model_part_lines.RemoveGeometry("GeometryLine1");
        model_part_lines.RemoveGeometryFromAllLevels(1);
        KRATOS_CHECK_EQUAL(model_part_lines.NumberOfGeometries(), 1);
        KRATOS_CHECK_EQUAL(model_part.NumberOfGeometries(), 2);
        KRATOS_CHECK_EQUAL(model_part.NumberOfGeometries(), 2);

        // check if correct geometries are removed
        KRATOS_CHECK_IS_FALSE(model_part_lines.HasGeometry("GeometryLine1"));
        KRATOS_CHECK(model_part.HasGeometry("GeometryLine1"));

        // check if correct geometries are removed
        KRATOS_CHECK_IS_FALSE(model_part_lines.HasGeometry(1));
        KRATOS_CHECK_IS_FALSE(model_part.HasGeometry(1));
    }
} // namespace Testing.
} // namespace Kratos.
