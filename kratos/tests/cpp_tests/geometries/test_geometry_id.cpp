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

    template<class TPointType>
    Quadrilateral2D4<TPointType>::Pointer
        GenerateQuadrilateralPointer(IndexType GeometryId)
    {
        Quadrilateral2D4<3, PointerVector<NodeType>>::Pointer p_quad =
            new Quadrilateral2D4<TPointType>(
                GeneratePoint<TPointType>(0.0, 0.0, 0.0),
                GeneratePoint<TPointType>(1.0, 0.0, 0.0),
                GeneratePoint<TPointType>(1.0, 1.0, 0.0),
                GeneratePoint<TPointType>(0.0, 1.0, 0.0)
                );
        p_quad.SetId(GeometryId);
        return p_quad;
    }

    ///// Test Geometry Id
    KRATOS_TEST_CASE_IN_SUITE(TestGeometryId, KratosCoreGeometriesFastSuite) {
        auto p_quad = GenerateQuadrilateralPointer<Point>(1);

        //// Check general information, input to ouput
        KRATOS_CHECK_EQUAL(p_quad->Id(), 1);
    }
} // namespace Testing.
} // namespace Kratos.
