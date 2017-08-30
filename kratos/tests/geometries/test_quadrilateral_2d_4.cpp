//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Masó Sotomayor
//
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/quadrilateral_2d_4.h"
#include "tests/geometries/test_geometry.h"

namespace Kratos {
  namespace Testing {

    typedef Node<3>                      PointType;
    typedef Node<3>::Pointer             PointPtrType;
    typedef Quadrilateral2D4<PointType>  GeometryType;
    typedef GeometryType::Pointer        GeometryPtrType;

    /** Generates a sample Quadrilateral2D4
     * Generates a regular quadrilateral on the origin
     * @return  Pointer to a Quadrilateral2D4
     */
    GeometryPtrType GenerateRegQuadrilateral2D4() {
      return GeometryPtrType(new GeometryType(
        GeneratePoint<PointType>( 0.0, 0.0, 0.0),
        GeneratePoint<PointType>( 1.0, 0.0, 0.0),
        GeneratePoint<PointType>( 1.1, 1.1, 0.0),
        GeneratePoint<PointType>( 0.0, 1.1, 0.0)
      ));
    }
    
    /** Generates a sample Quadrilateral2D4
     * Generates a regular quadrilateral centered on the origin and 45deg rotated
     * @return  Pointer to a Quadrilateral2D4
     */
    GeometryPtrType GenerateDiagQuadrilateral2D4() {
      return GeometryPtrType(new GeometryType(
        GeneratePoint<PointType>( 1.0, 0.0, 0.0),
        GeneratePoint<PointType>( 0.0, 1.0, 0.0),
        GeneratePoint<PointType>(-1.0, 0.0, 0.0),
        GeneratePoint<PointType>( 0.0,-1.0, 0.0)
      ));
    }

    /** Test a box and quadrilateral HasIntersection which should give true
     */
    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral2D4NodeIntersection, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateRegQuadrilateral2D4();
        PointType point_1 (-0.5, 0.5, 0.0);
        PointType point_2 ( 0.5, 1.5, 0.0);
        KRATOS_CHECK(geom->HasIntersection(point_1, point_2));
    }

    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral2D4EdgeIntersection, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateDiagQuadrilateral2D4();
        PointType point_1 ( 0.5, 0.5, 0.0);
        PointType point_2 ( 1.0, 1.0, 0.0);
        KRATOS_CHECK(geom->HasIntersection(point_1, point_2));
    }

  } // namespace Testing
} // namespace Kratos
