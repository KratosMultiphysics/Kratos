//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/quadrilateral_2d_4.h"
#include "tests/geometries/test_geometry.h"

// Utility includes
#include "utilities/geometry_utilities.h"

namespace Kratos
{
namespace Testing
{
    /// Factory functions

    /** Generates a sample Quadrilateral2D4
     * Generates a quadrilateral defined by four random
     * points in the space.
     * @return  Pointer to a Quadrilateral2D4
     */
    template<class TPointType>
    typename Quadrilateral2D4<TPointType>::Pointer GenerateQuadrilateral2D4(
        typename TPointType::Pointer PointA = GeneratePoint<TPointType>(),
        typename TPointType::Pointer PointB = GeneratePoint<TPointType>(),
        typename TPointType::Pointer PointC = GeneratePoint<TPointType>(),
        typename TPointType::Pointer PointD = GeneratePoint<TPointType>()
        )
    {
      return typename Quadrilateral2D4<TPointType>::Pointer(new Quadrilateral2D4<TPointType>(
        PointA,
        PointB,
        PointC,
        PointD
      ));
    }

    /** Generates a sample Quadrilateral2D4
     * Generates a regular quadrilateral on the origin
     * @return  Pointer to a Quadrilateral2D4
     */
    template<class TPointType>
    typename Quadrilateral2D4<TPointType>::Pointer GenerateRightQuadrilateral2D4()
    {
      return typename Quadrilateral2D4<TPointType>::Pointer(new Quadrilateral2D4<TPointType>(
        GeneratePoint<TPointType>( 0.0, 0.0, 0.0),
        GeneratePoint<TPointType>( 1.0, 0.0, 0.0),
        GeneratePoint<TPointType>( 1.1, 1.1, 0.0),
        GeneratePoint<TPointType>( 0.0, 1.1, 0.0)
      ));
    }

    /** Generates a sample Quadrilateral2D4
     * Generates a regular quadrilateral centered on the origin
     * and 45deg rotated
     * @return  Pointer to a Quadrilateral2D4
     */
    template<class TPointType>
    typename Quadrilateral2D4<TPointType>::Pointer GenerateDiagQuadrilateral2D4()
    {
      return typename Quadrilateral2D4<TPointType>::Pointer(new Quadrilateral2D4<TPointType>(
        GeneratePoint<TPointType>( 1.0, 0.0, 0.0),
        GeneratePoint<TPointType>( 0.0, 1.0, 0.0),
        GeneratePoint<TPointType>(-1.0, 0.0, 0.0),
        GeneratePoint<TPointType>( 0.0,-1.0, 0.0)
      ));
    }

    /** Test a box and quadrilateral HasIntersection which should give true
     */
    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral2D4NodeBoxIntersection, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateRightQuadrilateral2D4<Node<3>>();
        Point point_1 (-0.3, 0.8, 0.0);
        Point point_2 ( 0.2, 1.5, 0.0);
        KRATOS_CHECK(geom->HasIntersection(point_1, point_2));
        
        Point point_3 ( 0.9, 1.5, 0.0);
        Point point_4 ( 1.1, 0.8, 0.0);
        KRATOS_CHECK(geom->HasIntersection(point_3, point_4));
        
        Point point_5 ( 0.9,-0.8, 0.0);
        Point point_6 ( 1.1, 0.1, 0.0);
        KRATOS_CHECK(geom->HasIntersection(point_5, point_6));
        
        Point point_7 (-0.3, 0.1, 0.0);
        Point point_8 ( 0.2,-0.6, 0.0);
        KRATOS_CHECK(geom->HasIntersection(point_7, point_8));
    }

    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral2D4EdgeBoxIntersection, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateDiagQuadrilateral2D4<Node<3>>();
        Point point_1 ( 0.2, 0.2, 0.0 );
        Point point_2 ( 1.0, 1.0, 0.0 );
        KRATOS_CHECK(geom->HasIntersection(point_1, point_2));
        
        Point point_3 (-0.2, 0.2, 0.0 );
        Point point_4 (-0.9, 0.9, 0.0 );
        KRATOS_CHECK(geom->HasIntersection(point_3, point_4));
        
        Point point_5 (-0.2,-0.2, 0.0 );
        Point point_6 (-0.9,-0.9, 0.0 );
        KRATOS_CHECK(geom->HasIntersection(point_5, point_6));
        
        Point point_7 ( 0.2,-0.2, 0.0 );
        Point point_8 ( 1.0,-1.0, 0.0 );
        KRATOS_CHECK(geom->HasIntersection(point_7, point_8));
    }

    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral2D4BoxNoIntersection, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateDiagQuadrilateral2D4<Node<3>>();
        Point point_1 ( 0.7, 0.4, 0.0 );
        Point point_2 ( 1.0, 1.2, 0.0 );
        KRATOS_CHECK_IS_FALSE(geom->HasIntersection(point_1, point_2));
    }

} // namespace Testing
} // namespace Kratos
