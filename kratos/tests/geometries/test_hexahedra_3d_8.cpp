//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//
//

// System includes
#include <set>

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/hexahedra_3d_8.h"
#include "tests/geometries/test_geometry.h"
#include "tests/geometries/test_shape_function_derivatives.h"
#include "tests/geometries/cross_check_shape_functions_values.h"

namespace Kratos {
	namespace Testing {

    typedef Node<3>                   PointType;
    typedef Node<3>::Pointer          PointPtrType;
    typedef Hexahedra3D8<PointType>   HexaGeometryType;
    typedef HexaGeometryType::Pointer HexaGeometryPtrType;

    /** Generates a sample Hexahedra3D8.
     * Generates a hexahedra defined by eight random points in the space.
     * @return  Pointer to a Hexahedra3D8
     */
    HexaGeometryPtrType GenerateHexahedra3D8(
        PointPtrType PointA = GeneratePoint<PointType>(),
        PointPtrType PointB = GeneratePoint<PointType>(),
        PointPtrType PointC = GeneratePoint<PointType>(),
        PointPtrType PointD = GeneratePoint<PointType>(),
        PointPtrType PointE = GeneratePoint<PointType>(),
        PointPtrType PointF = GeneratePoint<PointType>(),
        PointPtrType PointG = GeneratePoint<PointType>(),
        PointPtrType PointH = GeneratePoint<PointType>()) {
      return HexaGeometryPtrType(new HexaGeometryType(PointA, PointB, PointC, PointD, PointE, PointF, PointG, PointH));
    }

    /** Generates a sample Hexahedra3D8.
     * Generates a hexahedra with center on the origin with positive volume and side 1.
     * @return  Pointer to a Hexahedra3D8
     */
    HexaGeometryPtrType GenerateOriginCenterLen1Hexahedra3D8() {
      return HexaGeometryPtrType(new HexaGeometryType(
        GeneratePoint<PointType>(-0.5, -0.5, -0.5),
        GeneratePoint<PointType>( 0.5, -0.5, -0.5),
        GeneratePoint<PointType>( 0.5,  0.5, -0.5),
        GeneratePoint<PointType>(-0.5,  0.5, -0.5),
        GeneratePoint<PointType>(-0.5, -0.5,  0.5),
        GeneratePoint<PointType>( 0.5, -0.5,  0.5),
        GeneratePoint<PointType>( 0.5,  0.5,  0.5),
        GeneratePoint<PointType>(-0.5,  0.5,  0.5)
      ));
    }

    /** Generates a sample Hexahedra3D8.
     * Generates a hexahedra with center on the origin with positive volume and side 2.
     * @return  Pointer to a Hexahedra3D8
     */
    HexaGeometryPtrType GenerateOriginCenterLen2Hexahedra3D8() {
      return HexaGeometryPtrType(new HexaGeometryType(
        GeneratePoint<PointType>(-1.0, -1.0, -1.0),
        GeneratePoint<PointType>( 1.0, -1.0, -1.0),
        GeneratePoint<PointType>( 1.0,  1.0, -1.0),
        GeneratePoint<PointType>(-1.0,  1.0, -1.0),
        GeneratePoint<PointType>(-1.0, -1.0,  1.0),
        GeneratePoint<PointType>( 1.0, -1.0,  1.0),
        GeneratePoint<PointType>( 1.0,  1.0,  1.0),
        GeneratePoint<PointType>(-1.0,  1.0,  1.0)
      ));
    }

    /** Checks if the number of edges is correct.
     * Checks if the number of edges is correct.
     */
    KRATOS_TEST_CASE_IN_SUITE(Hexahedra3D8EdgesNumber, KratosCoreGeometriesFastSuite) {
      auto geomRegLen1 = GenerateOriginCenterLen1Hexahedra3D8();
      auto geomRegLen2 = GenerateOriginCenterLen2Hexahedra3D8();

      KRATOS_CHECK_EQUAL(geomRegLen1->EdgesNumber(), 12);
      KRATOS_CHECK_EQUAL(geomRegLen2->EdgesNumber(), 12);
    }

    /** Checks if the number of faces is correct.
     * Checks if the number of faces is correct.
     */
    KRATOS_TEST_CASE_IN_SUITE(Hexahedra3D8FacesNumber, KratosCoreGeometriesFastSuite) {
      auto geomRegLen1 = GenerateOriginCenterLen1Hexahedra3D8();
      auto geomRegLen2 = GenerateOriginCenterLen2Hexahedra3D8();

      KRATOS_CHECK_EQUAL(geomRegLen1->FacesNumber(), 6);
      KRATOS_CHECK_EQUAL(geomRegLen2->FacesNumber(), 6);
    }

    /** Checks if the characteristic length of the hexahedra is calculated correctly.
     * Checks if the characteristic length of the hexahedra is calculated correctly.
     */
    KRATOS_TEST_CASE_IN_SUITE(Hexahedra3D8Length, KratosCoreGeometriesFastSuite) {
      auto geomRegLen1 = GenerateOriginCenterLen1Hexahedra3D8();
      auto geomRegLen2 = GenerateOriginCenterLen2Hexahedra3D8();

      KRATOS_CHECK_NEAR(geomRegLen1->Length(), 0.353553, TOLERANCE);
      KRATOS_CHECK_NEAR(geomRegLen2->Length(), 1.000000, TOLERANCE);
    }

    /** Checks if the area of the hexahedra is calculated correctly.
     * Checks if the area of the hexahedra is calculated correctly.
     */
    KRATOS_TEST_CASE_IN_SUITE(Hexahedra3D8Area, KratosCoreGeometriesFastSuite) {
      auto geomRegLen1 = GenerateOriginCenterLen1Hexahedra3D8();
      auto geomRegLen2 = GenerateOriginCenterLen2Hexahedra3D8();

      KRATOS_CHECK_NEAR(geomRegLen1->Area(), 1.0, TOLERANCE);
      KRATOS_CHECK_NEAR(geomRegLen2->Area(), 8.0, TOLERANCE);
    }

    /** Checks if the volume of the hexahedra is calculated correctly.
     * Checks if the volume of the hexahedra is calculated correctly.
     * For hexahedra 3D8 'volume()' call defaults to 'area()'
     */
    KRATOS_TEST_CASE_IN_SUITE(Hexahedra3D8Volume, KratosCoreGeometriesFastSuite) {
      auto geomRegLen1 = GenerateOriginCenterLen1Hexahedra3D8();
      auto geomRegLen2 = GenerateOriginCenterLen2Hexahedra3D8();

      KRATOS_CHECK_NEAR(geomRegLen1->Volume(), 1.0, TOLERANCE);
      KRATOS_CHECK_NEAR(geomRegLen2->Volume(), 8.0, TOLERANCE);
  	}

    /**
     * This test performs the check of the box intersection method
     */
    KRATOS_TEST_CASE_IN_SUITE(Hexahedra3D8BoxIntersection, KratosCoreGeometriesFastSuite) {
      auto hexahedron = GenerateOriginCenterLen1Hexahedra3D8();

      //hexahedron inside the box
      KRATOS_CHECK(hexahedron->HasIntersection(Point(-0.6,-0.6,-0.6), Point(0.6,0.6,0.6)));

      //hexahedron contains the box
      KRATOS_CHECK(hexahedron->HasIntersection(Point(-.25,-.25,-.25), Point(.25,.25,.25)));

      //hexahedron intersects the box
      KRATOS_CHECK(hexahedron->HasIntersection(Point(.25,.25,.25), Point(1.0,1.0,1.0)));

      //hexahedron not intersects the box
      KRATOS_CHECK_IS_FALSE(hexahedron->HasIntersection(Point(.51,.51,.51), Point(1.1,1.1,1.1)));
    }
    
    /** Checks the inside test for a given point respect to the hexahedra
    * Checks the inside test for a given point respect to the hexahedra
    * It performs 4 tests:
    * A Point inside the hexahedra: Expected result TRUE
    * A Point outside the hexahedra: Expected result FALSE
    * A Point over a vertex of the hexahedra: Expected result TRUE
    * A Point over an edge of the hexahedra: Expected result TRUE
    */
    KRATOS_TEST_CASE_IN_SUITE(Hexahedra3D8IsInside, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateOriginCenterLen1Hexahedra3D8();

        Point PointInside(0.4999, 0.4999, 0.4999);
        Point PointOutside(0.5001, 0.5001, 0.5001);
        Point PointInVertex(-0.5, -0.5, -0.5);
        Point PointInEdge(0.5, 0.5, 0.0);

        Point LocalCoords;
        
        KRATOS_CHECK(geom->IsInside(PointInside, LocalCoords, EPSILON));
        KRATOS_CHECK_IS_FALSE(geom->IsInside(PointOutside, LocalCoords, EPSILON));
        KRATOS_CHECK(geom->IsInside(PointInVertex, LocalCoords, EPSILON));
        KRATOS_CHECK(geom->IsInside(PointInEdge, LocalCoords, EPSILON));
    }

    /** Checks the point local coordinates for a given point respect to the
    * hexahedra. The baricentre of the hexahedra is selected due to its known
    * solution.
    */
    KRATOS_TEST_CASE_IN_SUITE(Hexahedra3D8PointLocalCoordinates, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateOriginCenterLen2Hexahedra3D8();

        // Compute the global coordinates of the baricentre
        auto points = geom->Points();
        Point baricentre = points[0] + points[1] + points[2] + points[3] + points[4] + points[5] + points[6] + points[7];
        baricentre /= 8.0;

        // Compute the baricentre local coordinates
        array_1d<double, 3> baricentre_local_coords;
        geom->PointLocalCoordinates(baricentre_local_coords, baricentre);

        KRATOS_CHECK_NEAR(baricentre_local_coords(0), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(baricentre_local_coords(1), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(baricentre_local_coords(2), 0.0, TOLERANCE);
        
        Point baricentre_face_1;
        baricentre_face_1.Coordinates()[0] = 1.0;
        baricentre_face_1.Coordinates()[1] = 0.0;
        baricentre_face_1.Coordinates()[2] = 0.0;

        // Compute the baricentre local coordinates
        array_1d<double, 3> baricentre_local_coords_face_1;
        geom->PointLocalCoordinates(baricentre_local_coords_face_1, baricentre_face_1);
        
        KRATOS_CHECK_NEAR(baricentre_local_coords_face_1(0), 1.0, TOLERANCE);
        KRATOS_CHECK_NEAR(baricentre_local_coords_face_1(1), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(baricentre_local_coords_face_1(2), 0.0, TOLERANCE);
        
        Point baricentre_face_2;
        baricentre_face_2.Coordinates()[0] =  0.0;
        baricentre_face_2.Coordinates()[1] = -1.0;
        baricentre_face_2.Coordinates()[2] =  0.0;

        // Compute the baricentre local coordinates
        array_1d<double, 3> baricentre_local_coords_face_2;
        geom->PointLocalCoordinates(baricentre_local_coords_face_2, baricentre_face_2);
        
        KRATOS_CHECK_NEAR(baricentre_local_coords_face_2(0), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(baricentre_local_coords_face_2(1),-1.0, TOLERANCE);
        KRATOS_CHECK_NEAR(baricentre_local_coords_face_2(2), 0.0, TOLERANCE);
        
        Point baricentre_face_3;
        baricentre_face_3.Coordinates()[0] = 0.0;
        baricentre_face_3.Coordinates()[1] = 0.3;
        baricentre_face_3.Coordinates()[2] = 1.0;

        // Compute the baricentre local coordinates
        array_1d<double, 3> baricentre_local_coords_face_3;
        geom->PointLocalCoordinates(baricentre_local_coords_face_3, baricentre_face_3);

        KRATOS_CHECK_NEAR(baricentre_local_coords_face_3(0), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(baricentre_local_coords_face_3(1), 0.3, TOLERANCE);
        KRATOS_CHECK_NEAR(baricentre_local_coords_face_3(2), 1.0, TOLERANCE);

        Point outside_point;
        outside_point.Coordinates()[0] = 1.0;
        outside_point.Coordinates()[1] = 1.0;
        outside_point.Coordinates()[2] = 1.0;

        // Compute the baricentre local coordinates
        array_1d<double, 3> local_coords_outside_point;
        geom->PointLocalCoordinates(local_coords_outside_point, outside_point);

        KRATOS_CHECK_NEAR(local_coords_outside_point(0), 1.0, TOLERANCE);
        KRATOS_CHECK_NEAR(local_coords_outside_point(1), 1.0, TOLERANCE);
        KRATOS_CHECK_NEAR(local_coords_outside_point(2), 1.0, TOLERANCE);
    }

  KRATOS_TEST_CASE_IN_SUITE(Hexahedra3D8ShapeFunctionsValues, KratosCoreGeometriesFastSuite) {
      auto geom = GenerateOriginCenterLen2Hexahedra3D8();
      array_1d<double, 3> coord(3);
      coord[0] = 1.0 / 2.0;
      coord[1] = 1.0 / 4.0;
      coord[2] = 1.0 / 16.0;
      KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(0, coord), 0.0439453125, TOLERANCE);
      KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(1, coord), 0.1318359375, TOLERANCE);
      KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(2, coord), 0.2197265625, TOLERANCE);
      KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(3, coord), 0.0732421875, TOLERANCE);
      KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(4, coord), 0.0498046875, TOLERANCE);
      KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(5, coord), 0.1494140625, TOLERANCE);
      KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(6, coord), 0.2490234375, TOLERANCE);
      KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(7, coord), 0.0830078125, TOLERANCE);
      CrossCheckShapeFunctionsValues(*geom);
  }

  KRATOS_TEST_CASE_IN_SUITE(Hexahedra3D8ShapeFunctionsLocalGradients, KratosCoreGeometriesFastSuite) {
      auto geom = GenerateOriginCenterLen2Hexahedra3D8();
      TestAllShapeFunctionsLocalGradients(*geom);
  }
	}
}  // namespace Kratos.
