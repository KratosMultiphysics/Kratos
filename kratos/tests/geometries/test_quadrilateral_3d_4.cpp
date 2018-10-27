//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/quadrilateral_3d_4.h"
#include "tests/geometries/test_geometry.h"
#include "tests/geometries/test_shape_function_derivatives.h"
#include "tests/geometries/cross_check_shape_functions_values.h"

// Utility includes
#include "utilities/geometry_utilities.h"

namespace Kratos 
{
namespace Testing 
{
    ///@name Type Definitions
    ///@{

    typedef Node<3> NodeType;

    ///@}

    /// Factory functions

    /** Generates a sample Quadrilateral3D4.
    * Generates a quadrilateral defined by three random points in the space.
    * @return  Pointer to a Quadrilateral3D4
    */
    template<class TPointType>
    typename Quadrilateral3D4<TPointType>::Pointer GenerateQuadrilateral3D4(
        typename TPointType::Pointer PointA = GeneratePoint<TPointType>(),
        typename TPointType::Pointer PointB = GeneratePoint<TPointType>(),
        typename TPointType::Pointer PointC = GeneratePoint<TPointType>(),
        typename TPointType::Pointer PointD = GeneratePoint<TPointType>()
        ) 
    {
        return typename Quadrilateral3D4<TPointType>::Pointer(new Quadrilateral3D4<TPointType>(
            PointA,
            PointB,
            PointC,
            PointD
            ));
    }

    /** Generates a sample Quadrilateral3D4.
    * Generates a right quadrilateral with origin in the origin and leg size 1.
    * @return  Pointer to a Quadrilateral3D4
    */
    template<class TPointType>
    typename Quadrilateral3D4<TPointType>::Pointer GenerateRightQuadrilateral3D4() 
    {
        return typename Quadrilateral3D4<TPointType>::Pointer(new Quadrilateral3D4<TPointType>(
        GeneratePoint<TPointType>(0.0, 0.0, 0.0),
        GeneratePoint<TPointType>(std::cos(Globals::Pi/4), 0.0, std::sin(Globals::Pi/4)),
        GeneratePoint<TPointType>(1.0, 1.0, 0.5),
        GeneratePoint<TPointType>(0.0, 1.0, 0.0)
        ));
    }

    /** Generates a sample Quadrilateral3D4.
    * Generates a right quadrilateral with origin in the origin and leg size 1.
    * @return  Pointer to a Quadrilateral3D4
    */
    template<class TPointType>
    typename Quadrilateral3D4<TPointType>::Pointer GenerateFlatQuadrilateral3D4() 
    {
        return typename Quadrilateral3D4<TPointType>::Pointer(new Quadrilateral3D4<TPointType>(
        GeneratePoint<TPointType>( 0.0, 0.0, 0.0),
        GeneratePoint<TPointType>( 1.0, 0.0, 0.0),
        GeneratePoint<TPointType>( 1.0, 1.0, 0.0),
        GeneratePoint<TPointType>( 0.0, 1.0, 0.0)
        ));
    }

    /// Tests

    /** Checks if the number of edges is correct.
    * Checks if the number of edges is correct.
    */
    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4EdgesNumber, KratosCoreFastSuite)
    {
        auto geom = GenerateRightQuadrilateral3D4<NodeType>();

        KRATOS_CHECK_EQUAL(geom->EdgesNumber(), 4);
    }

    /** Checks if the number of faces is correct.
    * Checks if the number of faces is correct.
    */
    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4FacesNumber, KratosCoreFastSuite)
    {
        auto geom = GenerateRightQuadrilateral3D4<NodeType>();

        // That for planar geometries it also return the number of edges.
        KRATOS_CHECK_EQUAL(geom->FacesNumber(), 4);
    }

    /** Checks if the area of the quadrilateral is calculated correctly.
    * Checks if the area of the quadrilateral is calculated correctly.
    */
    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4Area, KratosCoreFastSuite)
    {
        auto geom = GenerateRightQuadrilateral3D4<NodeType>();
        
        KRATOS_CHECK_NEAR(geom->Area(), 1.06947235, TOLERANCE);
//         KRATOS_CHECK_NEAR(geom->Area(), 1.08935, TOLERANCE); // NOTE: Solution from Mathematica
    }

    /** Tests the PointLocalCoordinates for Quadrilateral2D4.
     * Tests the PointLocalCoordinates for Quadrilateral2D4.
     */
    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4PointLocalCoordinates, KratosCoreFastSuite) {
        auto geom = GenerateFlatQuadrilateral3D4<NodeType>();

        Point TestPointA(1.0, 1.0, 0.0);
        Point TestPointB(0.5, 0.5, 0.0);
        Point TestResultA(0.0, 0.0, 0.0);
        Point TestResultB(0.0, 0.0, 0.0);

        geom->PointLocalCoordinates(TestResultA, TestPointA);
        geom->PointLocalCoordinates(TestResultB, TestPointB);

        // Test transformation in the edge
        KRATOS_CHECK_NEAR(TestResultA[0], 1.0, TOLERANCE);
        KRATOS_CHECK_NEAR(TestResultA[1], 1.0, TOLERANCE);
        KRATOS_CHECK_NEAR(TestResultA[2], 0.0, TOLERANCE);

        // Test transformation in the center
        KRATOS_CHECK_NEAR(TestResultB[0], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(TestResultB[1], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(TestResultB[2], 0.0, TOLERANCE);
    }

    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4ShapeFunctionsValues, KratosCoreFastSuite) {
      auto geom = GenerateFlatQuadrilateral3D4<NodeType>();
      array_1d<double, 3> coord(3);
      coord[0] = 1.0 / 2.0;
      coord[1] = 1.0 / 4.0;
      coord[2] = 0.0;
      KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(0, coord), 0.09375, TOLERANCE);
      KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(1, coord), 0.28125, TOLERANCE);
      KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(2, coord), 0.46875, TOLERANCE);
      KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(3, coord), 0.15625, TOLERANCE);
      CrossCheckShapeFunctionsValues(*geom);
    }

    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4ShapeFunctionsLocalGradients, KratosCoreFastSuite) {
        auto geom = GenerateFlatQuadrilateral3D4<NodeType>();
        TestAllShapeFunctionsLocalGradients(*geom);
    }

//     /** Checks if the volume of the quadrilateral is calculated correctly.
//     * Checks if the volume of the quadrilateral is calculated correctly.
//     * For quadrilateral 2D3 'volume()' call defaults to 'area()'
//     */
//     KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4Volume, KratosCoreFastSuite)
//     {
//         auto geom = GenerateRightQuadrilateral3D4<NodeType>();
// 
//         KRATOS_CHECK_EXCEPTION_IS_THROWN(geom->Volume(), "Calling base class 'Volume' method instead of derived class one.");
//     }
// 
//     /** Checks the inside test for a given point respect to the quadrilateral
//     * Checks the inside test for a given point respect to the quadrilateral
//     * It performs 4 tests:
//     * A Point inside the quadrilateral: Expected result TRUE
//     * A Point outside the quadrilateral: Expected result FALSE
//     * A Point over a vertex of the quadrilateral: Expected result TRUE
//     * A Point over an edge of the quadrilateral: Expected result TRUE
//     */
//     KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4IsInside, KratosCoreFastSuite)
//     {
//         auto geom = GenerateRightQuadrilateral3D4<NodeType>();
// 
//         Point<3> PointInside(1.0/3.0, 2.0/3.0, 1.0/6.0);
//         Point<3> PointOutside(2.0/3.0, 2.0/3.0, 0.0);
//         Point<3> PointInVertex(0.0, 0.0, 0.0);
//         Point<3> PointInEdge(0.5, 0.5, 0.0);
// 
//         Point<3> LocalCoords;
// 
//         // It appears that the function checks whether the PROJECTION of the point is inside the geometry.
//         KRATOS_CHECK(geom->IsInside(PointInside, LocalCoords, EPSILON));
//         KRATOS_CHECK_IS_FALSE(geom->IsInside(PointOutside, LocalCoords, EPSILON));
//         KRATOS_CHECK(geom->IsInside(PointInVertex, LocalCoords, EPSILON));
//         KRATOS_CHECK(geom->IsInside(PointInEdge, LocalCoords, EPSILON));
//     }
// 
//     /** Tests the Jacobian determinants using 'GI_GAUSS_1' integration method.
//     * Tests the Jacobian determinants using 'GI_GAUSS_1' integration method.
//     */
//     KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4DeterminantOfJacobianArray1, KratosCoreFastSuite)
//     {
//         auto geom = GenerateRightQuadrilateral3D4<NodeType>();
//         const double ExpectedJacobian = 1.0;
// 
//         Vector JacobianDeterminants;
//         geom->DeterminantOfJacobian( JacobianDeterminants, GeometryData::GI_GAUSS_1 );
// 
//         for (unsigned int i=0; i<JacobianDeterminants.size(); ++i)
//         {
//             KRATOS_CHECK_NEAR(JacobianDeterminants[i], ExpectedJacobian, TOLERANCE);
//         }
//     }
// 
//     /** Tests the Jacobian determinants using 'GI_GAUSS_2' integration method.
//     * Tests the Jacobian determinants using 'GI_GAUSS_2' integration method.
//     */
//     KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4DeterminantOfJacobianArray2, KratosCoreFastSuite)
//     {
//         auto geom = GenerateRightQuadrilateral3D4<NodeType>();
//         const double ExpectedJacobian = 1.0;
// 
//         Vector JacobianDeterminants;
//         geom->DeterminantOfJacobian( JacobianDeterminants, GeometryData::GI_GAUSS_2 );
// 
//         for (unsigned int i=0; i<JacobianDeterminants.size(); ++i)
//         {
//             KRATOS_CHECK_NEAR(JacobianDeterminants[i], ExpectedJacobian, TOLERANCE);
//         }
//     }
// 
//     /** Tests the Jacobian determinants using 'GI_GAUSS_3' integration method.
//     * Tests the Jacobian determinants using 'GI_GAUSS_3' integration method.
//     */
//     KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4DeterminantOfJacobianArray3, KratosCoreFastSuite)
//     {
//         auto geom = GenerateRightQuadrilateral3D4<NodeType>();
//         const double ExpectedJacobian = 1.0;
// 
//         Vector JacobianDeterminants;
//         geom->DeterminantOfJacobian( JacobianDeterminants, GeometryData::GI_GAUSS_3 );
// 
//         for (unsigned int i=0; i<JacobianDeterminants.size(); ++i)
//         {
//             KRATOS_CHECK_NEAR(JacobianDeterminants[i], ExpectedJacobian, TOLERANCE);
//         }
//     }
// 
//     /** Tests the Jacobian determinants using 'GI_GAUSS_4' integration method.
//     * Tests the Jacobian determinants using 'GI_GAUSS_4' integration method.
//     */
//     KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4DeterminantOfJacobianArray4, KratosCoreFastSuite)
//     {
//         auto geom = GenerateRightQuadrilateral3D4<NodeType>();
//         const double ExpectedJacobian = 1.0;
// 
//         Vector JacobianDeterminants;
//         geom->DeterminantOfJacobian( JacobianDeterminants, GeometryData::GI_GAUSS_4 );
// 
//         for (unsigned int i=0; i<JacobianDeterminants.size(); ++i)
//         {
//             KRATOS_CHECK_NEAR(JacobianDeterminants[i], ExpectedJacobian, TOLERANCE);
//         }
//     }
// 
//     /** Tests the Jacobian determinants using 'GI_GAUSS_5' integration method.
//     * Tests the Jacobian determinants using 'GI_GAUSS_5' integration method.
//     */
//     KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4DeterminantOfJacobianArray5, KratosCoreFastSuite)
//     {
//         auto geom = GenerateRightQuadrilateral3D4<NodeType>();
//         const double ExpectedJacobian = 1.0;
// 
//         Vector JacobianDeterminants;
//         geom->DeterminantOfJacobian( JacobianDeterminants, GeometryData::GI_GAUSS_5 );
// 
//         for (unsigned int i=0; i<JacobianDeterminants.size(); ++i)
//         {
//             KRATOS_CHECK_NEAR(JacobianDeterminants[i], ExpectedJacobian, TOLERANCE);
//         }
//     }
// 
//     /** Tests the Jacobian determinants using 'GI_GAUSS_1' integration method.
//     * Tests the Jacobian determinants using 'GI_GAUSS_1' integration method.
//     */
//     KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4DeterminantOfJacobianIndex1, KratosCoreFastSuite)
//     {
//         auto geom = GenerateRightQuadrilateral3D4<NodeType>();
//         const double ExpectedJacobian = 1.0;
// 
//         double JacobianDeterminant = geom->DeterminantOfJacobian( 1, GeometryData::GI_GAUSS_1 );
//         KRATOS_CHECK_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
//     }
// 
//     /** Tests the Jacobian determinants using 'GI_GAUSS_2' integration method.
//     * Tests the Jacobian determinants using 'GI_GAUSS_2' integration method.
//     */
//     KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4DeterminantOfJacobianIndex2, KratosCoreFastSuite)
//     {
//         auto geom = GenerateRightQuadrilateral3D4<NodeType>();
//         double JacobianDeterminant = 0.0;
//         const double ExpectedJacobian = 1.0;
// 
//         JacobianDeterminant = geom->DeterminantOfJacobian( 1, GeometryData::GI_GAUSS_2 );
//         KRATOS_CHECK_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
// 
//         JacobianDeterminant = geom->DeterminantOfJacobian( 2, GeometryData::GI_GAUSS_2 );
//         KRATOS_CHECK_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
//     }
// 
//     /** Tests the Jacobian determinants using 'GI_GAUSS_3' integration method.
//     * Tests the Jacobian determinants using 'GI_GAUSS_3' integration method.
//     */
//     KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4DeterminantOfJacobianIndex3, KratosCoreFastSuite)
//     {
//         auto geom = GenerateRightQuadrilateral3D4<NodeType>();
//         double JacobianDeterminant = 0.0;
//         const double ExpectedJacobian = 1.0;
// 
//         JacobianDeterminant = geom->DeterminantOfJacobian( 1, GeometryData::GI_GAUSS_3 );
//         KRATOS_CHECK_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
// 
//         JacobianDeterminant = geom->DeterminantOfJacobian( 2, GeometryData::GI_GAUSS_3 );
//         KRATOS_CHECK_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
// 
//         JacobianDeterminant = geom->DeterminantOfJacobian( 3, GeometryData::GI_GAUSS_3 );
//         KRATOS_CHECK_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
//     }
// 
//     /** Tests the Jacobian determinants using 'GI_GAUSS_4' integration method.
//     * Tests the Jacobian determinants using 'GI_GAUSS_4' integration method.
//     */
//     KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4DeterminantOfJacobianIndex4, KratosCoreFastSuite)
//     {
//         auto geom = GenerateRightQuadrilateral3D4<NodeType>();
//         double JacobianDeterminant = 0.0;
//         const double ExpectedJacobian = 1.0;
// 
//         JacobianDeterminant = geom->DeterminantOfJacobian( 1, GeometryData::GI_GAUSS_4 );
//         KRATOS_CHECK_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
// 
//         JacobianDeterminant = geom->DeterminantOfJacobian( 2, GeometryData::GI_GAUSS_4 );
//         KRATOS_CHECK_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
// 
//         JacobianDeterminant = geom->DeterminantOfJacobian( 3, GeometryData::GI_GAUSS_4 );
//         KRATOS_CHECK_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
// 
//         JacobianDeterminant = geom->DeterminantOfJacobian( 4, GeometryData::GI_GAUSS_4 );
//         KRATOS_CHECK_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
//     }
// 
//     /** Tests the Jacobian determinants using 'GI_GAUSS_4' integration method.
//     * Tests the Jacobian determinants using 'GI_GAUSS_4' integration method.
//     */
//     KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4DeterminantOfJacobianIndex5, KratosCoreFastSuite)
//     {
//         auto geom = GenerateRightQuadrilateral3D4<NodeType>();
//         double JacobianDeterminant = 0.0;
//         const double ExpectedJacobian = 1.0;
// 
//         JacobianDeterminant = geom->DeterminantOfJacobian( 1, GeometryData::GI_GAUSS_5 );
//         KRATOS_CHECK_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
// 
//         JacobianDeterminant = geom->DeterminantOfJacobian( 2, GeometryData::GI_GAUSS_5 );
//         KRATOS_CHECK_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
// 
//         JacobianDeterminant = geom->DeterminantOfJacobian( 3, GeometryData::GI_GAUSS_5 );
//         KRATOS_CHECK_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
// 
//         JacobianDeterminant = geom->DeterminantOfJacobian( 4, GeometryData::GI_GAUSS_5 );
//         KRATOS_CHECK_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
// 
//         JacobianDeterminant = geom->DeterminantOfJacobian( 5, GeometryData::GI_GAUSS_5 );
//         KRATOS_CHECK_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
//     }

} // namespace Testing.
} // namespace Kratos.
