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

// Utility includes
#include "utilities/geometry_utilities.h"

namespace Kratos 
{
namespace Testing 
{
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

    /// Tests

    /** Checks if the number of edges is correct.
    * Checks if the number of edges is correct.
    */
    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4EdgesNumber, KratosCoreGeometriesFastSuite) 
    {
        auto geom = GenerateRightQuadrilateral3D4<Node<3>>();

        KRATOS_CHECK_EQUAL(geom->EdgesNumber(), 4);
    }

    /** Checks if the number of faces is correct.
    * Checks if the number of faces is correct.
    */
    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4FacesNumber, KratosCoreGeometriesFastSuite) 
    {
        auto geom = GenerateRightQuadrilateral3D4<Node<3>>();

        // That for planar geometries it also return the number of edges.
        KRATOS_CHECK_EQUAL(geom->FacesNumber(), 4);
    }

    /** Checks if the area of the quadrilateral is calculated correctly.
    * Checks if the area of the quadrilateral is calculated correctly.
    */
    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4Area, KratosCoreGeometriesFastSuite) 
    {
        auto geom = GenerateRightQuadrilateral3D4<Node<3>>();
        
        KRATOS_CHECK_NEAR(geom->Area(), 1.06947235, TOLERANCE);
//         KRATOS_CHECK_NEAR(geom->Area(), 1.08935, TOLERANCE); // NOTE: Solution from Mathematica
    }

//     /** Checks if the volume of the quadrilateral is calculated correctly.
//     * Checks if the volume of the quadrilateral is calculated correctly.
//     * For quadrilateral 2D3 'volume()' call defaults to 'area()'
//     */
//     KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4Volume, KratosCoreGeometriesFastSuite) 
//     {
//         auto geom = GenerateRightQuadrilateral3D4<Node<3>>();
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
//     KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4IsInside, KratosCoreGeometriesFastSuite) 
//     {
//         auto geom = GenerateRightQuadrilateral3D4<Node<3>>();
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
//     KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4DeterminantOfJacobianArray1, KratosCoreGeometriesFastSuite) 
//     {
//         auto geom = GenerateRightQuadrilateral3D4<Node<3>>();
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
//     KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4DeterminantOfJacobianArray2, KratosCoreGeometriesFastSuite) 
//     {
//         auto geom = GenerateRightQuadrilateral3D4<Node<3>>();
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
//     KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4DeterminantOfJacobianArray3, KratosCoreGeometriesFastSuite) 
//     {
//         auto geom = GenerateRightQuadrilateral3D4<Node<3>>();
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
//     KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4DeterminantOfJacobianArray4, KratosCoreGeometriesFastSuite) 
//     {
//         auto geom = GenerateRightQuadrilateral3D4<Node<3>>();
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
//     KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4DeterminantOfJacobianArray5, KratosCoreGeometriesFastSuite) 
//     {
//         auto geom = GenerateRightQuadrilateral3D4<Node<3>>();
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
//     KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4DeterminantOfJacobianIndex1, KratosCoreGeometriesFastSuite) 
//     {
//         auto geom = GenerateRightQuadrilateral3D4<Node<3>>();
//         const double ExpectedJacobian = 1.0;
// 
//         double JacobianDeterminant = geom->DeterminantOfJacobian( 1, GeometryData::GI_GAUSS_1 );
//         KRATOS_CHECK_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
//     }
// 
//     /** Tests the Jacobian determinants using 'GI_GAUSS_2' integration method.
//     * Tests the Jacobian determinants using 'GI_GAUSS_2' integration method.
//     */
//     KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4DeterminantOfJacobianIndex2, KratosCoreGeometriesFastSuite) 
//     {
//         auto geom = GenerateRightQuadrilateral3D4<Node<3>>();
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
//     KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4DeterminantOfJacobianIndex3, KratosCoreGeometriesFastSuite) 
//     {
//         auto geom = GenerateRightQuadrilateral3D4<Node<3>>();
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
//     KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4DeterminantOfJacobianIndex4, KratosCoreGeometriesFastSuite) 
//     {
//         auto geom = GenerateRightQuadrilateral3D4<Node<3>>();
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
//     KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4DeterminantOfJacobianIndex5, KratosCoreGeometriesFastSuite)
//     {
//         auto geom = GenerateRightQuadrilateral3D4<Node<3>>();
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
