//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes
#include <set>

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/prism_3d_15.h"
#include "tests/cpp_tests/geometries/test_geometry.h"

namespace Kratos::Testing
{
using PointType = Node;
using PointPtrType = Node::Pointer;
using Prism15GeometryType = Prism3D15<PointType>;
using Prism15GeometryPtrType = Prism3D15<PointType>::Pointer;

/** Generates a sample Prism3D15.
 * Generates a trirectangular prism on the origin with positive volume and side 1.
 * @return  Pointer to a Prism3D15
 */
Prism15GeometryPtrType GenerateRegularPrism3D15() {
    return Prism15GeometryPtrType(new Prism15GeometryType(
    GeneratePoint<PointType>(0.0, 0.0, 0.0), //node0 //vertex of bottom
    GeneratePoint<PointType>(1.0, 0.0, 0.0), //node1
    GeneratePoint<PointType>(0.0, 1.0, 0.0), //node2
    GeneratePoint<PointType>(0.0, 0.0, 1.0), //node3 //vertex of top
    GeneratePoint<PointType>(1.0, 0.0, 1.0), //node4
    GeneratePoint<PointType>(0.0, 1.0, 1.0), //node5
    GeneratePoint<PointType>(0.5, 0.0, 0.0), //node6 //mid of bottom
    GeneratePoint<PointType>(0.5, 0.5, 0.0), //node7
    GeneratePoint<PointType>(0.0, 0.5, 0.0), //node8
    GeneratePoint<PointType>(0.0, 0.0, 0.5), //node9 //vertex of mid height
    GeneratePoint<PointType>(1.0, 0.0, 0.5), //node10
    GeneratePoint<PointType>(0.0, 1.0, 0.5), //node11
    GeneratePoint<PointType>(0.5, 0.0, 1.0), //node12 //mid of top
    GeneratePoint<PointType>(0.5, 0.5, 1.0), //node13
    GeneratePoint<PointType>(0.0, 0.5, 1.0)  //node14
    ));
}

/** Checks if the number of edges is correct.
 * Checks if the number of edges is correct.
 */
KRATOS_TEST_CASE_IN_SUITE(Prism3D15EdgesNumber, KratosCoreGeometriesFastSuite) {
    auto geomRegular = GenerateRegularPrism3D15();

    KRATOS_CHECK_EQUAL(geomRegular->EdgesNumber(), 9);
}

/** Checks if the number of faces is correct.
 * Checks if the number of faces is correct.
 */
KRATOS_TEST_CASE_IN_SUITE(Prism3D15FacesNumber, KratosCoreGeometriesFastSuite) {
    auto geomRegular = GenerateRegularPrism3D15();

    KRATOS_CHECK_EQUAL(geomRegular->FacesNumber(), 5);
}

/** Checks if the characteristic length of the prism is calculated correctly.
 * Checks if the characteristic length of the prism is calculated correctly.
 */
KRATOS_TEST_CASE_IN_SUITE(Prism3D15Length, KratosCoreGeometriesFastSuite) {
    auto geomRegular = GenerateRegularPrism3D15();

    KRATOS_CHECK_NEAR(geomRegular->Length(), 0.264567, TOLERANCE);
}

/** Checks if the area of the prism is calculated correctly.
 * Checks if the area of the prism is calculated correctly.
 */
KRATOS_TEST_CASE_IN_SUITE(Prism3D15Area, KratosCoreGeometriesFastSuite) {
    auto geomRegular = GenerateRegularPrism3D15();

    KRATOS_CHECK_NEAR(geomRegular->Area(),  0.5, TOLERANCE);
}

/** Checks if the volume of the prism is calculated correctly.
 * Checks if the volume of the prism is calculated correctly.
 * For prism 3D6 'volume()' call defaults to 'area()'
 */
KRATOS_TEST_CASE_IN_SUITE(Prism3D15Volume, KratosCoreGeometriesFastSuite) {
    auto geomRegular = GenerateRegularPrism3D15();
    KRATOS_CHECK_NEAR(geomRegular->Volume(),  0.5, TOLERANCE);
}

/** Checks the inside test for a given point respect to the prism
* Checks the inside test for a given point respect to the prism
* It performs 4 tests:
* A Point inside the prism: Expected result TRUE
* A Point outside the prism: Expected result FALSE
* A Point over a vertex of the prism: Expected result TRUE
* A Point over an edge of the prism: Expected result TRUE
*/
KRATOS_TEST_CASE_IN_SUITE(Prism3D15IsInside, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateRegularPrism3D15();

    Point PointInside(0.1666, 0.1666, 0.1666);
    Point PointOutside(0.66, 0.66, 0.66);
    Point PointInVertex(0.0, 0.0, 0.0);
    Point PointInEdge(0.33, 0.33, 0.33);

    Point LocalCoords;

    KRATOS_CHECK(geom->IsInside(PointInside, LocalCoords, EPSILON));
    KRATOS_CHECK_IS_FALSE(geom->IsInside(PointOutside, LocalCoords, EPSILON));
    KRATOS_CHECK(geom->IsInside(PointInVertex, LocalCoords, EPSILON));
    KRATOS_CHECK(geom->IsInside(PointInEdge, LocalCoords, EPSILON));
}

/** Checks the point local coordinates for a given point respect to the
* prism. The centre of the prism is selected due to its known
* solution.
*/
KRATOS_TEST_CASE_IN_SUITE(Prism3D15PointLocalCoordinates, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateRegularPrism3D15();

    // Compute the global coordinates of the centre
    auto points = geom->Points();
    auto centre = Point{points[0] + points[1] + points[2] + points[3] + points[4] + points[5]};
    centre /= 6.0;

    // Compute the centre local coordinates
    array_1d<double, 3> centre_local_coords;
    geom->PointLocalCoordinates(centre_local_coords, centre);

    KRATOS_CHECK_NEAR(centre_local_coords(0), 1.0/3.0, TOLERANCE);
    KRATOS_CHECK_NEAR(centre_local_coords(1), 1.0/3.0, TOLERANCE);
    KRATOS_CHECK_NEAR(centre_local_coords(2), 1.0/2.0, TOLERANCE);
}

/** Tests the area using 'GI_GAUSS_1' integration method.
* Tests the area using 'GI_GAUSS_1' integration method.
*/
KRATOS_TEST_CASE_IN_SUITE(Prism3D15GaussPoint1, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateRegularPrism3D15();

    const double expected_vol = 0.5;

    KRATOS_CHECK_NEAR(CalculateAreaByIntegration(*geom, GeometryData::IntegrationMethod::GI_GAUSS_1), expected_vol, TOLERANCE);
    VerifyStrainExactness(*geom, GeometryData::IntegrationMethod::GI_GAUSS_1);
}

/** Tests the area using 'GI_GAUSS_2' integration method.
* Tests the area using 'GI_GAUSS_2' integration method.
*/
KRATOS_TEST_CASE_IN_SUITE(Prism3D15GaussPoint2, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateRegularPrism3D15();

    const double expected_vol = 0.5;

    KRATOS_CHECK_NEAR(CalculateAreaByIntegration(*geom, GeometryData::IntegrationMethod::GI_GAUSS_2), expected_vol, TOLERANCE);
    VerifyStrainExactness(*geom, GeometryData::IntegrationMethod::GI_GAUSS_2);
}

/** Tests the area using 'GI_GAUSS_3' integration method.
* Tests the area using 'GI_GAUSS_3' integration method.
*/
KRATOS_TEST_CASE_IN_SUITE(Prism3D15GaussPoint3, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateRegularPrism3D15();

    const double expected_vol = 0.5;

    KRATOS_CHECK_NEAR(CalculateAreaByIntegration(*geom, GeometryData::IntegrationMethod::GI_GAUSS_3), expected_vol, TOLERANCE);
    VerifyStrainExactness(*geom, GeometryData::IntegrationMethod::GI_GAUSS_3);
}

/** Tests the area using 'GI_GAUSS_4' integration method.
* Tests the area using 'GI_GAUSS_4' integration method.
*/
KRATOS_TEST_CASE_IN_SUITE(Prism3D15GaussPoint4, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateRegularPrism3D15();

    const double expected_vol = 0.5;

    KRATOS_CHECK_NEAR(CalculateAreaByIntegration(*geom, GeometryData::IntegrationMethod::GI_GAUSS_4), expected_vol, TOLERANCE);
    VerifyStrainExactness(*geom, GeometryData::IntegrationMethod::GI_GAUSS_4);
}

/** Tests the area using 'GI_GAUSS_5' integration method.
* Tests the area using 'GI_GAUSS_5' integration method.
*/
KRATOS_TEST_CASE_IN_SUITE(Prism3D15GaussPoint5, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateRegularPrism3D15();

    const double expected_vol = 0.5;

    KRATOS_CHECK_NEAR(CalculateAreaByIntegration(*geom, GeometryData::IntegrationMethod::GI_GAUSS_5), expected_vol, TOLERANCE);
    VerifyStrainExactness(*geom, GeometryData::IntegrationMethod::GI_GAUSS_5);
}

void ShapeFunctionSanityCheck(
    Geometry<Node>::Pointer geom,
    const GeometryData::IntegrationMethod method,
    bool isFullIntegration)
{

    const unsigned int numPoints = geom->PointsNumber();
    const unsigned int numGauss = geom->IntegrationPointsNumber(method);

    Vector detJ;
    GeometryData::ShapeFunctionsGradientsType DN_DX;
    geom->ShapeFunctionsIntegrationPointsGradients(DN_DX, detJ, method);

    Matrix shapeFunctions = ZeroMatrix(numGauss, numPoints);
    noalias(shapeFunctions) = geom->ShapeFunctionsValues(method);

    // Shape functions are a partition of unity
    for (unsigned int g = 0; g < numGauss; g++) {
        double shapeFuncSum = 0.0;
        for (unsigned int i = 0; i < numPoints; i++) {
            shapeFuncSum += shapeFunctions(g, i);
        }
        KRATOS_CHECK_NEAR(shapeFuncSum, 1.0, TOLERANCE);
    }

    // Shape function gradients shoud sum to 0
    // Gradient of interpolated nodal coordinates should be dx_i/d_xj = kronecker(1, j)
    for (unsigned int g = 0; g < numGauss; g++) {
        double gradSumX = 0.0;
        double gradSumY = 0.0;
        double gradSumZ = 0.0;

        Matrix coordinateGradients = ZeroMatrix(3,3);
        const auto& gradients = DN_DX[g];
        for (unsigned int n = 0; n < numPoints; n++) {
            gradSumX += gradients(n, 0);
            gradSumY += gradients(n, 1);
            gradSumZ += gradients(n, 2);

            const auto& node = (*geom)[n].Coordinates();

            for (unsigned int i = 0; i < 3; i++) {
                for (unsigned int j = 0; j < 3; j++) {
                    coordinateGradients(i, j) += gradients(n, j) * node[i];
                }
            }
        }
        KRATOS_CHECK_NEAR(gradSumX, 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(gradSumY, 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(gradSumZ, 0.0, TOLERANCE);
        for (unsigned int i = 0; i < 3; i++) {
            for (unsigned int j = 0; j < 3; j++) {
                KRATOS_CHECK_NEAR(coordinateGradients(i, j), (i == j ? 1.0 : 0.0), TOLERANCE);
            }
        }
    }

    // Sum of weights times detJ is volume (only if fully integrated)
    if (isFullIntegration) {
        const auto& integrationPoints = geom->IntegrationPoints(method);

        double volume = 0.0;
        double totalWeight = 0.0;
        for (unsigned int g = 0; g < numGauss; g++) {
            volume += detJ[g] * integrationPoints[g].Weight();
            totalWeight += integrationPoints[g].Weight();
        }

        KRATOS_CHECK_NEAR(volume, 0.5, TOLERANCE);
        KRATOS_CHECK_NEAR(totalWeight, 0.5, TOLERANCE);
    }

}

KRATOS_TEST_CASE_IN_SUITE(Prism3D15ShapeFunctionsSanityCheckGauss1, KratosCoreGeometriesFastSuite) {
    ShapeFunctionSanityCheck(GenerateRegularPrism3D15(), GeometryData::IntegrationMethod::GI_GAUSS_1, false);
}

KRATOS_TEST_CASE_IN_SUITE(Prism3D15ShapeFunctionsSanityCheckGauss2, KratosCoreGeometriesFastSuite) {
    ShapeFunctionSanityCheck(GenerateRegularPrism3D15(), GeometryData::IntegrationMethod::GI_GAUSS_2, false);
}

KRATOS_TEST_CASE_IN_SUITE(Prism3D15ShapeFunctionsSanityCheckGauss3, KratosCoreGeometriesFastSuite) {
    ShapeFunctionSanityCheck(GenerateRegularPrism3D15(), GeometryData::IntegrationMethod::GI_GAUSS_3, true);
}

KRATOS_TEST_CASE_IN_SUITE(Prism3D15ShapeFunctionsSanityCheckGauss4, KratosCoreGeometriesFastSuite) {
    ShapeFunctionSanityCheck(GenerateRegularPrism3D15(), GeometryData::IntegrationMethod::GI_GAUSS_4, true);
}

KRATOS_TEST_CASE_IN_SUITE(Prism3D15ShapeFunctionsSanityCheckGauss5, KratosCoreGeometriesFastSuite) {
    ShapeFunctionSanityCheck(GenerateRegularPrism3D15(), GeometryData::IntegrationMethod::GI_GAUSS_5, true);
}

}  // namespace Kratos::Testing.
