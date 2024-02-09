//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Mohit Tyagi
//

// System includes
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_2d_9.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"
#include "geometries/prism_3d_6.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_27.h"
#include "tests/cpp_tests/geometries/test_geometry.h"
#include "utilities/element_size_calculator.h"

namespace Kratos
{
    namespace Testing
    {
        namespace
        {
            template<class TPrimalFunction, class TDerivativeFunction>
            void RunElementSizeCalculatorDerivativesTest(
                Geometry<NodeType>& rGeometry,
                const TPrimalFunction& rPrimalMethod,
                const TDerivativeFunction& rDerivativeMethod,
                const double Delta,
                const double Tolerance)
            {
                const double fd_ref_value = rPrimalMethod(rGeometry);

                for (unsigned int c = 0; c < rGeometry.PointsNumber(); ++c) {
                    for (unsigned int k = 0; k < rGeometry.WorkingSpaceDimension(); ++k) {
                        const double analytical_derivative = rDerivativeMethod(c, k, rGeometry);

                        rGeometry[c].Coordinates()[k] += Delta;
                        const double fd_value = rPrimalMethod(rGeometry);
                        rGeometry[c].Coordinates()[k] -= Delta;
                        const double fd_derivative = (fd_value - fd_ref_value) / Delta;

                        KRATOS_EXPECT_NEAR(analytical_derivative, fd_derivative, Tolerance);
                    }
                }
            }
        }

        KRATOS_TEST_CASE_IN_SUITE(Prism3D6ElementSizeCase1, KratosCoreFastSuite)
        {
            Geometry<NodeType>::PointsArrayType nodes;
            nodes.push_back(NodeType::Pointer(new NodeType(1, 0.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(2, 1.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(3, 0.0, 1.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(4, 0.0, 0.0, 1.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(5, 1.0, 0.0, 1.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(6, 0.0, 1.0, 1.0)));
            const auto prism1 = *GeometryType::Pointer(new Prism3D6<NodeType>(nodes));
            const double minimum_size = ElementSizeCalculator<3, 6>::MinimumElementSize(prism1);
            const double average_size = ElementSizeCalculator<3, 6>::AverageElementSize(prism1);
            KRATOS_EXPECT_NEAR(minimum_size, std::sqrt(0.5), TOLERANCE);
            KRATOS_EXPECT_NEAR(average_size, std::pow(0.5, 1.0 / 3.0), TOLERANCE);
        }

        KRATOS_TEST_CASE_IN_SUITE(Prism3D6ElementSizeCase2, KratosCoreFastSuite)
        {
            Geometry<NodeType>::PointsArrayType nodes;
            nodes.push_back(NodeType::Pointer(new NodeType(1, 0.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(2, 1.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(3, 0.0, 1.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(4, 0.0, 0.0, 0.5)));
            nodes.push_back(NodeType::Pointer(new NodeType(5, 1.0, 0.0, 0.5)));
            nodes.push_back(NodeType::Pointer(new NodeType(6, 0.0, 1.0, 0.5)));
            const auto prism2 = *GeometryType::Pointer(new Prism3D6<NodeType>(nodes));
            const double minimum_size = ElementSizeCalculator<3, 6>::MinimumElementSize(prism2);
            const double average_size = ElementSizeCalculator<3, 6>::AverageElementSize(prism2);
            KRATOS_EXPECT_NEAR(minimum_size, 0.5, TOLERANCE);
            KRATOS_EXPECT_NEAR(average_size, std::pow(0.25, 1.0 / 3.0), TOLERANCE); //
        }

        KRATOS_TEST_CASE_IN_SUITE(Triangle2D3MinimumElementSizeDerivatives, KratosCoreFastSuite)
        {
            Geometry<NodeType>::PointsArrayType nodes;
            nodes.push_back(NodeType::Pointer(new NodeType(1, 0.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(2, 1.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(3, 0.0, 1.0, 0.0)));
            auto p_geometry = GeometryType::Pointer(new Triangle2D3<NodeType>(nodes));
            RunElementSizeCalculatorDerivativesTest(
                *p_geometry, ElementSizeCalculator<2, 3>::MinimumElementSize,
                ElementSizeCalculator<2, 3>::MinimumElementSizeDerivative, 1e-8, 1e-7);
        }

        KRATOS_TEST_CASE_IN_SUITE(Triangle2D3AverageElementSizeDerivatives, KratosCoreFastSuite)
        {
            Geometry<NodeType>::PointsArrayType nodes;
            nodes.push_back(NodeType::Pointer(new NodeType(1, 0.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(2, 1.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(3, 0.0, 1.0, 0.0)));
            auto p_geometry = GeometryType::Pointer(new Triangle2D3<NodeType>(nodes));
            RunElementSizeCalculatorDerivativesTest(
                *p_geometry, ElementSizeCalculator<2, 3>::AverageElementSize,
                ElementSizeCalculator<2, 3>::AverageElementSizeDerivative, 1e-8, 1e-7);
        }

        KRATOS_TEST_CASE_IN_SUITE(Triangle2D6MinimumElementSizeDerivatives, KratosCoreFastSuite)
        {
            Geometry<NodeType>::PointsArrayType nodes;
            nodes.push_back(NodeType::Pointer(new NodeType(1, 0.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(2, 1.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(3, 0.0, 1.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(4, 0.5, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(5, 0.5, 0.5, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(6, 0.0, 0.5, 0.0)));
            auto p_geometry = GeometryType::Pointer(new Triangle2D6<NodeType>(nodes));
            RunElementSizeCalculatorDerivativesTest(
                *p_geometry, ElementSizeCalculator<2, 6>::MinimumElementSize,
                ElementSizeCalculator<2, 6>::MinimumElementSizeDerivative, 1e-8, 1e-7);
        }

        KRATOS_TEST_CASE_IN_SUITE(Triangle2D6AverageElementSizeDerivatives, KratosCoreFastSuite)
        {
            Geometry<NodeType>::PointsArrayType nodes;
            nodes.push_back(NodeType::Pointer(new NodeType(1, 0.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(2, 1.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(3, 0.0, 1.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(4, 0.5, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(5, 0.5, 0.5, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(6, 0.0, 0.5, 0.0)));
            auto p_geometry = GeometryType::Pointer(new Triangle2D6<NodeType>(nodes));
            RunElementSizeCalculatorDerivativesTest(
                *p_geometry, ElementSizeCalculator<2, 6>::AverageElementSize,
                ElementSizeCalculator<2, 6>::AverageElementSizeDerivative, 1e-8, 1e-7);
        }

        KRATOS_TEST_CASE_IN_SUITE(Quadrilateral2D4MinimumElementSizeDerivatives, KratosCoreFastSuite)
        {
            Geometry<NodeType>::PointsArrayType nodes;
            nodes.push_back(NodeType::Pointer(new NodeType(1, 0.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(2, 1.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(3, 1.0, 1.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(4, 0.0, 1.0, 0.0)));
            auto p_geometry = GeometryType::Pointer(new Quadrilateral2D4<NodeType>(nodes));
            RunElementSizeCalculatorDerivativesTest(
                *p_geometry, ElementSizeCalculator<2, 4>::MinimumElementSize,
                ElementSizeCalculator<2, 4>::MinimumElementSizeDerivative, 1e-8, 1e-7);
        }

        KRATOS_TEST_CASE_IN_SUITE(Quadrilateral2D4AverageElementSizeDerivatives, KratosCoreFastSuite)
        {
            Geometry<NodeType>::PointsArrayType nodes;
            nodes.push_back(NodeType::Pointer(new NodeType(1, 0.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(2, 1.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(3, 1.0, 1.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(4, 0.0, 1.0, 0.0)));
            auto p_geometry = GeometryType::Pointer(new Quadrilateral2D4<NodeType>(nodes));
            RunElementSizeCalculatorDerivativesTest(
                *p_geometry, ElementSizeCalculator<2, 4>::AverageElementSize,
                ElementSizeCalculator<2, 4>::AverageElementSizeDerivative, 1e-8, 1e-7);
        }

        KRATOS_TEST_CASE_IN_SUITE(Quadrilateral2D9MinimumElementSizeDerivatives, KratosCoreFastSuite)
        {
            Geometry<NodeType>::PointsArrayType nodes;
            nodes.push_back(NodeType::Pointer(new NodeType(1, 0.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(2, 1.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(3, 1.0, 1.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(4, 0.0, 1.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(5, 0.5, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(6, 1.0, 0.5, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(7, 0.5, 1.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(8, 0.0, 0.5, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(9, 0.5, 0.5, 0.0)));
            auto p_geometry = GeometryType::Pointer(new Quadrilateral2D9<NodeType>(nodes));
            RunElementSizeCalculatorDerivativesTest(
                *p_geometry, ElementSizeCalculator<2, 9>::MinimumElementSize,
                ElementSizeCalculator<2, 9>::MinimumElementSizeDerivative, 1e-8, 1e-7);
        }

        KRATOS_TEST_CASE_IN_SUITE(Quadrilateral2D9AverageElementSizeDerivatives, KratosCoreFastSuite)
        {
            Geometry<NodeType>::PointsArrayType nodes;
            nodes.push_back(NodeType::Pointer(new NodeType(1, 0.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(2, 1.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(3, 1.0, 1.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(4, 0.0, 1.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(5, 0.5, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(6, 1.0, 0.5, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(7, 0.5, 1.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(8, 0.0, 0.5, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(9, 0.5, 0.5, 0.0)));
            auto p_geometry = GeometryType::Pointer(new Quadrilateral2D9<NodeType>(nodes));
            RunElementSizeCalculatorDerivativesTest(
                *p_geometry, ElementSizeCalculator<2, 9>::AverageElementSize,
                ElementSizeCalculator<2, 9>::AverageElementSizeDerivative, 1e-8, 1e-7);
        }

        KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4MinimumElementSizeDerivatives, KratosCoreFastSuite)
        {
            Geometry<NodeType>::PointsArrayType nodes;
            nodes.push_back(NodeType::Pointer(new NodeType(1, 0.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(2, 1.0, 2.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(3, 0.5, 0.5, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(4, 0.2, 0.3, 1.0)));
            auto p_geometry = GeometryType::Pointer(new Tetrahedra3D4<NodeType>(nodes));
            RunElementSizeCalculatorDerivativesTest(
                *p_geometry, ElementSizeCalculator<3, 4>::MinimumElementSize,
                ElementSizeCalculator<3, 4>::MinimumElementSizeDerivative, 1e-8, 1e-7);
        }

        KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4AverageElementSizeDerivatives, KratosCoreFastSuite)
        {
            Geometry<NodeType>::PointsArrayType nodes;
            nodes.push_back(NodeType::Pointer(new NodeType(1, 0.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(2, 1.0, 2.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(3, 0.2, 0.3, 1.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(4, 0.5, 0.5, 0.0)));
            auto p_geometry = GeometryType::Pointer(new Tetrahedra3D4<NodeType>(nodes));
            RunElementSizeCalculatorDerivativesTest(
                *p_geometry, ElementSizeCalculator<3, 4>::AverageElementSize,
                ElementSizeCalculator<3, 4>::AverageElementSizeDerivative, 1e-8, 1e-7);
        }

        KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D10MinimumElementSizeDerivatives, KratosCoreFastSuite)
        {
            Geometry<NodeType>::PointsArrayType nodes;
            nodes.push_back(NodeType::Pointer(new NodeType(1, 0.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(2, 0.5, 0.5, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(3, 1.0, 2.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(4, 0.2, 0.3, 1.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(5, 0.25, 0.25, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(6, 0.75, 1.25, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(7, 0.5, 1.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(8, 0.1, 0.15, 0.5)));
            nodes.push_back(NodeType::Pointer(new NodeType(9, 0.35, 0.4, 0.5)));
            nodes.push_back(NodeType::Pointer(new NodeType(10, 0.6, 1.15, 0.5)));
            auto p_geometry = GeometryType::Pointer(new Tetrahedra3D10<NodeType>(nodes));
            RunElementSizeCalculatorDerivativesTest(
                *p_geometry, ElementSizeCalculator<3, 10>::MinimumElementSize,
                ElementSizeCalculator<3, 10>::MinimumElementSizeDerivative, 1e-8, 1e-7);
        }

        KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D10AverageElementSizeDerivatives, KratosCoreFastSuite)
        {
            Geometry<NodeType>::PointsArrayType nodes;
            nodes.push_back(NodeType::Pointer(new NodeType(1, 0.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(2, 0.5, 0.5, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(3, 1.0, 2.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(4, 0.2, 0.3, 1.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(5, 0.25, 0.25, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(6, 0.75, 1.25, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(7, 0.5, 1.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(8, 0.1, 0.15, 0.5)));
            nodes.push_back(NodeType::Pointer(new NodeType(9, 0.35, 0.4, 0.5)));
            nodes.push_back(NodeType::Pointer(new NodeType(10, 0.6, 1.15, 0.5)));
            auto p_geometry = GeometryType::Pointer(new Tetrahedra3D10<NodeType>(nodes));
            RunElementSizeCalculatorDerivativesTest(
                *p_geometry, ElementSizeCalculator<3, 10>::AverageElementSize,
                ElementSizeCalculator<3, 10>::AverageElementSizeDerivative, 1e-8, 1e-7);
        }

        KRATOS_TEST_CASE_IN_SUITE(Prism3D6MinimumElementSizeDerivatives, KratosCoreFastSuite)
        {
            Geometry<NodeType>::PointsArrayType nodes;
            nodes.push_back(NodeType::Pointer(new NodeType(1, 0.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(2, 1.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(3, 0.0, 1.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(4, 0.0, 0.0, 1.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(5, 1.0, 0.0, 1.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(6, 0.0, 1.0, 1.0)));
            auto p_geometry = GeometryType::Pointer(new Prism3D6<NodeType>(nodes));
            RunElementSizeCalculatorDerivativesTest(
                *p_geometry, ElementSizeCalculator<3, 6>::MinimumElementSize,
                ElementSizeCalculator<3, 6>::MinimumElementSizeDerivative, 1e-8, 1e-7);
        }

        KRATOS_TEST_CASE_IN_SUITE(Prism3D6AverageElementSizeDerivatives, KratosCoreFastSuite)
        {
            Geometry<NodeType>::PointsArrayType nodes;
            nodes.push_back(NodeType::Pointer(new NodeType(1, 0.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(2, 1.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(3, 0.0, 1.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(4, 0.0, 0.0, 1.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(5, 1.0, 0.0, 1.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(6, 0.0, 1.0, 1.0)));
            auto p_geometry = GeometryType::Pointer(new Prism3D6<NodeType>(nodes));
            RunElementSizeCalculatorDerivativesTest(
                *p_geometry, ElementSizeCalculator<3, 6>::AverageElementSize,
                ElementSizeCalculator<3, 6>::AverageElementSizeDerivative, 1e-8, 1e-7);
        }

        KRATOS_TEST_CASE_IN_SUITE(Hexahedra3D8MinimumElementSizeDerivatives, KratosCoreFastSuite)
        {
            Geometry<NodeType>::PointsArrayType nodes;
            nodes.push_back(NodeType::Pointer(new NodeType(1, 0.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(2, 1.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(3, 1.0, 1.1, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(4, 0.0, 1.1, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(5, 0.1, 0.1, 0.5)));
            nodes.push_back(NodeType::Pointer(new NodeType(6, 0.7, 0.1, 0.4)));
            nodes.push_back(NodeType::Pointer(new NodeType(7, 0.7, 0.8, 0.3)));
            nodes.push_back(NodeType::Pointer(new NodeType(8, 0.1, 0.8, 0.6)));
            auto p_geometry = GeometryType::Pointer(new Hexahedra3D8<NodeType>(nodes));
            RunElementSizeCalculatorDerivativesTest(
                *p_geometry, ElementSizeCalculator<3, 8>::MinimumElementSize,
                ElementSizeCalculator<3, 8>::MinimumElementSizeDerivative, 1e-8, 1e-7);
        }

        KRATOS_TEST_CASE_IN_SUITE(Hexahedra3D8AverageElementSizeDerivatives, KratosCoreFastSuite)
        {
            Geometry<NodeType>::PointsArrayType nodes;
            nodes.push_back(NodeType::Pointer(new NodeType(1, 0.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(2, 1.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(3, 1.0, 1.1, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(4, 0.0, 1.1, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(5, 0.1, 0.1, 0.5)));
            nodes.push_back(NodeType::Pointer(new NodeType(6, 0.7, 0.1, 0.4)));
            nodes.push_back(NodeType::Pointer(new NodeType(7, 0.7, 0.8, 0.3)));
            nodes.push_back(NodeType::Pointer(new NodeType(8, 0.1, 0.8, 0.6)));
            auto p_geometry = GeometryType::Pointer(new Hexahedra3D8<NodeType>(nodes));
            RunElementSizeCalculatorDerivativesTest(
                *p_geometry, ElementSizeCalculator<3, 8>::AverageElementSize,
                ElementSizeCalculator<3, 8>::AverageElementSizeDerivative, 1e-8, 1e-7);
        }

        KRATOS_TEST_CASE_IN_SUITE(Hexahedra3D27MinimumElementSizeDerivatives, KratosCoreFastSuite)
        {
            Geometry<NodeType>::PointsArrayType nodes;
            nodes.push_back(NodeType::Pointer(new NodeType(1, 0.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(2, 1.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(3, 1.0, 1.1, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(4, 0.0, 1.1, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(5, 0.1, 0.1, 0.5)));
            nodes.push_back(NodeType::Pointer(new NodeType(6, 0.7, 0.1, 0.4)));
            nodes.push_back(NodeType::Pointer(new NodeType(7, 0.7, 0.8, 0.3)));
            nodes.push_back(NodeType::Pointer(new NodeType(8, 0.1, 0.8, 0.6)));
            nodes.push_back(NodeType::Pointer(new NodeType(9, 0.5, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(10, 1.0, 0.55, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(11, 0.5, 1.1, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(12, 0.0, 0.55, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(13, 0.05, 0.05, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(14, 0.85, 0.05, 0.2)));
            nodes.push_back(NodeType::Pointer(new NodeType(15, 0.85, 0.95, 0.15)));
            nodes.push_back(NodeType::Pointer(new NodeType(16, 0.05, 0.95, 0.3)));
            nodes.push_back(NodeType::Pointer(new NodeType(17, 0.4, 0.1, 0.45)));
            nodes.push_back(NodeType::Pointer(new NodeType(18, 0.7, 0.35, 0.35)));
            nodes.push_back(NodeType::Pointer(new NodeType(19, 0.4, 0.8, 0.45)));
            nodes.push_back(NodeType::Pointer(new NodeType(20, 0.7, 0.45, 0.55)));
            nodes.push_back(NodeType::Pointer(new NodeType(21, 0.4, 0.05, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(22, 0.45, 0.05, 0.1)));
            nodes.push_back(NodeType::Pointer(new NodeType(23, 0.85, 0.325, 0.175)));
            nodes.push_back(NodeType::Pointer(new NodeType(24, 0.6, 0.95, 0.225)));
            nodes.push_back(NodeType::Pointer(new NodeType(25, 0.05, 0.55, 0.15)));
            nodes.push_back(NodeType::Pointer(new NodeType(26, 0.7, 0.4, 0.45)));
            nodes.push_back(NodeType::Pointer(new NodeType(27, 0.45, 0.5, 0.1625)));
            auto p_geometry = GeometryType::Pointer(new Hexahedra3D27<NodeType>(nodes));
            RunElementSizeCalculatorDerivativesTest(
                *p_geometry, ElementSizeCalculator<3, 27>::MinimumElementSize,
                ElementSizeCalculator<3, 27>::MinimumElementSizeDerivative, 1e-8, 1e-7);
        }

        KRATOS_TEST_CASE_IN_SUITE(Hexahedra3D27AverageElementSizeDerivatives, KratosCoreFastSuite)
        {
            Geometry<NodeType>::PointsArrayType nodes;
            nodes.push_back(NodeType::Pointer(new NodeType(1, 0.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(2, 1.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(3, 1.0, 1.1, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(4, 0.0, 1.1, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(5, 0.1, 0.1, 0.5)));
            nodes.push_back(NodeType::Pointer(new NodeType(6, 0.7, 0.1, 0.4)));
            nodes.push_back(NodeType::Pointer(new NodeType(7, 0.7, 0.8, 0.3)));
            nodes.push_back(NodeType::Pointer(new NodeType(8, 0.1, 0.8, 0.6)));
            nodes.push_back(NodeType::Pointer(new NodeType(9, 0.5, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(10, 1.0, 0.55, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(11, 0.5, 1.1, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(12, 0.0, 0.55, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(13, 0.05, 0.05, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(14, 0.85, 0.05, 0.2)));
            nodes.push_back(NodeType::Pointer(new NodeType(15, 0.85, 0.95, 0.15)));
            nodes.push_back(NodeType::Pointer(new NodeType(16, 0.05, 0.95, 0.3)));
            nodes.push_back(NodeType::Pointer(new NodeType(17, 0.4, 0.1, 0.45)));
            nodes.push_back(NodeType::Pointer(new NodeType(18, 0.7, 0.35, 0.35)));
            nodes.push_back(NodeType::Pointer(new NodeType(19, 0.4, 0.8, 0.45)));
            nodes.push_back(NodeType::Pointer(new NodeType(20, 0.7, 0.45, 0.55)));
            nodes.push_back(NodeType::Pointer(new NodeType(21, 0.4, 0.05, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(22, 0.45, 0.05, 0.1)));
            nodes.push_back(NodeType::Pointer(new NodeType(23, 0.85, 0.325, 0.175)));
            nodes.push_back(NodeType::Pointer(new NodeType(24, 0.6, 0.95, 0.225)));
            nodes.push_back(NodeType::Pointer(new NodeType(25, 0.05, 0.55, 0.15)));
            nodes.push_back(NodeType::Pointer(new NodeType(26, 0.7, 0.4, 0.45)));
            nodes.push_back(NodeType::Pointer(new NodeType(27, 0.45, 0.5, 0.1625)));
            auto p_geometry = GeometryType::Pointer(new Hexahedra3D27<NodeType>(nodes));
            RunElementSizeCalculatorDerivativesTest(
                *p_geometry, ElementSizeCalculator<3, 27>::AverageElementSize,
                ElementSizeCalculator<3, 27>::AverageElementSizeDerivative, 1e-8, 1e-7);
        }

    } // namespace Testing

} // namespace Kratos
