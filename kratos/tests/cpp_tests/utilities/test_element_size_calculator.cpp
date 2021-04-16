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
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/prism_3d_6.h"
#include "geometries/hexahedra_3d_8.h"
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

                        KRATOS_CHECK_NEAR(analytical_derivative, fd_derivative, Tolerance);
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
            KRATOS_CHECK_NEAR(minimum_size, std::sqrt(0.5), TOLERANCE);
            KRATOS_CHECK_NEAR(average_size, std::pow(0.5, 1.0 / 3.0), TOLERANCE);
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
            KRATOS_CHECK_NEAR(minimum_size, 0.5, TOLERANCE);
            KRATOS_CHECK_NEAR(average_size, std::pow(0.25, 1.0 / 3.0), TOLERANCE); //
        }

        KRATOS_TEST_CASE_IN_SUITE(Triangle2D3MinimumElementSizeDerivatives, KratosCoreFastSuite)
        {
            Geometry<NodeType>::PointsArrayType nodes;
            nodes.push_back(NodeType::Pointer(new NodeType(1, 0.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(2, 1.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(3, 0.0, 1.0, 0.0)));
            auto geometry = *GeometryType::Pointer(new Triangle2D3<NodeType>(nodes));
            RunElementSizeCalculatorDerivativesTest(
                geometry, ElementSizeCalculator<2, 3>::MinimumElementSize,
                ElementSizeCalculator<2, 3>::MinimumElementSizeDerivative, 1e-8, 1e-7);
        }

        KRATOS_TEST_CASE_IN_SUITE(Triangle2D3AverageElementSizeDerivatives, KratosCoreFastSuite)
        {
            Geometry<NodeType>::PointsArrayType nodes;
            nodes.push_back(NodeType::Pointer(new NodeType(1, 0.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(2, 1.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(3, 0.0, 1.0, 0.0)));
            auto geometry = *GeometryType::Pointer(new Triangle2D3<NodeType>(nodes));
            RunElementSizeCalculatorDerivativesTest(
                geometry, ElementSizeCalculator<2, 3>::AverageElementSize,
                ElementSizeCalculator<2, 3>::AverageElementSizeDerivative, 1e-8, 1e-7);
        }

        KRATOS_TEST_CASE_IN_SUITE(Quadrilateral2D4MinimumElementSizeDerivatives, KratosCoreFastSuite)
        {
            Geometry<NodeType>::PointsArrayType nodes;
            nodes.push_back(NodeType::Pointer(new NodeType(1, 0.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(2, 1.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(3, 1.0, 1.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(4, 0.0, 1.0, 0.0)));
            auto geometry = *GeometryType::Pointer(new Quadrilateral2D4<NodeType>(nodes));
            RunElementSizeCalculatorDerivativesTest(
                geometry, ElementSizeCalculator<2, 4>::MinimumElementSize,
                ElementSizeCalculator<2, 4>::MinimumElementSizeDerivative, 1e-8, 1e-7);
        }

        KRATOS_TEST_CASE_IN_SUITE(Quadrilateral2D4AverageElementSizeDerivatives, KratosCoreFastSuite)
        {
            Geometry<NodeType>::PointsArrayType nodes;
            nodes.push_back(NodeType::Pointer(new NodeType(1, 0.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(2, 1.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(3, 1.0, 1.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(4, 0.0, 1.0, 0.0)));
            auto geometry = *GeometryType::Pointer(new Quadrilateral2D4<NodeType>(nodes));
            RunElementSizeCalculatorDerivativesTest(
                geometry, ElementSizeCalculator<2, 4>::AverageElementSize,
                ElementSizeCalculator<2, 4>::AverageElementSizeDerivative, 1e-8, 1e-7);
        }

        KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4MinimumElementSizeDerivatives, KratosCoreFastSuite)
        {
            Geometry<NodeType>::PointsArrayType nodes;
            nodes.push_back(NodeType::Pointer(new NodeType(1, 0.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(2, 1.0, 2.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(3, 0.5, 0.5, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(4, 0.2, 0.3, 1.0)));
            auto geometry = *GeometryType::Pointer(new Tetrahedra3D4<NodeType>(nodes));
            RunElementSizeCalculatorDerivativesTest(
                geometry, ElementSizeCalculator<3, 4>::MinimumElementSize,
                ElementSizeCalculator<3, 4>::MinimumElementSizeDerivative, 1e-8, 1e-7);
        }

        KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4AverageElementSizeDerivatives, KratosCoreFastSuite)
        {
            Geometry<NodeType>::PointsArrayType nodes;
            nodes.push_back(NodeType::Pointer(new NodeType(1, 0.0, 0.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(2, 1.0, 2.0, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(3, 0.5, 0.5, 0.0)));
            nodes.push_back(NodeType::Pointer(new NodeType(4, 0.2, 0.3, 1.0)));
            auto geometry = *GeometryType::Pointer(new Tetrahedra3D4<NodeType>(nodes));
            RunElementSizeCalculatorDerivativesTest(
                geometry, ElementSizeCalculator<3, 4>::AverageElementSize,
                ElementSizeCalculator<3, 4>::AverageElementSizeDerivative, 1e-8, 1e-7);
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
            auto geometry = *GeometryType::Pointer(new Prism3D6<NodeType>(nodes));
            RunElementSizeCalculatorDerivativesTest(
                geometry, ElementSizeCalculator<3, 6>::MinimumElementSize,
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
            auto geometry = *GeometryType::Pointer(new Prism3D6<NodeType>(nodes));
            RunElementSizeCalculatorDerivativesTest(
                geometry, ElementSizeCalculator<3, 6>::AverageElementSize,
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
            auto geometry = *GeometryType::Pointer(new Hexahedra3D8<NodeType>(nodes));
            RunElementSizeCalculatorDerivativesTest(
                geometry, ElementSizeCalculator<3, 8>::MinimumElementSize,
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
            auto geometry = *GeometryType::Pointer(new Hexahedra3D8<NodeType>(nodes));
            RunElementSizeCalculatorDerivativesTest(
                geometry, ElementSizeCalculator<3, 8>::AverageElementSize,
                ElementSizeCalculator<3, 8>::AverageElementSizeDerivative, 1e-8, 1e-7);
        }

    } // namespace Testing

} // namespace Kratos
