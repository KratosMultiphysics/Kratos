//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "includes/model_part.h"

// Application includes
#include "custom_utilities/fluid_calculation_utilities.h"
#include "tests/cpp_tests/fluid_dynamics_fast_suite.h"

namespace Kratos {
namespace Testing {

    Element& GenerateTestModelPart(ModelPart& rModelPart, const unsigned int Dim) {

        // Set buffer size
        rModelPart.SetBufferSize(2);

        // Variables addition
        rModelPart.AddNodalSolutionStepVariable(DENSITY);
        rModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);

        // Process info creation
        rModelPart.GetProcessInfo()[DOMAIN_SIZE] =  Dim;
        auto p_elem_prop = rModelPart.CreateNewProperties(1);

        // Element creation
        auto p_node_1 = rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
        auto p_node_2 = rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
        auto p_node_3 = rModelPart.CreateNewNode(3, 0.0, 1.0, 0.0);
        auto p_node_4 = rModelPart.CreateNewNode(4, 0.0, 1.0, 1.0);

        // set variable values
        for (auto& r_node : rModelPart.Nodes()) {
            const double node_id = r_node.Id();

            r_node.FastGetSolutionStepValue(DENSITY) = node_id * 2.1;
            r_node.FastGetSolutionStepValue(DISPLACEMENT) = array_1d<double, 3>({node_id * 3.4, node_id * node_id + 3, node_id * 2 + 4 });

            r_node.FastGetSolutionStepValue(DENSITY, 1) = node_id * 4.9;
            r_node.FastGetSolutionStepValue(DISPLACEMENT, 1) = array_1d<double, 3>({node_id * 6.1, node_id * node_id * 5 + 3, node_id * node_id * 2 + 4 });

            r_node.SetValue(DENSITY, node_id * 5);
            r_node.SetValue(DISPLACEMENT, array_1d<double, 3>({node_id, node_id * 6.7, node_id * 4 + 6}));
        }

        if (Dim == 2) {
            return *rModelPart.CreateNewElement("Element2D3N", 1, {1, 2, 3}, p_elem_prop);
        } else {
            return *rModelPart.CreateNewElement("Element3D4N", 1, {1, 2, 3, 4}, p_elem_prop);
        }
    }

    void CalculateGeometryData(
        const Element::GeometryType& rGeometry,
        Matrix& rNContainer,
        Element::GeometryType::ShapeFunctionsGradientsType& rDN_DX)
    {
        const auto number_of_nodes = rGeometry.PointsNumber();
        const auto integration_method = GeometryData::IntegrationMethod::GI_GAUSS_2;
        const IndexType number_of_gauss_points = rGeometry.IntegrationPointsNumber(integration_method);

        Vector DetJ;
        rGeometry.ShapeFunctionsIntegrationPointsGradients(rDN_DX, DetJ, integration_method);

        if (rNContainer.size1() != number_of_gauss_points || rNContainer.size2() != number_of_nodes) {
            rNContainer.resize(number_of_gauss_points, number_of_nodes, false);
        }

        rNContainer = rGeometry.ShapeFunctionsValues(integration_method);
    }

    template<unsigned int TDim, class EvaluationMethodType, class NodalDataRetrievalType>
    void TestEvaluateInPoint(EvaluationMethodType&& rEvaluationMethod, NodalDataRetrievalType&& rNodalDataRetrievalType)
    {
        // Create a test element inside a modelpart
        Model model;
        ModelPart& model_part = model.CreateModelPart("Main");
        const auto& r_geometry = GenerateTestModelPart(model_part, TDim).GetGeometry();

        // calculate gauss point data
        Matrix Ns;
        Element::GeometryType::ShapeFunctionsGradientsType dNdXs;
        CalculateGeometryData(r_geometry, Ns, dNdXs);

        const Vector& N = row(Ns, 0);

        double density;
        array_1d<double, TDim> displacement;

        rEvaluationMethod(r_geometry, N, density, displacement);

        double check_density{0.0};
        array_1d<double, TDim> check_displacement = ZeroVector(TDim);

        for (std::size_t i = 0; i < r_geometry.PointsNumber(); ++i) {
            const auto& r_node = r_geometry[i];

            check_density += rNodalDataRetrievalType(r_node, DENSITY) * N[i];

            array_1d<double, 3> nodal_displacement;
            nodal_displacement[0] = rNodalDataRetrievalType(r_node, DISPLACEMENT_X);
            nodal_displacement[1] = rNodalDataRetrievalType(r_node, DISPLACEMENT_Y);
            nodal_displacement[2] = rNodalDataRetrievalType(r_node, DISPLACEMENT_Z);

            for (std::size_t dim = 0; dim < TDim; ++dim) {
                check_displacement[dim] += nodal_displacement[dim] * N[i];
            }
        }

        KRATOS_EXPECT_NEAR(density, check_density, 1e-16);
        KRATOS_EXPECT_VECTOR_NEAR(displacement, check_displacement, 1e-16);
    }

    template<unsigned int TDim, class EvaluationMethodType, class NodalDataRetrievalType>
    void TestEvaluateGradientInPoint(EvaluationMethodType&& rEvaluationMethod, NodalDataRetrievalType&& rNodalDataRetrievalType)
    {
        // Create a test element inside a modelpart
        Model model;
        ModelPart& model_part = model.CreateModelPart("Main");
        const auto& r_geometry = GenerateTestModelPart(model_part, TDim).GetGeometry();

        // calculate gauss point data
        Matrix Ns;
        Element::GeometryType::ShapeFunctionsGradientsType dNdXs;
        CalculateGeometryData(r_geometry, Ns, dNdXs);

        const Matrix& dNdX = dNdXs[0];
        array_1d<double, TDim> density_gradient;
        BoundedMatrix<double, TDim, TDim> displacement_gradient;

        rEvaluationMethod(r_geometry, dNdX, density_gradient, displacement_gradient);


        array_1d<double, TDim> check_density_gradient = ZeroVector(TDim);
        BoundedMatrix<double, TDim, TDim> check_displacement_gradient = ZeroMatrix(TDim, TDim);

        for (std::size_t i = 0; i < r_geometry.PointsNumber(); ++i) {
            const auto& r_node = r_geometry[i];

            const auto density = rNodalDataRetrievalType(r_node, DENSITY);

            array_1d<double, 3> nodal_displacement;
            nodal_displacement[0] = rNodalDataRetrievalType(r_node, DISPLACEMENT_X);
            nodal_displacement[1] = rNodalDataRetrievalType(r_node, DISPLACEMENT_Y);
            nodal_displacement[2] = rNodalDataRetrievalType(r_node, DISPLACEMENT_Z);

            const Vector& r_shape_derivatives = row(dNdX, i);
            for (std::size_t j = 0; j < TDim; ++j) {
                check_density_gradient[j] += density * r_shape_derivatives[j];


                for (std::size_t k = 0; k < TDim; ++k) {
                    check_displacement_gradient(k, j) += nodal_displacement[k] * r_shape_derivatives[j];
                }
            }
        }

        KRATOS_EXPECT_VECTOR_NEAR(density_gradient, check_density_gradient, 1e-16);
        KRATOS_EXPECT_MATRIX_NEAR(displacement_gradient, check_displacement_gradient, 1e-16);
    }

    KRATOS_TEST_CASE_IN_SUITE(FluidCalculationUtilitiesEvaluateInPoint2D, FluidDynamicsApplicationFastSuite)
    {
        TestEvaluateInPoint<2>([](const Element::GeometryType& rGeometry, const Vector& rN, double& rDensity, array_1d<double, 2>& rDisplacement) {
            FluidCalculationUtilities::EvaluateInPoint(rGeometry, rN,
                                                        std::tie(rDensity, DENSITY),
                                                        std::tie(rDisplacement, DISPLACEMENT));
        }, [](const ModelPart::NodeType& rNode, const Variable<double>& rVariable) {
            return rNode.FastGetSolutionStepValue(rVariable);
        });

        TestEvaluateInPoint<2>([](const Element::GeometryType& rGeometry, const Vector& rN, double& rDensity, array_1d<double, 2>& rDisplacement) {
            FluidCalculationUtilities::EvaluateInPoint(rGeometry, rN, 1,
                                                        std::tie(rDensity, DENSITY),
                                                        std::tie(rDisplacement, DISPLACEMENT));
        }, [](const ModelPart::NodeType& rNode, const Variable<double>& rVariable) {
            return rNode.FastGetSolutionStepValue(rVariable, 1);
        });
    }

    KRATOS_TEST_CASE_IN_SUITE(FluidCalculationUtilitiesEvaluateInPoint3D, FluidDynamicsApplicationFastSuite)
    {
        TestEvaluateInPoint<3>([](const Element::GeometryType& rGeometry, const Vector& rN, double& rDensity, array_1d<double, 3>& rDisplacement) {
            FluidCalculationUtilities::EvaluateInPoint(rGeometry, rN,
                                                        std::tie(rDensity, DENSITY),
                                                        std::tie(rDisplacement, DISPLACEMENT));
        }, [](const ModelPart::NodeType& rNode, const Variable<double>& rVariable) {
            return rNode.FastGetSolutionStepValue(rVariable);
        });

        TestEvaluateInPoint<3>([](const Element::GeometryType& rGeometry, const Vector& rN, double& rDensity, array_1d<double, 3>& rDisplacement) {
            FluidCalculationUtilities::EvaluateInPoint(rGeometry, rN, 1,
                                                        std::tie(rDensity, DENSITY),
                                                        std::tie(rDisplacement, DISPLACEMENT));
        }, [](const ModelPart::NodeType& rNode, const Variable<double>& rVariable) {
            return rNode.FastGetSolutionStepValue(rVariable, 1);
        });
    }

    KRATOS_TEST_CASE_IN_SUITE(FluidCalculationUtilitiesEvaluateNonHistoricalInPoint2D, FluidDynamicsApplicationFastSuite)
    {
        TestEvaluateInPoint<2>([](const Element::GeometryType& rGeometry, const Vector& rN, double& rDensity, array_1d<double, 2>& rDisplacement) {
            FluidCalculationUtilities::EvaluateNonHistoricalInPoint(rGeometry, rN,
                                                        std::tie(rDensity, DENSITY),
                                                        std::tie(rDisplacement, DISPLACEMENT));
        }, [](const ModelPart::NodeType& rNode, const Variable<double>& rVariable) {
            return rNode.GetValue(rVariable);
        });

    }

    KRATOS_TEST_CASE_IN_SUITE(FluidCalculationUtilitiesEvaluateNonHistoricalInPoint3D, FluidDynamicsApplicationFastSuite)
    {
        TestEvaluateInPoint<3>([](const Element::GeometryType& rGeometry, const Vector& rN, double& rDensity, array_1d<double, 3>& rDisplacement) {
            FluidCalculationUtilities::EvaluateNonHistoricalInPoint(rGeometry, rN,
                                                        std::tie(rDensity, DENSITY),
                                                        std::tie(rDisplacement, DISPLACEMENT));
        }, [](const ModelPart::NodeType& rNode, const Variable<double>& rVariable) {
            return rNode.GetValue(rVariable);
        });
    }


    KRATOS_TEST_CASE_IN_SUITE(FluidCalculationUtilitiesEvaluateGradientInPoint2D, FluidDynamicsApplicationFastSuite)
    {
        TestEvaluateGradientInPoint<2>([](const Element::GeometryType& rGeometry, const Matrix& rdNdX, array_1d<double, 2>& rDensityGradient, BoundedMatrix<double, 2, 2>& rDisplacementGradient) {
            FluidCalculationUtilities::EvaluateGradientInPoint(rGeometry, rdNdX,
                                                                std::tie(rDensityGradient, DENSITY),
                                                                std::tie(rDisplacementGradient, DISPLACEMENT));
        }, [](const ModelPart::NodeType& rNode, const Variable<double>& rVariable) {
            return rNode.FastGetSolutionStepValue(rVariable);
        });

        TestEvaluateGradientInPoint<2>([](const Element::GeometryType& rGeometry, const Matrix& rdNdX, array_1d<double, 2>& rDensityGradient, BoundedMatrix<double, 2, 2>& rDisplacementGradient) {
            FluidCalculationUtilities::EvaluateGradientInPoint(rGeometry, rdNdX, 1,
                                                                std::tie(rDensityGradient, DENSITY),
                                                                std::tie(rDisplacementGradient, DISPLACEMENT));
        }, [](const ModelPart::NodeType& rNode, const Variable<double>& rVariable) {
            return rNode.FastGetSolutionStepValue(rVariable, 1);
        });
    }

    KRATOS_TEST_CASE_IN_SUITE(FluidCalculationUtilitiesEvaluateGradientInPoint3D, FluidDynamicsApplicationFastSuite)
    {
        TestEvaluateGradientInPoint<3>([](const Element::GeometryType& rGeometry, const Matrix& rdNdX, array_1d<double, 3>& rDensityGradient, BoundedMatrix<double, 3, 3>& rDisplacementGradient) {
            FluidCalculationUtilities::EvaluateGradientInPoint(rGeometry, rdNdX,
                                                                std::tie(rDensityGradient, DENSITY),
                                                                std::tie(rDisplacementGradient, DISPLACEMENT));
        }, [](const ModelPart::NodeType& rNode, const Variable<double>& rVariable) {
            return rNode.FastGetSolutionStepValue(rVariable);
        });

        TestEvaluateGradientInPoint<3>([](const Element::GeometryType& rGeometry, const Matrix& rdNdX, array_1d<double, 3>& rDensityGradient, BoundedMatrix<double, 3, 3>& rDisplacementGradient) {
            FluidCalculationUtilities::EvaluateGradientInPoint(rGeometry, rdNdX, 1,
                                                                std::tie(rDensityGradient, DENSITY),
                                                                std::tie(rDisplacementGradient, DISPLACEMENT));
        }, [](const ModelPart::NodeType& rNode, const Variable<double>& rVariable) {
            return rNode.FastGetSolutionStepValue(rVariable, 1);
        });
    }

    KRATOS_TEST_CASE_IN_SUITE(FluidCalculationUtilitiesEvaluateNonHistoricalGradientInPoint2D, FluidDynamicsApplicationFastSuite)
    {
        TestEvaluateGradientInPoint<2>([](const Element::GeometryType& rGeometry, const Matrix& rdNdX, array_1d<double, 2>& rDensityGradient, BoundedMatrix<double, 2, 2>& rDisplacementGradient) {
            FluidCalculationUtilities::EvaluateNonHistoricalGradientInPoint(rGeometry, rdNdX,
                                                                std::tie(rDensityGradient, DENSITY),
                                                                std::tie(rDisplacementGradient, DISPLACEMENT));
        }, [](const ModelPart::NodeType& rNode, const Variable<double>& rVariable) {
            return rNode.GetValue(rVariable);
        });
    }

    KRATOS_TEST_CASE_IN_SUITE(FluidCalculationUtilitiesEvaluateNonHistoricalGradientInPoint3D, FluidDynamicsApplicationFastSuite)
    {
        TestEvaluateGradientInPoint<3>([](const Element::GeometryType& rGeometry, const Matrix& rdNdX, array_1d<double, 3>& rDensityGradient, BoundedMatrix<double, 3, 3>& rDisplacementGradient) {
            FluidCalculationUtilities::EvaluateNonHistoricalGradientInPoint(rGeometry, rdNdX,
                                                                std::tie(rDensityGradient, DENSITY),
                                                                std::tie(rDisplacementGradient, DISPLACEMENT));
        }, [](const ModelPart::NodeType& rNode, const Variable<double>& rVariable) {
            return rNode.GetValue(rVariable);
        });
    }

} // namespace Testing
}  // namespace Kratos.
