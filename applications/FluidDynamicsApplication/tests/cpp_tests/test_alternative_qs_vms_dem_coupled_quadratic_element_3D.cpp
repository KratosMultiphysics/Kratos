//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Joaquin Gonzalez-Usua
//

#include <iomanip> // for std::setprecision

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "includes/cfd_variables.h"
#include "fluid_dynamics_application_variables.h"
#include "custom_elements/alternative_qs_vms_dem_coupled.h"
#include "custom_utilities/qsvms_dem_coupled_data.h"
// Application includes
#include "custom_constitutive/newtonian_3d_law.h"

namespace Kratos {

namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(SecondShapeDerivativesInterpolation3D, FluidDynamicsApplicationFastSuite)
{
    Model model;
    const double buffer_size = 2;
    const unsigned int Dim = 3;
    ModelPart& model_part = model.CreateModelPart("Main",buffer_size);

    // Set the element properties
    Properties::Pointer p_properties = model_part.CreateNewProperties(0);

    // Square creation
    model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    model_part.CreateNewNode(2, 0.0, 0.0, 0.05);
    model_part.CreateNewNode(3, 0.05,0.0,0.0);
    model_part.CreateNewNode(4, 0.0, 0.05, 0.0);
    model_part.CreateNewNode(5, 0.05, 0.0, 0.05);
    model_part.CreateNewNode(6, 0.05, 0.05, 0.0);
    model_part.CreateNewNode(7, 0.0, 0.05, 0.05);
    model_part.CreateNewNode(8, 0.05, 0.05, 0.05);
    model_part.CreateNewNode(9, 0.0, 0.1, 0.0);
    model_part.CreateNewNode(10, 0.0, 0.0, 0.1);
    model_part.CreateNewNode(11, 0.1, 0.0, 0.0);
    model_part.CreateNewNode(12, 0.05, 0.1, 0.0);
    model_part.CreateNewNode(13, 0.0, 0.1, 0.05);
    model_part.CreateNewNode(14, 0.0, 0.05, 0.1);
    model_part.CreateNewNode(15, 0.05, 0.0, 0.1);
    model_part.CreateNewNode(16, 0.1, 0.0, 0.05);
    model_part.CreateNewNode(17, 0.1, 0.05, 0.0);
    model_part.CreateNewNode(18, 0.05, 0.1, 0.05);
    model_part.CreateNewNode(19, 0.05, 0.05, 0.1);
    model_part.CreateNewNode(20, 0.1, 0.05, 0.05);
    model_part.CreateNewNode(21, 0.1, 0.1, 0.0);
    model_part.CreateNewNode(22, 0.1, 0.0, 0.1);
    model_part.CreateNewNode(23, 0.0, 0.1, 0.1);
    model_part.CreateNewNode(25, 0.1, 0.1, 0.05);
    model_part.CreateNewNode(26, 0.05, 0.1, 0.1);
    model_part.CreateNewNode(29, 0.1, 0.05, 0.1);
    model_part.CreateNewNode(39, 0.1, 0.1, 0.1);


    std::vector<ModelPart::IndexType> element_nodes {39, 23, 9, 21, 22, 10, 1, 11, 26, 13, 12, 25, 29, 14, 4, 17, 15, 2, 3, 16, 18, 19, 7, 6, 20, 5, 8};
    model_part.CreateNewElement("AlternativeQSVMSDEMCoupled3D27N", 1, element_nodes, p_properties);


    for (ModelPart::ElementIterator i_elem = model_part.ElementsBegin(); i_elem != model_part.ElementsEnd(); i_elem++) {
        auto& rElement = *i_elem;

        Matrix exact_hess(Dim,Dim);

        DenseVector<DenseVector<Matrix>> DDN_DDX;

        const unsigned int NumNodes = rElement.GetGeometry().PointsNumber();

        AlternativeQSVMSDEMCoupled<QSVMSDEMCoupledData<Dim, 27>>* p_element = dynamic_cast<AlternativeQSVMSDEMCoupled<QSVMSDEMCoupledData<Dim, 27>>*>(&rElement);
        const Geometry<Node>::IntegrationMethod integration_method = p_element->GetIntegrationMethod();

        p_element->GetShapeSecondDerivatives(DDN_DDX);

        Matrix NContainer = p_element->GetGeometry().ShapeFunctionsValues(integration_method);

        const Geometry<Node>::IntegrationPointsArrayType integration_points = p_element->GetGeometry().IntegrationPoints(integration_method);
        const unsigned int number_of_integration_points = integration_points.size();

        Matrix gauss_point_coordinates = ZeroMatrix(number_of_integration_points,Dim);

        for(unsigned int g = 0; g < number_of_integration_points; g++){
            Matrix hess = ZeroMatrix(Dim, Dim);
            Matrix Hessian(Dim, Dim);
            for (unsigned int i = 0; i < NumNodes; ++i){
                const array_1d<double, 3>& r_coordinates = p_element->GetGeometry()[i].Coordinates();
                const double x1 = r_coordinates[0];
                const double x2 = r_coordinates[1];
                const double x3 = r_coordinates[2];

                const auto r_pressure = std::pow(x1,2)*std::pow(x2,2) + std::pow(x1,2)*std::pow(x3,2) + std::pow(x2,2)*std::pow(x3,2) + std::pow(x1,2)*x2 + std::pow(x1,2)*x3 + std::pow(x2,2)*x1 + std::pow(x2,2)*x3 + std::pow(x3,2)*x1 + std::pow(x3,2)*x2 + std::pow(x1,2) + std::pow(x2,2) + std::pow(x3,2) + x1*x2*x3 + x1*x2 + x2*x3 + x1*x3 + x1 + x2 + x3 + 1.0;
                for (unsigned int d = 0; d < Dim; ++d){
                    gauss_point_coordinates(g,d) += NContainer(g,i) * r_coordinates[d];
                    for (unsigned int e = 0; e < Dim; ++e){
                        hess(d,e) += DDN_DDX[g][i](d,e) * r_pressure;
                        }
                    }
            }

            const double x1 = gauss_point_coordinates(g,0);
            const double x2 = gauss_point_coordinates(g,1);
            const double x3 = gauss_point_coordinates(g,2);

            exact_hess(0,0) = 2.0 * std::pow(x2,2) + 2.0 * std::pow(x3,2) + 2.0 * x2 + 2.0 * x3 + 2.0;
            exact_hess(0,1) = 4.0 * x1 * x2 + 2.0 * x1 + 2.0 * x2 + x3 + 1.0;
            exact_hess(0,2) = 4.0 * x1 * x3 + 2.0 * x1 + 2.0 * x3 + x2 + 1.0;
            exact_hess(1,0) = 4.0 * x1 * x2 + 2.0 * x1 + 2.0 * x2 + x3 + 1.0;
            exact_hess(1,1) = 2.0 * std::pow(x1,2) + 2.0 * std::pow(x3,2) + 2.0 * x1 + 2.0 * x3 + 2.0;
            exact_hess(1,2) = 4.0 * x2 * x3 + 2.0 * x2 + 2.0 * x3 + x1 + 1.0;
            exact_hess(2,0) = 4.0 * x1 * x3 + 2.0 * x1 + 2.0 * x3 + x2 + 1.0;
            exact_hess(2,1) = 4.0 * x3 * x2 + 2.0 * x3 + 2.0 * x2 + x1 + 1.0;
            exact_hess(2,2) = 2.0 * std::pow(x1,2) + 2.0 * std::pow(x2,2) + 2.0 * x1 + 2.0 * x2 + 2.0;

            for (unsigned int d = 0; d < Dim; ++d){
                for (unsigned int e = 0; e < Dim; ++e){
                    KRATOS_EXPECT_NEAR(hess(d,e), exact_hess(d,e),1e-10);
                }
            }
        }

    }
}

}  // namespace Testing
}  // namespace Kratos