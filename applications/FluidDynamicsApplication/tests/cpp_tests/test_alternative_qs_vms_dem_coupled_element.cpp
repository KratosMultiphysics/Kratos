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

// Application includes
#include "custom_constitutive/newtonian_2d_law.h"
#include "custom_elements/alternative_qs_vms_dem_coupled.h"
#include "custom_utilities/qsvms_dem_coupled_data.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(AlternativeQSVMSDEMCoupled2D4N, FluidDynamicsApplicationFastSuite)
{
    Model model;
    unsigned int buffer_size = 2;
    unsigned int Dim = 2;
    ModelPart& model_part = model.CreateModelPart("Main",buffer_size);

    // Variables addition
    model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    model_part.AddNodalSolutionStepVariable(PRESSURE);
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    model_part.AddNodalSolutionStepVariable(ACCELERATION);
    model_part.AddNodalSolutionStepVariable(NODAL_AREA);
    model_part.AddNodalSolutionStepVariable(ADVPROJ);
    model_part.AddNodalSolutionStepVariable(DIVPROJ);
    model_part.AddNodalSolutionStepVariable(REACTION);
    model_part.AddNodalSolutionStepVariable(REACTION_WATER_PRESSURE);
    model_part.AddNodalSolutionStepVariable(MASS_SOURCE);
    model_part.AddNodalSolutionStepVariable(FLUID_FRACTION);
    model_part.AddNodalSolutionStepVariable(FLUID_FRACTION_RATE);
    model_part.AddNodalSolutionStepVariable(FLUID_FRACTION_GRADIENT);
    model_part.AddNodalSolutionStepVariable(PERMEABILITY);

    // Process info creation
    double delta_time = 0.1;
    model_part.GetProcessInfo().SetValue(DELTA_TIME, delta_time);

    // Set the element properties
    Properties::Pointer p_properties = model_part.CreateNewProperties(0);
    p_properties->SetValue(DENSITY, 1000.0);
    p_properties->SetValue(DYNAMIC_VISCOSITY, 1.0e-5);
    p_properties->SetValue(DYNAMIC_TAU, 0.0);
    ConstitutiveLaw::Pointer pConsLaw(new Newtonian2DLaw());
    p_properties->SetValue(CONSTITUTIVE_LAW, pConsLaw);

    // Geometry creation
    model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    model_part.CreateNewNode(2, 1.0, 0.1, 0.0);
    model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    model_part.CreateNewNode(4, 1.0, 0.9, 0.0);

    for (ModelPart::NodeIterator it_node=model_part.NodesBegin(); it_node<model_part.NodesEnd(); ++it_node){
        it_node->AddDof(VELOCITY_X,REACTION_X);
        it_node->AddDof(VELOCITY_Y,REACTION_Y);
        it_node->AddDof(VELOCITY_Z,REACTION_Z);
        it_node->AddDof(PRESSURE,REACTION_WATER_PRESSURE);
        double& r_fluid_fraction = it_node->FastGetSolutionStepValue(FLUID_FRACTION);
        r_fluid_fraction = 1.0;
        Matrix& r_permeability = it_node->FastGetSolutionStepValue(PERMEABILITY);
        r_permeability = ZeroMatrix(Dim, Dim);
        for (unsigned int d = 0; d < Dim; ++d){
            r_permeability(d,d) = 0.0;
        }
    }

    std::vector<ModelPart::IndexType> element_nodes {1, 2, 4, 3};
    model_part.CreateNewElement("AlternativeQSVMSDEMCoupled2D4N", 1, element_nodes, p_properties);

    // Loop starts at 1 because you need one less clone than time steps (JC)
    for (unsigned int i = 1; i < buffer_size; i++) {
        model_part.CloneTimeStep(i * delta_time);
    }

    // Define the nodal values
    Matrix reference_velocity(4,2);
    reference_velocity(0,0) = 0.0; reference_velocity(0,1) = 0.1;
    reference_velocity(1,0) = 0.1; reference_velocity(1,1) = 0.2;
    reference_velocity(2,0) = 0.2; reference_velocity(2,1) = 0.3;
    reference_velocity(3,0) = 0.3; reference_velocity(3,1) = 0.4;


    Geometry<Node>& r_geometry = model_part.ElementsBegin()->GetGeometry();


    for(unsigned int i=0; i<4; i++){
        r_geometry[i].FastGetSolutionStepValue(PRESSURE)    = 0.0;
        r_geometry[i].FastGetSolutionStepValue(PRESSURE, 1) = 0.0;
        for(unsigned int k=0; k<2; k++){
            r_geometry[i].FastGetSolutionStepValue(VELOCITY)[k]    = reference_velocity(i,k);
            r_geometry[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.9*reference_velocity(i,k);
            r_geometry[i].FastGetSolutionStepValue(MESH_VELOCITY)[k]    = 0.0;
            r_geometry[i].FastGetSolutionStepValue(MESH_VELOCITY, 1)[k] = 0.0;
        }
    }

    // RHS and LHS
    Vector RHS = ZeroVector(12);
    Matrix LHS = ZeroMatrix(12,12);

    std::vector<double> output = {-0.7903859428,0.1422482204,-0.01715539615,-10.404201,-4.187299126,-0.05064119836,-21.69732287,-21.35183401,-0.07314514572,-14.77475686,-22.26978176,-0.05905825977}; // QSVMSDEMCoupled2D4N

    for (ModelPart::ElementIterator i = model_part.ElementsBegin(); i != model_part.ElementsEnd(); i++) {
        const auto& r_process_info = model_part.GetProcessInfo();
        i->Initialize(r_process_info); // Initialize constitutive law

        const auto& rElem = *i;
        rElem.Check(r_process_info);
        i->CalculateLocalVelocityContribution(LHS, RHS, r_process_info);
        // std::cout << i->Info() << std::setprecision(10) << std::endl;
        // KRATOS_WATCH(RHS);

        for (unsigned int j = 0; j < output.size(); j++) {
            KRATOS_EXPECT_NEAR(RHS[j], output[j], 1e-4);
        }
    }
    double porosity = 0.5;
    for (ModelPart::NodeIterator it_node=model_part.NodesBegin(); it_node<model_part.NodesEnd(); ++it_node){
        double& r_fluid_fraction = it_node->FastGetSolutionStepValue(FLUID_FRACTION);
        r_fluid_fraction = porosity;
        Matrix& r_permeability = it_node->FastGetSolutionStepValue(PERMEABILITY);
        r_permeability = ZeroMatrix(Dim, Dim);
    }

    for (ModelPart::ElementIterator i = model_part.ElementsBegin(); i != model_part.ElementsEnd(); i++) {
        const auto& r_process_info = model_part.GetProcessInfo();
        i->Initialize(r_process_info); // Initialize constitutive law
        const auto& rElem = *i;
        rElem.Check(r_process_info);
        i->CalculateLocalVelocityContribution(LHS, RHS, r_process_info);

        for (unsigned int j = 0; j < output.size(); j++) {
            KRATOS_CHECK_NEAR(RHS[j], porosity*output[j], 1e-5);
        }
    }
}

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

        AlternativeQSVMSDEMCoupled<QSVMSDEMCoupledData<3, 27>>* p_element = dynamic_cast<AlternativeQSVMSDEMCoupled<QSVMSDEMCoupledData<3, 27>>*>(&rElement);
        const Geometry<Node>::IntegrationMethod integration_method = p_element->GetIntegrationMethod();


        GeometryUtils::ShapeFunctionsSecondDerivativesTransformOnAllIntegrationPoints(
            DDN_DDX,rElement.GetGeometry(),rElement.GetIntegrationMethod());

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
                    KRATOS_CHECK_NEAR(hess(d,e), exact_hess(d,e),1e-10);
                }
            }
        }

    }
}

KRATOS_TEST_CASE_IN_SUITE(SecondShapeDerivativesInterpolation, FluidDynamicsApplicationFastSuite)
{
    Model model;
    const double buffer_size = 2;
    const unsigned int Dim = 2;
    ModelPart& model_part = model.CreateModelPart("Main",buffer_size);

    // Set the element properties
    Properties::Pointer p_properties = model_part.CreateNewProperties(0);

    // Square creation
    model_part.CreateNewNode(1, 0.0, 0.0125, 0.0);
    model_part.CreateNewNode(2, 0.0, 0.00625, 0.0);
    model_part.CreateNewNode(3, 0.00625,0.0125,0.0);
    model_part.CreateNewNode(4, 0.00625, 0.00625, 0.0);
    model_part.CreateNewNode(5, 0.0, 0.0, 0.0);
    model_part.CreateNewNode(6, 0.00625, 0.0, 0.0);
    model_part.CreateNewNode(7, 0.0125, 0.0125, 0.0);
    model_part.CreateNewNode(8, 0.0125, 0.00625, 0.0);
    model_part.CreateNewNode(9, 0.0125, 0.0, 0.0);


    std::vector<ModelPart::IndexType> element_nodes {9, 7, 1, 5, 8, 3, 2, 6, 4};
    model_part.CreateNewElement("AlternativeQSVMSDEMCoupled2D9N", 1, element_nodes, p_properties);


    for (ModelPart::ElementIterator i_elem = model_part.ElementsBegin(); i_elem != model_part.ElementsEnd(); i_elem++) {
        auto& rElement = *i_elem;

        Matrix exact_hess(Dim,Dim);

        DenseVector<DenseVector<Matrix>> DDN_DDX;

        const unsigned int NumNodes = rElement.GetGeometry().PointsNumber();

        AlternativeQSVMSDEMCoupled<QSVMSDEMCoupledData<Dim, 9>>* p_element = dynamic_cast<AlternativeQSVMSDEMCoupled<QSVMSDEMCoupledData<Dim, 9>>*>(&rElement);
        const Geometry<Node>::IntegrationMethod integration_method = p_element->GetIntegrationMethod();

        GeometryUtils::ShapeFunctionsSecondDerivativesTransformOnAllIntegrationPoints(
            DDN_DDX, p_element->GetGeometry(), p_element->GetIntegrationMethod());

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

                const auto r_pressure = std::pow(x1,2)*std::pow(x2,2) + std::pow(x1,2)*x2 + std::pow(x2,2)*x1 + std::pow(x1,2) + std::pow(x2,2) + x1*x2 + x1 + x2 + 1.0;
                for (unsigned int d = 0; d < Dim; ++d){
                    gauss_point_coordinates(g,d) += NContainer(g,i) * r_coordinates[d];
                    for (unsigned int e = 0; e < Dim; ++e){
                        hess(d,e) += DDN_DDX[g][i](d,e) * r_pressure;
                        }
                    }
            }

            const double x1 = gauss_point_coordinates(g,0);
            const double x2 = gauss_point_coordinates(g,1);

            exact_hess(0,0) = 2.0 * std::pow(x2,2) + 2.0 * x2 + 2.0;
            exact_hess(0,1) = 4.0 * x1 * x2 + 2.0 * x1 + 2.0 * x2 + 1.0;
            exact_hess(1,0) = 4.0 * x1 * x2 + 2.0 * x1 + 2.0 * x2 + 1.0;
            exact_hess(1,1) = 2.0 * std::pow(x1,2) + 2.0 * x1 + 2.0;

            for (unsigned int d = 0; d < Dim; ++d){
                for (unsigned int e = 0; e < Dim; ++e){
                    KRATOS_CHECK_NEAR(hess(d,e), exact_hess(d,e),1e-10);
                }
            }
        }

    }
}

}  // namespace Testing
}  // namespace Kratos