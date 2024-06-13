// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Klaus B. Sautter
//

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "structural_mechanics_application_variables.h"

#include "custom_elements/truss_element_3D2N.hpp"
#include "custom_elements/truss_element_linear_3D2N.hpp"

namespace
{

using namespace Kratos;

class StubBilinearLaw : public ConstitutiveLaw
{
public:
    // Only implement the interface that is needed by the test
  StubBilinearLaw(double Strain, double TangentModulus1, double TangentModulus2) :
 mStrain{Strain}, mTangentModuli{TangentModulus1, TangentModulus2}
  {}

    ConstitutiveLaw::Pointer Clone() const override
    {
        return std::make_shared<StubBilinearLaw>(*this);
    }

    SizeType GetStrainSize() const override
    {
        return 3;
    }

    double& CalculateValue(Parameters& rParameterValues, const Variable<double>& rThisVariable, double& rValue) override
    {
        KRATOS_ERROR_IF_NOT(rThisVariable == TANGENT_MODULUS);

        rValue = rParameterValues.GetStrainVector()[0] < mStrain ? mTangentModuli[0] : mTangentModuli[1];
        return rValue;
    }

private:
    double mStrain = 0.0;
    array_1d<double, 2> mTangentModuli{2.0, 1.0};
};

ModelPart& CreateTestModelPart(Model& rModel)
{
  auto&r_result = rModel.CreateModelPart("ModelPart", 1);
  r_result.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);
  r_result.AddNodalSolutionStepVariable(DISPLACEMENT);
  return r_result;
}

}


namespace Kratos
{
namespace Testing
{

    void AddDisplacementDofsElement(ModelPart& rModelPart){
        for (auto& r_node : rModelPart.Nodes()){
            r_node.AddDof(DISPLACEMENT_X);
            r_node.AddDof(DISPLACEMENT_Y);
            r_node.AddDof(DISPLACEMENT_Z);
        }
    }

    // Tests the mass matrix of the TrussElement3D2N
    KRATOS_TEST_CASE_IN_SUITE(TrussElement3D2NMassMatrix, KratosStructuralMechanicsFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart",1);
        r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);

        const double density = 7850.0;
        const double length  = 2.0;
        const double area = 0.01;

        // Set the element properties
        auto p_elem_prop = r_model_part.CreateNewProperties(0);
        p_elem_prop->SetValue(YOUNG_MODULUS, 2.0e+06);
        p_elem_prop->SetValue(DENSITY, density);
        p_elem_prop->SetValue(CROSS_AREA, area);
        const auto &r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("TrussConstitutiveLaw");
        p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

        // Create the test element
        auto p_node_1 = r_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
        auto p_node_2 = r_model_part.CreateNewNode(2, length , 0.0 , 0.0);

        AddDisplacementDofsElement(r_model_part);

        std::vector<ModelPart::IndexType> element_nodes {1,2};
        auto p_element = r_model_part.CreateNewElement("TrussElement3D2N", 1, element_nodes, p_elem_prop);

        const auto& r_process_info = r_model_part.GetProcessInfo();

        p_element->Initialize(r_process_info); // Initialize the element to initialize the constitutive law
        const auto& r_const_elem_ref = *p_element;
        r_const_elem_ref.Check(r_process_info);

        const unsigned int number_of_nodes = p_element->GetGeometry().size();
        const unsigned int dimension = p_element->GetGeometry().WorkingSpaceDimension();
        const unsigned int number_of_dofs = number_of_nodes * dimension;



        p_elem_prop->SetValue(COMPUTE_LUMPED_MASS_MATRIX,true);
        Matrix mm_lumped = ZeroMatrix(number_of_dofs,number_of_dofs);
        p_element->CalculateMassMatrix(mm_lumped,r_process_info);

        p_elem_prop->SetValue(COMPUTE_LUMPED_MASS_MATRIX,false);
        Matrix mm_consistent = ZeroMatrix(number_of_dofs,number_of_dofs);
        p_element->CalculateMassMatrix(mm_consistent,r_process_info);

        for (unsigned int i=0; i<number_of_dofs;++i){
            double diagonal_entry = 0.0;
            for (unsigned int j=0; j<number_of_dofs;++j) diagonal_entry += mm_consistent(i,j);
            KRATOS_EXPECT_NEAR(mm_lumped(i,i),diagonal_entry,1.0e-10);
        }


        Matrix mm_consistent_analytical = ZeroMatrix(number_of_dofs);
        mm_consistent_analytical(0, 0) = 2.0;
        mm_consistent_analytical(0, 3) = 1.0;
        mm_consistent_analytical(1, 1) = 2.0;
        mm_consistent_analytical(1, 4) = 1.0;
        mm_consistent_analytical(2, 2) = 2.0;
        mm_consistent_analytical(2, 5) = 1.0;

        mm_consistent_analytical(3, 0) = 1.0;
        mm_consistent_analytical(3, 3) = 2.0;
        mm_consistent_analytical(4, 1) = 1.0;
        mm_consistent_analytical(4, 4) = 2.0;
        mm_consistent_analytical(5, 2) = 1.0;
        mm_consistent_analytical(5, 5) = 2.0;

        mm_consistent_analytical *= density*area*length / 6.0;

        KRATOS_EXPECT_MATRIX_NEAR(mm_consistent_analytical, mm_consistent, 1e-10);


        Vector lumped_mass_vector = ZeroVector(number_of_dofs);
        p_element->CalculateLumpedMassVector(lumped_mass_vector,r_process_info);

        const double lumped_mass = area*length*density*0.5;
        for (unsigned int i=0;i<number_of_dofs;++i)
        {
            KRATOS_EXPECT_NEAR(mm_lumped(i,i),lumped_mass_vector[i],1.0e-10);
            KRATOS_EXPECT_NEAR(lumped_mass,lumped_mass_vector[i],1.0e-10);
        }

    }

    // Tests the dead load of the TrussElement3D2N
    KRATOS_TEST_CASE_IN_SUITE(TrussElement3D2NDeadLoad, KratosStructuralMechanicsFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart",1);
        r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);

        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);

        const double density = 7850.0;
        const double length  = 2.0;
        const double area = 0.01;

        // Set the element properties
        auto p_elem_prop = r_model_part.CreateNewProperties(0);
        p_elem_prop->SetValue(YOUNG_MODULUS, 2.0e+06);
        p_elem_prop->SetValue(DENSITY, density);
        p_elem_prop->SetValue(CROSS_AREA, area);
        array_1d<double, 3> gravity = ZeroVector(3);
        gravity[0] = 1.0;
        gravity[1] = 2.0;
        gravity[2] = 3.0;
        p_elem_prop->SetValue(VOLUME_ACCELERATION,gravity);
        const auto &r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("TrussConstitutiveLaw");
        p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

        // Create the test element
        auto p_node_1 = r_model_part.CreateNewNode(1, 0.0 , 0.0 , 0.0);
        auto p_node_2 = r_model_part.CreateNewNode(2, length , length , length);

        AddDisplacementDofsElement(r_model_part);

        std::vector<ModelPart::IndexType> element_nodes {1,2};
        auto p_element = r_model_part.CreateNewElement("TrussElement3D2N", 1, element_nodes, p_elem_prop);

        const auto& r_process_info = r_model_part.GetProcessInfo();

        p_element->Initialize(r_process_info); // Initialize the element to initialize the constitutive law
        const auto& r_const_elem_ref = *p_element;
        r_const_elem_ref.Check(r_process_info);

        const unsigned int number_of_nodes = p_element->GetGeometry().size();
        const unsigned int dimension = p_element->GetGeometry().WorkingSpaceDimension();
        const unsigned int number_of_dofs = number_of_nodes * dimension;

        for (unsigned int i=0;i<number_of_nodes;++i){
            array_1d<double, 3>& r_current_acceleration = p_element->GetGeometry()[i].FastGetSolutionStepValue(VOLUME_ACCELERATION);
            r_current_acceleration = gravity;
        }

        Vector rhs = ZeroVector(number_of_dofs);
        p_element->CalculateRightHandSide(rhs,r_process_info);


        const double m1 = 0.50 * length * area * density * sqrt(3.0);

        for (unsigned int i=0;i<number_of_nodes;++i){
            for (unsigned int j=0;j<dimension;++j){
                const unsigned int index = (i*dimension)+j;
                KRATOS_EXPECT_NEAR(
                    rhs[index],
                    m1*gravity[j],
                    1.0e-10);
            }
        }


    }

    KRATOS_TEST_CASE_IN_SUITE(TangentModulusOfTrussElement3D2NUsesGreenLagrangeStrain, KratosStructuralMechanicsFastSuite)
    {
        Model current_model;
        auto& r_model_part = CreateTestModelPart(current_model);

        // Create two nodes and a truss element
        constexpr auto length = 2.0;
        auto p_bottom_node = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        auto p_top_node = r_model_part.CreateNewNode(2, 0.0, 0.0, length);
        AddDisplacementDofsElement(r_model_part);

        auto p_elem_prop = r_model_part.CreateNewProperties(0);
        constexpr auto elongation = 0.01;
        constexpr auto linear_strain = elongation / length;
        constexpr auto tangent_modulus_1 = 2.0e+03;
        constexpr auto tangent_modulus_2 = 1.0e+03;
        auto p_test_law = std::make_shared<StubBilinearLaw>(1.0001 * linear_strain, tangent_modulus_1, tangent_modulus_2);
        p_elem_prop->SetValue(CONSTITUTIVE_LAW, p_test_law);

        std::vector<ModelPart::IndexType> element_nodes {1,2};
        auto p_element = r_model_part.CreateNewElement("TrussElement3D2N", 1, element_nodes, p_elem_prop);
        p_element->Initialize(r_model_part.GetProcessInfo());

        // Set the displacements for testing
        p_bottom_node->FastGetSolutionStepValue(DISPLACEMENT, 0) = array_1d<double, 3>{0.0, 0.0, 0.0};
        p_top_node->FastGetSolutionStepValue(DISPLACEMENT, 0) = array_1d<double, 3>{0.0, 0.0, elongation};

        auto p_truss_element = dynamic_cast<TrussElement3D2N*>(p_element.get());
        KRATOS_EXPECT_NE(p_truss_element, nullptr);
        KRATOS_EXPECT_NEAR(tangent_modulus_2, p_truss_element->ReturnTangentModulus1D(r_model_part.GetProcessInfo()), 1.0e-8);
    }

    KRATOS_TEST_CASE_IN_SUITE(TangentModulusOfTrussElementLinear3D2NUsesLinearStrain, KratosStructuralMechanicsFastSuite)
    {
        Model current_model;
        auto& r_model_part = CreateTestModelPart(current_model);

        // Create two nodes and a truss element
        constexpr auto length = 2.0;
        auto p_bottom_node = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        auto p_top_node = r_model_part.CreateNewNode(2, 0.0, 0.0, length);
        AddDisplacementDofsElement(r_model_part);

        auto p_elem_prop = r_model_part.CreateNewProperties(0);
        constexpr auto elongation = 0.01;
        constexpr auto linear_strain = elongation / length;
        constexpr auto tangent_modulus_1 = 2.0e+03;
        constexpr auto tangent_modulus_2 = 1.0e+03;
        auto p_test_law = std::make_shared<StubBilinearLaw>(1.0001 * linear_strain, tangent_modulus_1, tangent_modulus_2);
        p_elem_prop->SetValue(CONSTITUTIVE_LAW, p_test_law);

        std::vector<ModelPart::IndexType> element_nodes {1,2};
        auto p_element = r_model_part.CreateNewElement("TrussLinearElement3D2N", 1, element_nodes, p_elem_prop);
        p_element->Initialize(r_model_part.GetProcessInfo());

        // Set the displacements for testing
        p_bottom_node->FastGetSolutionStepValue(DISPLACEMENT, 0) = array_1d<double, 3>{0.0, 0.0, 0.0};
        p_top_node->FastGetSolutionStepValue(DISPLACEMENT, 0) = array_1d<double, 3>{0.0, 0.0, elongation};

        auto p_truss_element = dynamic_cast<TrussElementLinear3D2N*>(p_element.get());
        KRATOS_EXPECT_NE(p_truss_element, nullptr);
        KRATOS_EXPECT_NEAR(tangent_modulus_1, p_truss_element->ReturnTangentModulus1D(r_model_part.GetProcessInfo()), 1.0e-10);
    }

}
}
