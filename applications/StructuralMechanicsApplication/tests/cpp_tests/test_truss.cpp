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
#include "structural_mechanics_fast_suite.h"
#include "structural_mechanics_application_variables.h"

#include "custom_elements/truss_elements/truss_element_3D2N.hpp"
#include "custom_elements/truss_elements/truss_element_linear_3D2N.hpp"

#include <utility>

namespace
{

using namespace Kratos;

class StubBilinearLaw : public ConstitutiveLaw
{
public:
    // Only implement the interface that is needed by the tests
    StubBilinearLaw(double Strain, double TangentModulus1, double TangentModulus2) :
        mStrain{Strain}, mTangentModuli{TangentModulus1, TangentModulus2}
    {}

    StubBilinearLaw() = default;

    [[nodiscard]] ConstitutiveLaw::Pointer Clone() const override
    {
        return std::make_shared<StubBilinearLaw>(*this);
    }

    [[nodiscard]] SizeType GetStrainSize() const override
    {
        return 1;
    }

    double& CalculateValue(Parameters& rParameterValues, const Variable<double>& rThisVariable, double& rValue) override
    {
        KRATOS_ERROR_IF_NOT(rThisVariable == TANGENT_MODULUS);

        rValue = rParameterValues.GetStrainVector()[0] < mStrain ? mTangentModuli[0] : mTangentModuli[1];
        return rValue;
    }

    using ConstitutiveLaw::CalculateValue;

    void InitializeMaterial(const Properties&,
                            const GeometryType&,
                            const Vector&) override
    {
        mIsInitialized = true;
    }

    [[nodiscard]] bool IsInitialized() const
    {
        return mIsInitialized;
    }

private:
    double mStrain = 0.0;
    array_1d<double, 2> mTangentModuli{2.0, 1.0};
    bool mIsInitialized = false;
};


ModelPart& CreateTestModelPart(Model& rModel)
{
  auto&r_result = rModel.CreateModelPart("ModelPart", 1);
  r_result.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);
  r_result.AddNodalSolutionStepVariable(DISPLACEMENT);
  return r_result;
}

std::pair<Node::Pointer, Node::Pointer> CreateEndNodes(ModelPart& rModelPart, double VerticalDistance)
{
  auto p_bottom_node = rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
  auto p_top_node = rModelPart.CreateNewNode(2, 0.0, 0.0, VerticalDistance);
  return std::make_pair(p_bottom_node, p_top_node);
}

std::shared_ptr<StubBilinearLaw> CreateStubBilinearLaw(double Elongation, double Length, double TangentModulus1, double TangentModulus2)
{
    const auto linear_strain = Elongation / Length;
    const auto new_length = Length + Elongation;
    const auto green_lagrange_strain = (new_length * new_length - Length * Length) / (2.0 * Length * Length);
    const auto threshold_strain = 0.5 * (linear_strain + green_lagrange_strain);
    return std::make_shared<StubBilinearLaw>(threshold_strain, TangentModulus1, TangentModulus2);
}

}


namespace Kratos::Testing
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

        constexpr auto length = 2.0;
        auto [p_bottom_node, p_top_node] = CreateEndNodes(r_model_part, length);
        AddDisplacementDofsElement(r_model_part);

        constexpr auto elongation = 0.01;
        constexpr auto tangent_modulus_1 = 2.0e+03;
        constexpr auto tangent_modulus_2 = 1.0e+03;
        auto p_elem_prop = r_model_part.CreateNewProperties(0);
        p_elem_prop->SetValue(CONSTITUTIVE_LAW, CreateStubBilinearLaw(elongation, length, tangent_modulus_1, tangent_modulus_2));

        std::vector<ModelPart::IndexType> element_nodes {p_bottom_node->Id(), p_top_node->Id()};
        auto p_element = r_model_part.CreateNewElement("TrussElement3D2N", 1, element_nodes, p_elem_prop);
        p_element->Initialize(r_model_part.GetProcessInfo());

        p_bottom_node->FastGetSolutionStepValue(DISPLACEMENT) = array_1d<double, 3>{0.0, 0.0, 0.0};
        p_top_node->FastGetSolutionStepValue(DISPLACEMENT) = array_1d<double, 3>{0.0, 0.0, elongation};

        auto p_truss_element = dynamic_cast<TrussElement3D2N*>(p_element.get());
        KRATOS_EXPECT_NE(p_truss_element, nullptr);
        KRATOS_EXPECT_DOUBLE_EQ(tangent_modulus_2, p_truss_element->ReturnTangentModulus1D(r_model_part.GetProcessInfo()));
    }

    KRATOS_TEST_CASE_IN_SUITE(TangentModulusOfTrussElementLinear3D2NUsesLinearStrain, KratosStructuralMechanicsFastSuite)
    {
        Model current_model;
        auto& r_model_part = CreateTestModelPart(current_model);

        constexpr auto length = 2.0;
        auto [p_bottom_node, p_top_node] = CreateEndNodes(r_model_part, length);
        AddDisplacementDofsElement(r_model_part);

        constexpr auto elongation = 0.01;
        constexpr auto tangent_modulus_1 = 2.0e+03;
        constexpr auto tangent_modulus_2 = 1.0e+03;
        auto p_elem_prop = r_model_part.CreateNewProperties(0);
        p_elem_prop->SetValue(CONSTITUTIVE_LAW, CreateStubBilinearLaw(elongation, length, tangent_modulus_1, tangent_modulus_2));

        std::vector<ModelPart::IndexType> element_nodes {p_bottom_node->Id(), p_top_node->Id()};
        auto p_element = r_model_part.CreateNewElement("TrussLinearElement3D2N", 1, element_nodes, p_elem_prop);
        p_element->Initialize(r_model_part.GetProcessInfo());

        p_bottom_node->FastGetSolutionStepValue(DISPLACEMENT) = array_1d<double, 3>{0.0, 0.0, 0.0};
        p_top_node->FastGetSolutionStepValue(DISPLACEMENT) = array_1d<double, 3>{0.0, 0.0, elongation};

        auto p_truss_element = dynamic_cast<TrussElement3D2N*>(p_element.get());
        KRATOS_EXPECT_NE(p_truss_element, nullptr);
        KRATOS_EXPECT_DOUBLE_EQ(tangent_modulus_1, p_truss_element->ReturnTangentModulus1D(r_model_part.GetProcessInfo()));
    }

    KRATOS_TEST_CASE_IN_SUITE(TrussElementLinear3D2NInitializesConstitutiveLaw, KratosStructuralMechanicsFastSuite)
    {
        Model current_model;
        auto& r_model_part = CreateTestModelPart(current_model);
        constexpr auto length = 2.0;
        auto [p_bottom_node, p_top_node] = CreateEndNodes(r_model_part, length);
        AddDisplacementDofsElement(r_model_part);
        auto p_elem_prop = r_model_part.CreateNewProperties(0);
        p_elem_prop->SetValue(CONSTITUTIVE_LAW, std::make_shared<StubBilinearLaw>());
        const std::vector<ModelPart::IndexType> element_nodes {p_bottom_node->Id(), p_top_node->Id()};
        auto p_element = r_model_part.CreateNewElement("TrussLinearElement3D2N", 1, element_nodes, p_elem_prop);

        p_element->Initialize(r_model_part.GetProcessInfo());
        std::vector<ConstitutiveLaw::Pointer> constitutive_laws;
        p_element->CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, constitutive_laws, r_model_part.GetProcessInfo());
        auto p_constitutive_law = dynamic_cast<const StubBilinearLaw*>(constitutive_laws[0].get());
        KRATOS_EXPECT_TRUE(p_constitutive_law->IsInitialized())
    }

    KRATOS_TEST_CASE_IN_SUITE(TrussElementLinear3D2N_CalculatesPK2Stress, KratosStructuralMechanicsFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart",1);
        r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);
        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

        // Set the element properties
        auto p_elem_prop = r_model_part.CreateNewProperties(0);
        constexpr auto youngs_modulus = 2.0e+06;
        p_elem_prop->SetValue(YOUNG_MODULUS, youngs_modulus);
        const auto &r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("TrussConstitutiveLaw");
        p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

        // Create the test element
        constexpr double directional_length = 2.0;
        auto p_node_1 = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        auto p_node_2 = r_model_part.CreateNewNode(2, directional_length, directional_length, directional_length);

        AddDisplacementDofsElement(r_model_part);

        std::vector<ModelPart::IndexType> element_nodes {1,2};
        auto p_element = r_model_part.CreateNewElement("TrussLinearElement3D2N", 1, element_nodes, p_elem_prop);
        const auto& r_process_info = r_model_part.GetProcessInfo();
        p_element->Initialize(r_process_info); // Initialize the element to initialize the constitutive law

        constexpr auto induced_strain = 0.1;
        p_element->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT) += ScalarVector(3, induced_strain * directional_length);

        std::vector<Vector> stress_vector;
        p_element->CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, stress_vector, r_process_info);

        constexpr auto expected_stress = induced_strain * youngs_modulus;
        KRATOS_EXPECT_DOUBLE_EQ(expected_stress, stress_vector[0][0]);

        constexpr auto pre_stress = 1.0e5;
        p_element->GetProperties().SetValue(TRUSS_PRESTRESS_PK2, pre_stress);
        p_element->CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, stress_vector, r_process_info);
        KRATOS_EXPECT_DOUBLE_EQ(expected_stress + pre_stress, stress_vector[0][0]);
    }
KRATOS_TEST_CASE_IN_SUITE(LienarTrussElement2D_CalculatesPK2Stress, KratosStructuralMechanicsFastSuite)
    {
        Model current_model;
        auto &r_model_part = current_model.CreateModelPart("ModelPart",1);
        r_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);
        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

        // Set the element properties
        auto p_elem_prop = r_model_part.CreateNewProperties(0);
        constexpr auto youngs_modulus = 2.0e+06;
        p_elem_prop->SetValue(YOUNG_MODULUS, youngs_modulus);
        const auto &r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("TrussConstitutiveLaw");
        p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

        // Create the test element
        constexpr double directional_length = 2.0;
        auto p_node_1 = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        auto p_node_2 = r_model_part.CreateNewNode(2, directional_length, directional_length, directional_length);

        AddDisplacementDofsElement(r_model_part);

        std::vector<ModelPart::IndexType> element_nodes {1,2};
        auto p_element = r_model_part.CreateNewElement("LinearTrussElement2D2N", 1, element_nodes, p_elem_prop);
        const auto& r_process_info = r_model_part.GetProcessInfo();
        p_element->Initialize(r_process_info); // Initialize the element to initialize the constitutive law

        constexpr auto induced_strain = 0.1;
        p_element->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT) += ScalarVector(3, induced_strain * directional_length);

        std::vector<Vector> stress_vector;
        p_element->CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, stress_vector, r_process_info);

        constexpr auto expected_stress = induced_strain * youngs_modulus;
        KRATOS_EXPECT_DOUBLE_EQ(expected_stress, stress_vector[0][0]);

        constexpr auto pre_stress = 1.0e5;
        p_element->GetProperties().SetValue(TRUSS_PRESTRESS_PK2, pre_stress);
        p_element->CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, stress_vector, r_process_info);
        KRATOS_EXPECT_DOUBLE_EQ(expected_stress + pre_stress, stress_vector[0][0]);
    }
}
