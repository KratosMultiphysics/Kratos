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

void CreateTrussModel2N_and_CheckPK2Stress(std::string TrussElementName)
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
        auto p_element = r_model_part.CreateNewElement(std::move(TrussElementName), 1, element_nodes, p_elem_prop);
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

void CreateTrussModel_and_CheckInternalAndExternalForces(std::string TrussElementName,
                                                         std::size_t NumberOfNodes,
                                                         const array_1d<double, 3>& rUnitDirection,
                                                         const array_1d<double, 3>& rVolumeAcceleration)
    {
        Model current_model;
        auto& r_model_part = CreateTestModelPart(current_model);

        // Set the element properties
        constexpr auto youngs_modulus = 2.0e+06;
        constexpr auto area           = 0.01;
        constexpr auto density        = 7850.0;
        constexpr auto pre_stress     = 1.0e+03;
        auto p_elem_prop = r_model_part.CreateNewProperties(0);
        p_elem_prop->SetValue(YOUNG_MODULUS, youngs_modulus);
        p_elem_prop->SetValue(CROSS_AREA, area);
        p_elem_prop->SetValue(DENSITY, density);
        p_elem_prop->SetValue(TRUSS_PRESTRESS_PK2, pre_stress);
        p_elem_prop->SetValue(VOLUME_ACCELERATION, rVolumeAcceleration);
        p_elem_prop->SetValue(CONSTITUTIVE_LAW, KratosComponents<ConstitutiveLaw>::Get("TrussConstitutiveLaw").Clone());

        // Create the test element, uniformly stretched along its axis (the third node, if any, is the mid node)
        constexpr auto length     = 3.0;
        constexpr auto elongation = 0.03;
        const auto relative_node_positions = NumberOfNodes == 2 ? std::vector<double>{0.0, 1.0}
                                                                : std::vector<double>{0.0, 1.0, 0.5};
        std::vector<ModelPart::IndexType> element_nodes;
        for (std::size_t i = 0; i < NumberOfNodes; ++i) {
            const array_1d<double, 3> coordinates = relative_node_positions[i] * length * rUnitDirection;
            auto p_node = r_model_part.CreateNewNode(i + 1, coordinates[0], coordinates[1], coordinates[2]);
            p_node->FastGetSolutionStepValue(DISPLACEMENT) = relative_node_positions[i] * elongation * rUnitDirection;
            element_nodes.push_back(p_node->Id());
        }

        auto p_element = r_model_part.CreateNewElement(std::move(TrussElementName), 1, element_nodes, p_elem_prop);
        const auto& r_process_info = r_model_part.GetProcessInfo();
        p_element->Initialize(r_process_info); // Initialize the element to initialize the constitutive law

        Vector internal_forces, external_forces, rhs;
        p_element->Calculate(INTERNAL_FORCES_VECTOR, internal_forces, r_process_info);
        p_element->Calculate(EXTERNAL_FORCES_VECTOR, external_forces, r_process_info);
        p_element->CalculateRightHandSide(rhs, r_process_info);

        // The internal forces are the axial force acting along the axis on the end nodes, the external forces
        // are the self weight distributed over the nodes by the (linear or quadratic) shape functions
        const auto axial_force  = (youngs_modulus * elongation / length + pre_stress) * area;
        const auto total_mass   = density * area * length;
        const auto axial_force_signs = NumberOfNodes == 2 ? std::vector<double>{-1.0, 1.0}
                                                          : std::vector<double>{-1.0, 1.0, 0.0};
        const auto mass_fractions    = NumberOfNodes == 2 ? std::vector<double>{0.5, 0.5}
                                                          : std::vector<double>{1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0};

        const auto dimension = p_element->GetGeometry().WorkingSpaceDimension();
        Vector expected_internal_forces(NumberOfNodes * dimension), expected_external_forces(NumberOfNodes * dimension);
        for (std::size_t i = 0; i < NumberOfNodes; ++i) {
            for (std::size_t j = 0; j < dimension; ++j) {
                expected_internal_forces[i * dimension + j] = axial_force_signs[i] * axial_force * rUnitDirection[j];
                expected_external_forces[i * dimension + j] = mass_fractions[i] * total_mass * rVolumeAcceleration[j];
            }
        }

        constexpr auto tolerance = 1.0e-8;
        KRATOS_EXPECT_VECTOR_NEAR(internal_forces, expected_internal_forces, tolerance);
        KRATOS_EXPECT_VECTOR_NEAR(external_forces, expected_external_forces, tolerance);

        // The residual must be consistent with the separately computed forces
        const Vector expected_rhs = external_forces - internal_forces;
        KRATOS_EXPECT_VECTOR_NEAR(rhs, expected_rhs, tolerance);
    }

// Data of the prestressed TrussElement3D2N-based test trusses
constexpr auto truss_3D2N_youngs_modulus = 2.0e+06;
constexpr auto truss_3D2N_area           = 0.01;
constexpr auto truss_3D2N_pre_stress     = 1.0e+03;
constexpr auto truss_3D2N_length         = 3.0;

void CreateTruss3D2NModel_and_CheckInternalAndExternalForces(std::string TrussElementName,
                                                             double Elongation,
                                                             double ExpectedAxialForce)
    {
        Model current_model;
        auto& r_model_part = CreateTestModelPart(current_model);
        r_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);

        // Set the element properties
        constexpr auto density = 7850.0;
        auto p_elem_prop = r_model_part.CreateNewProperties(0);
        p_elem_prop->SetValue(YOUNG_MODULUS, truss_3D2N_youngs_modulus);
        p_elem_prop->SetValue(CROSS_AREA, truss_3D2N_area);
        p_elem_prop->SetValue(DENSITY, density);
        p_elem_prop->SetValue(TRUSS_PRESTRESS_PK2, truss_3D2N_pre_stress);
        p_elem_prop->SetValue(CONSTITUTIVE_LAW, KratosComponents<ConstitutiveLaw>::Get("TrussConstitutiveLaw").Clone());

        // Create the test element along an inclined axis, loaded by self weight and stretched along its axis
        const array_1d<double, 3> unit_direction{1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0};
        const array_1d<double, 3> volume_acceleration{1.0, 2.0, -9.81};
        const array_1d<double, 3> end_coordinates = truss_3D2N_length * unit_direction;
        r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        auto p_end_node = r_model_part.CreateNewNode(2, end_coordinates[0], end_coordinates[1], end_coordinates[2]);
        for (auto& r_node : r_model_part.Nodes()) {
            r_node.FastGetSolutionStepValue(VOLUME_ACCELERATION) = volume_acceleration;
        }
        p_end_node->FastGetSolutionStepValue(DISPLACEMENT) = Elongation * unit_direction;

        const std::vector<ModelPart::IndexType> element_nodes {1, 2};
        auto p_element = r_model_part.CreateNewElement(std::move(TrussElementName), 1, element_nodes, p_elem_prop);
        const auto& r_process_info = r_model_part.GetProcessInfo();
        p_element->Initialize(r_process_info); // Initialize the element to initialize the constitutive law

        Vector internal_forces, external_forces, rhs;
        p_element->Calculate(INTERNAL_FORCES_VECTOR, internal_forces, r_process_info);
        p_element->Calculate(EXTERNAL_FORCES_VECTOR, external_forces, r_process_info);
        p_element->CalculateRightHandSide(rhs, r_process_info);

        // The internal forces are the axial force acting along the axis on both nodes,
        // the external forces are the self weight lumped equally into both nodes
        const auto nodal_mass = 0.5 * density * truss_3D2N_area * truss_3D2N_length;
        Vector expected_internal_forces(6), expected_external_forces(6);
        for (std::size_t j = 0; j < 3; ++j) {
            expected_internal_forces[j]     = -ExpectedAxialForce * unit_direction[j];
            expected_internal_forces[3 + j] =  ExpectedAxialForce * unit_direction[j];
            expected_external_forces[j]     = nodal_mass * volume_acceleration[j];
            expected_external_forces[3 + j] = nodal_mass * volume_acceleration[j];
        }

        constexpr auto tolerance = 1.0e-8;
        KRATOS_EXPECT_VECTOR_NEAR(internal_forces, expected_internal_forces, tolerance);
        KRATOS_EXPECT_VECTOR_NEAR(external_forces, expected_external_forces, tolerance);

        // The residual must be consistent with the separately computed forces
        const Vector expected_rhs = external_forces - internal_forces;
        KRATOS_EXPECT_VECTOR_NEAR(rhs, expected_rhs, tolerance);
    }

double CalculateExpectedAxialForceOfTrussElement3D2N(double Elongation)
    {
        // The PK2 stress follows from the Green-Lagrange strain, the axial force is pushed forward to the current configuration
        const auto current_length = truss_3D2N_length + Elongation;
        const auto green_lagrange_strain = (current_length * current_length - truss_3D2N_length * truss_3D2N_length) /
                                           (2.0 * truss_3D2N_length * truss_3D2N_length);
        return (truss_3D2N_youngs_modulus * green_lagrange_strain + truss_3D2N_pre_stress) * truss_3D2N_area * current_length / truss_3D2N_length;
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
        p_top_node->Coordinates()[2] += elongation; // to account for Green-Lagrange strain

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
        p_top_node->Coordinates()[2] += elongation;

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
        CreateTrussModel2N_and_CheckPK2Stress("TrussLinearElement3D2N");
    }

    KRATOS_TEST_CASE_IN_SUITE(LinearTrussElement2D2N_CalculatesPK2Stress, KratosStructuralMechanicsFastSuite)
    {
        CreateTrussModel2N_and_CheckPK2Stress("LinearTrussElement2D2N");
    }

    KRATOS_TEST_CASE_IN_SUITE(LinearTrussElement3D2N_CalculatesPK2Stress, KratosStructuralMechanicsFastSuite)
    {
        CreateTrussModel2N_and_CheckPK2Stress("LinearTrussElement3D2N");
    }

    KRATOS_TEST_CASE_IN_SUITE(LinearTrussElement2D2N_CalculatesInternalAndExternalForces, KratosStructuralMechanicsFastSuite)
    {
        CreateTrussModel_and_CheckInternalAndExternalForces("LinearTrussElement2D2N", 2,
            array_1d<double, 3>{0.6, 0.8, 0.0}, array_1d<double, 3>{1.0, -9.81, 0.0});
    }

    KRATOS_TEST_CASE_IN_SUITE(LinearTrussElement2D3N_CalculatesInternalAndExternalForces, KratosStructuralMechanicsFastSuite)
    {
        CreateTrussModel_and_CheckInternalAndExternalForces("LinearTrussElement2D3N", 3,
            array_1d<double, 3>{0.6, 0.8, 0.0}, array_1d<double, 3>{1.0, -9.81, 0.0});
    }

    KRATOS_TEST_CASE_IN_SUITE(LinearTrussElement3D2N_CalculatesInternalAndExternalForces, KratosStructuralMechanicsFastSuite)
    {
        CreateTrussModel_and_CheckInternalAndExternalForces("LinearTrussElement3D2N", 2,
            array_1d<double, 3>{1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0}, array_1d<double, 3>{1.0, 2.0, -9.81});
    }

    KRATOS_TEST_CASE_IN_SUITE(LinearTrussElement3D3N_CalculatesInternalAndExternalForces, KratosStructuralMechanicsFastSuite)
    {
        CreateTrussModel_and_CheckInternalAndExternalForces("LinearTrussElement3D3N", 3,
            array_1d<double, 3>{1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0}, array_1d<double, 3>{1.0, 2.0, -9.81});
    }

    KRATOS_TEST_CASE_IN_SUITE(LinearTrussElement3D2N_CalculateThrowsForUnsupportedVectorVariable, KratosStructuralMechanicsFastSuite)
    {
        Model current_model;
        auto& r_model_part = CreateTestModelPart(current_model);
        auto [p_bottom_node, p_top_node] = CreateEndNodes(r_model_part, 2.0);
        auto p_elem_prop = r_model_part.CreateNewProperties(0);
        const std::vector<ModelPart::IndexType> element_nodes {p_bottom_node->Id(), p_top_node->Id()};
        auto p_element = r_model_part.CreateNewElement("LinearTrussElement3D2N", 1, element_nodes, p_elem_prop);

        Vector output;
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Calculate(PK2_STRESS_VECTOR, output, r_model_part.GetProcessInfo()),
                                          "Variable PK2_STRESS_VECTOR is not supported in element")
    }

    KRATOS_TEST_CASE_IN_SUITE(TrussElement3D2N_CalculatesInternalAndExternalForces, KratosStructuralMechanicsFastSuite)
    {
        constexpr auto elongation = 0.03;
        CreateTruss3D2NModel_and_CheckInternalAndExternalForces("TrussElement3D2N", elongation,
            CalculateExpectedAxialForceOfTrussElement3D2N(elongation));
    }

    KRATOS_TEST_CASE_IN_SUITE(TrussElementLinear3D2N_CalculatesInternalAndExternalForces, KratosStructuralMechanicsFastSuite)
    {
        constexpr auto elongation = 0.03;
        constexpr auto expected_axial_force =
            (truss_3D2N_youngs_modulus * elongation / truss_3D2N_length + truss_3D2N_pre_stress) * truss_3D2N_area;
        CreateTruss3D2NModel_and_CheckInternalAndExternalForces("TrussLinearElement3D2N", elongation, expected_axial_force);
    }

    KRATOS_TEST_CASE_IN_SUITE(CableElement3D2N_CalculatesInternalAndExternalForcesInTension, KratosStructuralMechanicsFastSuite)
    {
        constexpr auto elongation = 0.03;
        CreateTruss3D2NModel_and_CheckInternalAndExternalForces("CableElement3D2N", elongation,
            CalculateExpectedAxialForceOfTrussElement3D2N(elongation));
    }

    KRATOS_TEST_CASE_IN_SUITE(CableElement3D2N_HasNoInternalForcesInCompression, KratosStructuralMechanicsFastSuite)
    {
        constexpr auto elongation = -0.03;
        constexpr auto expected_axial_force = 0.0; // a compressed cable does not carry any load
        CreateTruss3D2NModel_and_CheckInternalAndExternalForces("CableElement3D2N", elongation, expected_axial_force);
    }

    KRATOS_TEST_CASE_IN_SUITE(TrussElement3D2N_CalculateThrowsForUnsupportedVectorVariable, KratosStructuralMechanicsFastSuite)
    {
        Model current_model;
        auto& r_model_part = CreateTestModelPart(current_model);
        auto [p_bottom_node, p_top_node] = CreateEndNodes(r_model_part, 2.0);
        auto p_elem_prop = r_model_part.CreateNewProperties(0);
        const std::vector<ModelPart::IndexType> element_nodes {p_bottom_node->Id(), p_top_node->Id()};
        auto p_element = r_model_part.CreateNewElement("TrussElement3D2N", 1, element_nodes, p_elem_prop);

        Vector output;
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Calculate(PK2_STRESS_VECTOR, output, r_model_part.GetProcessInfo()),
                                          "Variable PK2_STRESS_VECTOR is not supported in element")
    }
}
