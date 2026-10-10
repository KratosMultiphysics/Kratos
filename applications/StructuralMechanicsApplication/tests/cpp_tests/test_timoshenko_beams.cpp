// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Gennady Markelov
//

// System includes
#include <utility>

// Project includes
#include "containers/model.h"
#include "structural_mechanics_fast_suite.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos::Testing
{
void FillModelPartWithVariablesNodesAndDoF(ModelPart& rModelPart, std::size_t number_of_nodes, double XEnd, double YEnd)
{
    rModelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);
    rModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
    rModelPart.AddNodalSolutionStepVariable(ROTATION_Z);

    // Create the test geometry
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, XEnd, YEnd, 0.0);
    if (number_of_nodes==3) {
        rModelPart.CreateNewNode(3, XEnd/2, YEnd/2, 0.0);
    }

    for (auto& r_node : rModelPart.Nodes()){
        r_node.AddDof(DISPLACEMENT_X);
        r_node.AddDof(DISPLACEMENT_Y);
        r_node.AddDof(ROTATION_Z);
    }
}

std::vector<ModelPart::IndexType> GetElementNodesFromModelPart(ModelPart& rModelPart)
{
    std::vector<ModelPart::IndexType> element_node_ids;
    const auto& r_nodes = rModelPart.Nodes();
    element_node_ids.reserve(r_nodes.size());
    std::transform(r_nodes.ptr_begin(), r_nodes.ptr_end(), std::back_inserter(element_node_ids),
                   [](const auto& rNodePtr){ return rNodePtr->Id(); });
    return element_node_ids;
}

template<SizeType TNNodes>
void Create2DBeamModel_and_CheckPK2Stress(const std::string & TimoshenkoBeamElementName, const std::vector<double>& rExpectedShearStress, const std::vector<double>& rExpectedBendingMoment)
{
    Model current_model;
    auto &r_model_part = current_model.CreateModelPart("ModelPart",1);
    constexpr double directional_length = 2.0;
    FillModelPartWithVariablesNodesAndDoF(r_model_part, TNNodes, directional_length, directional_length);

    // Set the element properties
    auto p_elem_prop = r_model_part.CreateNewProperties(0);
    constexpr auto youngs_modulus = 2.0e+06;
    p_elem_prop->SetValue(YOUNG_MODULUS, youngs_modulus);
    p_elem_prop->SetValue(CROSS_AREA, 1.0);
    p_elem_prop->SetValue(I33, 1.0);
    p_elem_prop->SetValue(AREA_EFFECTIVE_Y, 5.0/6.0);

    const auto &r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("TimoshenkoBeamElasticConstitutiveLaw");
    p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

    auto element_node_ids = GetElementNodesFromModelPart(r_model_part);
    auto p_element = r_model_part.CreateNewElement(TimoshenkoBeamElementName, 1, element_node_ids, p_elem_prop);
    const auto& r_process_info = r_model_part.GetProcessInfo();
    p_element->Initialize(r_process_info); // Initialize the element to initialize the constitutive law

    //constexpr auto induced_strain = 0.1;
    constexpr auto induced_strain = 0.1;
    p_element->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT) += ScalarVector(3, induced_strain * directional_length);
    p_element->GetGeometry()[1].FastGetSolutionStepValue(ROTATION_Z) += 0.1;
    if constexpr(TNNodes==3) {
        p_element->GetGeometry()[2].FastGetSolutionStepValue(DISPLACEMENT) += ScalarVector(3, induced_strain * directional_length/2.0);
        p_element->GetGeometry()[2].FastGetSolutionStepValue(ROTATION_Z) += 0.1/2.0;
    }

    std::vector<Vector> stress_vectors;
    p_element->CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, stress_vectors, r_process_info);

    constexpr auto expected_stress = induced_strain * youngs_modulus;
    constexpr auto tolerance = 1.0e-5;
    Vector expected_stress_vector(3);
    expected_stress_vector[0] = expected_stress; expected_stress_vector[1] = rExpectedBendingMoment[0]; expected_stress_vector[2] = rExpectedShearStress[0];
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_stress_vector, stress_vectors[0], tolerance);
    expected_stress_vector[0] = expected_stress; expected_stress_vector[1] = rExpectedBendingMoment[1]; expected_stress_vector[2] = rExpectedShearStress[1];
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_stress_vector, stress_vectors[1], tolerance);
    if (stress_vectors.size()>2) {
        expected_stress_vector[0] = expected_stress; expected_stress_vector[1] = rExpectedBendingMoment[2]; expected_stress_vector[2] = rExpectedShearStress[2];
        KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_stress_vector, stress_vectors[2], tolerance);
    }

    expected_stress_vector[0] = expected_stress; expected_stress_vector[1] = rExpectedBendingMoment[0]; expected_stress_vector[2] = rExpectedShearStress[0];
    Vector pre_stress(3);
    pre_stress[0] = 1.0e5; pre_stress[1] = 1.0e4; pre_stress[2] = 1.0e3;
    p_element->GetProperties().SetValue(BEAM_PRESTRESS_PK2, pre_stress);
    p_element->CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, stress_vectors, r_process_info);
    expected_stress_vector += pre_stress;
    KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(expected_stress_vector, stress_vectors[0], tolerance);
}

template<SizeType TNNodes>
void Create2DPlaneStrainBeamModel_and_CheckPK2Stress(const std::string& TimoshenkoBeamElementName)
{
    Model current_model;
    auto &r_model_part = current_model.CreateModelPart("ModelPart",1);
    constexpr double directional_length = 2.0;
    FillModelPartWithVariablesNodesAndDoF(r_model_part, TNNodes, directional_length, 0.0);

    // Set the element properties
    auto p_elem_prop = r_model_part.CreateNewProperties(0);
    constexpr auto youngs_modulus = 2.0e+06;
    constexpr auto poissons_ratio = 0.2;
    constexpr auto thickness = 1.0;
    p_elem_prop->SetValue(YOUNG_MODULUS, youngs_modulus);
    p_elem_prop->SetValue(POISSON_RATIO, poissons_ratio);
    p_elem_prop->SetValue(THICKNESS, thickness);
    constexpr auto effective_shear_thickness = 5.0 / 6.0;
    p_elem_prop->SetValue(THICKNESS_EFFECTIVE_Y, effective_shear_thickness);

    const auto &r_clone_cl = KratosComponents<ConstitutiveLaw>::Get("TimoshenkoBeamPlaneStrainElasticConstitutiveLaw");
    p_elem_prop->SetValue(CONSTITUTIVE_LAW, r_clone_cl.Clone());

    // Create the test element
    auto element_node_ids = GetElementNodesFromModelPart(r_model_part);
    auto p_element = r_model_part.CreateNewElement(TimoshenkoBeamElementName, 1, element_node_ids, p_elem_prop);
    const auto& r_process_info = r_model_part.GetProcessInfo();
    p_element->Initialize(r_process_info); // Initialize the element to initialize the constitutive law
    constexpr auto induced_strain = 0.2;
    constexpr auto induced_rotation = 0.1;
    constexpr auto transverse_bending_displacement = induced_rotation * directional_length / 2.0;
    constexpr auto shear_displacement = 0.03;
    Vector end_displacement(3);
    end_displacement[0] = induced_strain * directional_length; end_displacement[1] = transverse_bending_displacement + shear_displacement; end_displacement[2] = 0.0;
    p_element->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT) += end_displacement;
    p_element->GetGeometry()[1].FastGetSolutionStepValue(ROTATION_Z) += induced_rotation;
    constexpr auto mid_induced_rotation = induced_rotation / 2.0;
    constexpr auto mid_transverse_bending_displacement = mid_induced_rotation * directional_length / 4.0;
    Vector mid_displacement(3);
    mid_displacement[0] = induced_strain * directional_length / 2.0; mid_displacement[1] = mid_transverse_bending_displacement + shear_displacement / 2.0; mid_displacement[2] = 0.0;
    p_element->GetGeometry()[2].FastGetSolutionStepValue(DISPLACEMENT) += mid_displacement;
    p_element->GetGeometry()[2].FastGetSolutionStepValue(ROTATION_Z) += mid_induced_rotation;

    std::vector<Vector> stress_vectors;
    p_element->CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, stress_vectors, r_process_info);
    KRATOS_INFO("Hallo 3") << std::endl;

    constexpr auto expected_distributed_normal_force = induced_strain * youngs_modulus * thickness / (1.0 - poissons_ratio *poissons_ratio);
    constexpr auto expected_distributed_moment =  induced_rotation * youngs_modulus * (thickness * thickness * thickness / 12.0) / (directional_length * ( 1.0 - poissons_ratio *poissons_ratio));
    constexpr auto expected_distributed_shear_force = -(youngs_modulus * effective_shear_thickness / (2.0 * (1.0 + poissons_ratio) )) * (shear_displacement / directional_length);
    constexpr auto tolerance = 1.0e-5;
    Vector expected_stress_vector(5);
    expected_stress_vector[0] = expected_distributed_normal_force; expected_stress_vector[1] = expected_distributed_moment; expected_stress_vector[2] = expected_distributed_shear_force; expected_stress_vector[3] = poissons_ratio*expected_distributed_normal_force; expected_stress_vector[4] = poissons_ratio*expected_distributed_moment;
    KRATOS_EXPECT_VECTOR_NEAR(expected_stress_vector, stress_vectors[0], tolerance);
    KRATOS_EXPECT_VECTOR_NEAR(expected_stress_vector, stress_vectors[1], tolerance);

    Vector pre_stress(5);
    pre_stress[0] = 1.0e5; pre_stress[1] = 1.0e4; pre_stress[2] = 1.0e3; pre_stress[3] = 2.e4; pre_stress[4] = 2.e3;
    p_element->GetProperties().SetValue(BEAM_PRESTRESS_PK2, pre_stress);
    p_element->CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, stress_vectors, r_process_info);
    expected_stress_vector += pre_stress;
    KRATOS_EXPECT_VECTOR_NEAR(expected_stress_vector, stress_vectors[0], tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(LinearTimoshenkoBeam2D2N_CalculatesPK2Stress, KratosStructuralMechanicsFastSuite)
{
    const std::vector expected_shear_stress{-32608.7, -32608.7, -32608.7};
    const std::vector expected_bending_moment{34989.6, 70710.7, 106432.0};

    Create2DBeamModel_and_CheckPK2Stress<2>("LinearTimoshenkoBeamElement2D2N", expected_shear_stress, expected_bending_moment);
}

KRATOS_TEST_CASE_IN_SUITE(LinearTimoshenkodBeam2D3N_CalculatesPK2Stress, KratosStructuralMechanicsFastSuite)
{
    const std::vector expected_shear_stress{2604.74, -21673.6, -48084.9};
    const std::vector expected_bending_moment{86720.4, 68516.6, 61527.5};
    Create2DBeamModel_and_CheckPK2Stress<3>("LinearTimoshenkoBeamElement2D3N", expected_shear_stress, expected_bending_moment);
}

KRATOS_TEST_CASE_IN_SUITE(LinearTimoshenkodCurvedBeam2D3N_CalculatesPK2Stress, KratosStructuralMechanicsFastSuite)
{
    const std::vector expected_shear_stress{17610.4, 65722.9};
    const std::vector expected_bending_moment{70710.7, 70710.7};
    Create2DBeamModel_and_CheckPK2Stress<3>("LinearTimoshenkoCurvedBeamElement2D3N", expected_shear_stress, expected_bending_moment);
}

KRATOS_TEST_CASE_IN_SUITE(LinearTimoshenkodCurvedBeam2D3N_CalculatesPK2StressPlaneStrain, KratosStructuralMechanicsFastSuite)
{
    Create2DPlaneStrainBeamModel_and_CheckPK2Stress<3>("LinearTimoshenkoCurvedBeamElement2D3N");
}

class ParametrizedFinalizeSolutionStepForTimoshenkoBeams
    : public ::testing::TestWithParam<std::tuple<std::size_t, std::string>>
{
public:
    void SetUp() override
    {
        mpStructuralMechanicsApp = std::make_shared<KratosStructuralMechanicsApplication>();
        mpStructuralMechanicsApp->Register();
    }
private:
    KratosStructuralMechanicsApplication::Pointer mpStructuralMechanicsApp;
};

TEST_P(ParametrizedFinalizeSolutionStepForTimoshenkoBeams, FinalizeSolutionStepIsCalledForTimoshenkoBeams)
{
    class MockConstitutiveLaw : public ConstitutiveLaw
    {
    public:
        [[nodiscard]] ConstitutiveLaw::Pointer Clone() const override {return std::make_shared<MockConstitutiveLaw>();}
        [[nodiscard]] SizeType GetStrainSize() const override {return 6;} // should be 3 for 2D, but for those elements it seems no problem
        [[nodiscard]] bool RequiresFinalizeMaterialResponse() override { return true; }
        [[nodiscard]] StressMeasure GetStressMeasure() override {return StressMeasure_PK2;}

        MOCK_METHOD(void, FinalizeMaterialResponsePK2, (Parameters&), (override));
    };

    // Arrange
    const auto& [number_of_nodes, element_type] = GetParam();

    Model current_model;
    auto& r_model_part = current_model.CreateModelPart("ModelPart",2);
    constexpr auto directional_length = 2.0;
    FillModelPartWithVariablesNodesAndDoF(r_model_part, number_of_nodes, directional_length, directional_length);

    // Set the element properties
    auto p_elem_prop = r_model_part.CreateNewProperties(0);
    p_elem_prop->SetValue(YOUNG_MODULUS, 2.0e+06);
    p_elem_prop->SetValue(CROSS_AREA, 1.0);
    p_elem_prop->SetValue(I33, 1.0);
    p_elem_prop->SetValue(I22, 1.0);
    p_elem_prop->SetValue(AREA_EFFECTIVE_Y, 5.0/6.0);
    p_elem_prop->SetValue(AREA_EFFECTIVE_Z, 5.0/6.0);

    // mock CL that counts calls to FinalizeMaterialResponsePK2
    auto p_mockconstitutivelaw = std::make_shared<MockConstitutiveLaw>();
    p_elem_prop->SetValue(CONSTITUTIVE_LAW, p_mockconstitutivelaw);

    auto element_node_ids = GetElementNodesFromModelPart(r_model_part);
    auto p_element = r_model_part.CreateNewElement(element_type, 1, element_node_ids, p_elem_prop);

    const auto& r_process_info = r_model_part.GetProcessInfo();
    p_element->Initialize(r_process_info); // Initialize the element to initialize the constitutive law
    std::vector<ConstitutiveLaw::Pointer> constitutive_laws;
    p_element->CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, constitutive_laws, r_process_info);
    // the deed
    for( auto& rp_constitutive_law : constitutive_laws )
    {
        auto p_mock_law = dynamic_cast<MockConstitutiveLaw*>(rp_constitutive_law.get());
        EXPECT_CALL(*p_mock_law, FinalizeMaterialResponsePK2).Times(1);
    }
    p_element->FinalizeSolutionStep(r_process_info);
}

INSTANTIATE_TEST_SUITE_P(
    KratosStructuralMechanicsFastSuite,
    ParametrizedFinalizeSolutionStepForTimoshenkoBeams,
    ::testing::Values(
        std::make_tuple(std::size_t{2}, "LinearTimoshenkoBeamElement2D2N"),
        std::make_tuple(std::size_t{2}, "LinearTimoshenkoBeamElement3D2N"),
        std::make_tuple(std::size_t{3}, "LinearTimoshenkoBeamElement2D3N"),
        std::make_tuple(std::size_t{3}, "LinearTimoshenkoCurvedBeamElement2D3N"),
        std::make_tuple(std::size_t{3}, "LinearTimoshenkoCurvedBeamElement3D3N")
    )
);

class ParametrizedInternalExternalForcesForTimoshenkoBeams
    : public ::testing::TestWithParam<std::tuple<std::size_t, std::size_t, std::string>>
{
public:
    void SetUp() override
    {
        mpStructuralMechanicsApp = std::make_shared<KratosStructuralMechanicsApplication>();
        mpStructuralMechanicsApp->Register();
    }

protected:
    static constexpr double mDensity   = 2.0;
    static constexpr double mCrossArea = 0.5;

    // Creates an inclined straight beam of length 5 (2D) or 7 (3D), loaded by a body force
    Element::Pointer CreateBeamElement(ModelPart& rModelPart)
    {
        const auto& [dimension, number_of_nodes, element_type] = GetParam();

        rModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
        rModelPart.AddNodalSolutionStepVariable(ROTATION);

        const auto end = (dimension == 2) ? array_1d<double, 3>{3.0, 4.0, 0.0} : array_1d<double, 3>{2.0, 3.0, 6.0};
        rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
        rModelPart.CreateNewNode(2, end[0], end[1], end[2]);
        if (number_of_nodes == 3) {
            rModelPart.CreateNewNode(3, end[0] / 2.0, end[1] / 2.0, end[2] / 2.0);
        }

        auto p_elem_prop = rModelPart.CreateNewProperties(0);
        p_elem_prop->SetValue(YOUNG_MODULUS, 2.0e+06);
        p_elem_prop->SetValue(POISSON_RATIO, 0.2);
        p_elem_prop->SetValue(CROSS_AREA, mCrossArea);
        p_elem_prop->SetValue(I33, 0.01);
        p_elem_prop->SetValue(I22, 0.02);
        p_elem_prop->SetValue(IT, 0.03);
        p_elem_prop->SetValue(AREA_EFFECTIVE_Y, 0.4);
        p_elem_prop->SetValue(AREA_EFFECTIVE_Z, 0.4);
        p_elem_prop->SetValue(DENSITY, mDensity);
        p_elem_prop->SetValue(VOLUME_ACCELERATION, GetVolumeAcceleration());

        const auto cl_name = (dimension == 2) ? "TimoshenkoBeamElasticConstitutiveLaw" : "TimoshenkoBeamElasticConstitutiveLaw3D";
        p_elem_prop->SetValue(CONSTITUTIVE_LAW, KratosComponents<ConstitutiveLaw>::Get(cl_name).Clone());

        auto p_element = rModelPart.CreateNewElement(element_type, 1, GetElementNodesFromModelPart(rModelPart), p_elem_prop);
        p_element->Initialize(rModelPart.GetProcessInfo());
        return p_element;
    }

    array_1d<double, 3> GetVolumeAcceleration() const
    {
        const auto dimension = std::get<0>(GetParam());
        return (dimension == 2) ? array_1d<double, 3>{1.0, -9.81, 0.0} : array_1d<double, 3>{1.0, -9.81, 0.5};
    }

    static void ApplyArbitraryDeformation(ModelPart& rModelPart)
    {
        for (auto& r_node : rModelPart.Nodes()) {
            const double i = static_cast<double>(r_node.Id());
            r_node.FastGetSolutionStepValue(DISPLACEMENT) = array_1d<double, 3>{0.01 * i, -0.02 * i + 0.005, 0.003 * i * i};
            r_node.FastGetSolutionStepValue(ROTATION)     = array_1d<double, 3>{0.01 * i, -0.005 * i * i, 0.02 * i};
        }
    }

    // Nodal displacements and rotations ordered as the element DoFs
    Vector GetElementDofValues(const Element& rElement) const
    {
        const auto dimension = std::get<0>(GetParam());
        const auto& r_geometry = rElement.GetGeometry();
        const std::size_t dofs_per_node = (dimension == 2) ? 3 : 6;
        Vector values(dofs_per_node * r_geometry.size());
        for (std::size_t i_node = 0; i_node < r_geometry.size(); ++i_node) {
            const auto& r_displacement = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
            const auto& r_rotation     = r_geometry[i_node].FastGetSolutionStepValue(ROTATION);
            const std::size_t index = i_node * dofs_per_node;
            if (dimension == 2) {
                values[index]     = r_displacement[0];
                values[index + 1] = r_displacement[1];
                values[index + 2] = r_rotation[2];
            } else {
                for (std::size_t i = 0; i < 3; ++i) {
                    values[index + i]     = r_displacement[i];
                    values[index + 3 + i] = r_rotation[i];
                }
            }
        }
        return values;
    }

private:
    KratosStructuralMechanicsApplication::Pointer mpStructuralMechanicsApp;
};

TEST_P(ParametrizedInternalExternalForcesForTimoshenkoBeams, ExternalMinusInternalForcesEqualsRightHandSide)
{
    Model current_model;
    auto& r_model_part = current_model.CreateModelPart("ModelPart", 1);
    auto p_element = CreateBeamElement(r_model_part);
    ApplyArbitraryDeformation(r_model_part);
    const auto& r_process_info = r_model_part.GetProcessInfo();

    Vector internal_forces, external_forces, rhs;
    p_element->Calculate(INTERNAL_FORCES_VECTOR, internal_forces, r_process_info);
    p_element->Calculate(EXTERNAL_FORCES_VECTOR, external_forces, r_process_info);
    p_element->CalculateRightHandSide(rhs, r_process_info);

    KRATOS_EXPECT_GT(norm_2(internal_forces), 0.0);
    KRATOS_EXPECT_GT(norm_2(external_forces), 0.0);
    KRATOS_EXPECT_VECTOR_NEAR(rhs, Vector(external_forces - internal_forces), 1.0e-8);
}

TEST_P(ParametrizedInternalExternalForcesForTimoshenkoBeams, InternalForcesEqualStiffnessTimesDisplacements)
{
    Model current_model;
    auto& r_model_part = current_model.CreateModelPart("ModelPart", 1);
    auto p_element = CreateBeamElement(r_model_part);
    ApplyArbitraryDeformation(r_model_part);
    const auto& r_process_info = r_model_part.GetProcessInfo();

    Vector internal_forces;
    Matrix lhs;
    p_element->Calculate(INTERNAL_FORCES_VECTOR, internal_forces, r_process_info);
    p_element->CalculateLeftHandSide(lhs, r_process_info);

    // Linear elastic law, hence the internal forces are K * u
    const Vector expected_internal_forces = prod(lhs, GetElementDofValues(*p_element));
    KRATOS_EXPECT_GT(norm_2(expected_internal_forces), 0.0);
    KRATOS_EXPECT_VECTOR_NEAR(internal_forces, expected_internal_forces, 1.0e-10 * norm_2(expected_internal_forces));
}

TEST_P(ParametrizedInternalExternalForcesForTimoshenkoBeams, ExternalForcesSumUpToTotalBodyForce)
{
    Model current_model;
    auto& r_model_part = current_model.CreateModelPart("ModelPart", 1);
    auto p_element = CreateBeamElement(r_model_part);
    const auto& r_process_info = r_model_part.GetProcessInfo();

    Vector external_forces;
    p_element->Calculate(EXTERNAL_FORCES_VECTOR, external_forces, r_process_info);

    const auto& r_geometry = p_element->GetGeometry();
    const auto dimension = std::get<0>(GetParam());
    const std::size_t dofs_per_node = external_forces.size() / r_geometry.size();
    Vector resultant = ZeroVector(dimension);
    for (std::size_t i_node = 0; i_node < r_geometry.size(); ++i_node) {
        for (std::size_t i_dim = 0; i_dim < dimension; ++i_dim) {
            resultant[i_dim] += external_forces[i_node * dofs_per_node + i_dim];
        }
    }

    const double length = norm_2(r_geometry[1].Coordinates() - r_geometry[0].Coordinates());
    Vector expected_resultant(dimension);
    for (std::size_t i_dim = 0; i_dim < dimension; ++i_dim) {
        expected_resultant[i_dim] = mDensity * mCrossArea * length * GetVolumeAcceleration()[i_dim];
    }
    KRATOS_EXPECT_VECTOR_NEAR(resultant, expected_resultant, 1.0e-10);
}


TEST_P(ParametrizedInternalExternalForcesForTimoshenkoBeams, CalculateThrowsForUnsupportedVectorVariable)
{
    Model current_model;
    auto& r_model_part = current_model.CreateModelPart("ModelPart", 1);
    auto p_element = CreateBeamElement(r_model_part);

    Vector output;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(p_element->Calculate(PK2_STRESS_VECTOR, output, r_model_part.GetProcessInfo()),
                                      "Variable PK2_STRESS_VECTOR not supported in element");
}

INSTANTIATE_TEST_SUITE_P(
    KratosStructuralMechanicsFastSuite,
    ParametrizedInternalExternalForcesForTimoshenkoBeams,
    ::testing::Values(
        std::make_tuple(std::size_t{2}, std::size_t{2}, "LinearTimoshenkoBeamElement2D2N"),
        std::make_tuple(std::size_t{2}, std::size_t{3}, "LinearTimoshenkoBeamElement2D3N"),
        std::make_tuple(std::size_t{2}, std::size_t{3}, "LinearTimoshenkoCurvedBeamElement2D3N"),
        std::make_tuple(std::size_t{3}, std::size_t{2}, "LinearTimoshenkoBeamElement3D2N"),
        std::make_tuple(std::size_t{3}, std::size_t{3}, "LinearTimoshenkoCurvedBeamElement3D3N")
    )
);

Element::Pointer CreateLinearTimoshenkoBeam3D2N(ModelPart& rModelPart, const array_1d<double, 3>& rEnd)
{
    rModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
    rModelPart.AddNodalSolutionStepVariable(ROTATION);
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, rEnd[0], rEnd[1], rEnd[2]);

    auto p_elem_prop = rModelPart.CreateNewProperties(0);
    p_elem_prop->SetValue(YOUNG_MODULUS, 2.0e+06);
    p_elem_prop->SetValue(POISSON_RATIO, 0.2);
    p_elem_prop->SetValue(CROSS_AREA, 0.5);
    p_elem_prop->SetValue(I33, 0.01);
    p_elem_prop->SetValue(I22, 0.02);
    p_elem_prop->SetValue(IT, 0.03);
    p_elem_prop->SetValue(AREA_EFFECTIVE_Y, 0.4);
    p_elem_prop->SetValue(AREA_EFFECTIVE_Z, 0.3);
    p_elem_prop->SetValue(CONSTITUTIVE_LAW, KratosComponents<ConstitutiveLaw>::Get("TimoshenkoBeamElasticConstitutiveLaw3D").Clone());

    auto p_element = rModelPart.CreateNewElement("LinearTimoshenkoBeamElement3D2N", 1, std::vector<ModelPart::IndexType>{1, 2}, p_elem_prop);
    p_element->Initialize(rModelPart.GetProcessInfo());
    return p_element;
}

// Applies the same deformation w.r.t. the local axes, which are the columns of rRotation
void ApplyLocalDeformation(Element& rBeam, const BoundedMatrix<double, 3, 3>& rRotation)
{
    const std::vector local_displacements{array_1d<double, 3>{0.01, -0.02, 0.03}, array_1d<double, 3>{0.04, 0.05, -0.06}};
    const std::vector local_rotations{array_1d<double, 3>{0.002, -0.003, 0.004}, array_1d<double, 3>{-0.005, 0.006, 0.007}};
    for (std::size_t i_node = 0; i_node < 2; ++i_node) {
        auto& r_node = rBeam.GetGeometry()[i_node];
        r_node.FastGetSolutionStepValue(DISPLACEMENT) = prod(rRotation, local_displacements[i_node]);
        r_node.FastGetSolutionStepValue(ROTATION)     = prod(rRotation, local_rotations[i_node]);
    }
}

KRATOS_TEST_CASE_IN_SUITE(LinearTimoshenkoBeam3D2N_InternalForcesOfBeamAlongGlobalX, KratosStructuralMechanicsFastSuite)
{
    // Along the global X axis the local axes are the global axes
    Model current_model;
    auto& r_model_part = current_model.CreateModelPart("Aligned", 1);
    auto p_beam = CreateLinearTimoshenkoBeam3D2N(r_model_part, array_1d<double, 3>{7.0, 0.0, 0.0});
    ApplyLocalDeformation(*p_beam, IdentityMatrix(3));

    Vector internal_forces;
    p_beam->Calculate(INTERNAL_FORCES_VECTOR, internal_forces, r_model_part.GetProcessInfo());

    KRATOS_EXPECT_EQ(internal_forces.size(), 12);
    KRATOS_EXPECT_GT(norm_2(internal_forces), 0.0);

    // Axial force N = EA * (u2x - u1x) / L = 2e6 * 0.5 * 0.03 / 7
    const double axial_force = 3.0e+04 / 7.0;
    KRATOS_EXPECT_NEAR(internal_forces[0], -axial_force, 1.0e-10 * axial_force);
    KRATOS_EXPECT_NEAR(internal_forces[6],  axial_force, 1.0e-10 * axial_force);

    // Torque T = G * IT * (rot2x - rot1x) / L = 2e6 / (2 * 1.2) * 0.03 * (-0.007) / 7
    const double torque = -25.0;
    KRATOS_EXPECT_NEAR(internal_forces[3], -torque, 1.0e-10 * std::abs(torque));
    KRATOS_EXPECT_NEAR(internal_forces[9],  torque, 1.0e-10 * std::abs(torque));
}

KRATOS_TEST_CASE_IN_SUITE(LinearTimoshenkoBeam3D2N_InternalForcesAreFrameInvariant, KratosStructuralMechanicsFastSuite)
{
    // The same beam, deformed in the same way w.r.t. its local axes, is created along the
    // global X axis (local axes == global axes) and along an inclined axis. The internal forces
    // of the inclined beam must be the ones of the aligned beam, rotated to the inclined axes.
    Model current_model;
    auto& r_aligned_model_part  = current_model.CreateModelPart("Aligned", 1);
    auto& r_inclined_model_part = current_model.CreateModelPart("Inclined", 1);
    auto p_aligned_beam  = CreateLinearTimoshenkoBeam3D2N(r_aligned_model_part,  array_1d<double, 3>{7.0, 0.0, 0.0});
    auto p_inclined_beam = CreateLinearTimoshenkoBeam3D2N(r_inclined_model_part, array_1d<double, 3>{2.0, 3.0, 6.0});

    // Columns of the rotation are the local axes of the inclined beam: the beam axis, the default
    // local axis 2 (global Y) made orthogonal to the beam axis, and the cross product of both
    const double sqrt_10 = std::sqrt(10.0);
    BoundedMatrix<double, 3, 3> rotation;
    rotation(0, 0) = 2.0 / 7.0; rotation(0, 1) = -3.0 / (7.0 * sqrt_10); rotation(0, 2) = -3.0 / sqrt_10;
    rotation(1, 0) = 3.0 / 7.0; rotation(1, 1) = 20.0 / (7.0 * sqrt_10); rotation(1, 2) = 0.0;
    rotation(2, 0) = 6.0 / 7.0; rotation(2, 1) = -9.0 / (7.0 * sqrt_10); rotation(2, 2) = 1.0 / sqrt_10;
    ApplyLocalDeformation(*p_aligned_beam, IdentityMatrix(3));
    ApplyLocalDeformation(*p_inclined_beam, rotation);

    Vector aligned_internal_forces;
    Vector inclined_internal_forces;
    p_aligned_beam->Calculate(INTERNAL_FORCES_VECTOR, aligned_internal_forces, r_aligned_model_part.GetProcessInfo());
    p_inclined_beam->Calculate(INTERNAL_FORCES_VECTOR, inclined_internal_forces, r_inclined_model_part.GetProcessInfo());

    // Rotate each (force or moment) block of the aligned beam to the inclined axes
    Vector expected_internal_forces(12);
    for (std::size_t block = 0; block < 4; ++block) {
        array_1d<double, 3> aligned_block;
        for (std::size_t i = 0; i < 3; ++i) {
            aligned_block[i] = aligned_internal_forces[3 * block + i];
        }
        const array_1d<double, 3> rotated_block = prod(rotation, aligned_block);
        for (std::size_t i = 0; i < 3; ++i) {
            expected_internal_forces[3 * block + i] = rotated_block[i];
        }
    }

    KRATOS_EXPECT_VECTOR_NEAR(inclined_internal_forces, expected_internal_forces, 1.0e-10 * norm_2(expected_internal_forces));
}

}
