// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                     license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo & Sergio Jimenez
//

// System includes
#include <limits>
#include <string>
#include <typeinfo>

// External includes

// Project includes
#include "includes/kratos_components.h"
#include "includes/process_info.h"
#include "includes/stream_serializer.h"
#include "containers/model.h"

// Application includes
#include "constitutive_laws_application_variables.h"
#include "tests/cpp_tests/constitutive_laws_fast_suite.h"

// Integrator
#include "custom_constitutive/auxiliary_files/cl_integrators/generic_cl_integrator_plasticity.h"
#include "custom_constitutive/auxiliary_files/cl_integrators/generic_cl_integrator_damage.h"

// Yield surfaces
#include "custom_constitutive/auxiliary_files/yield_surfaces/generic_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/von_mises_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/modified_mohr_coulomb_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/rankine_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/simo_ju_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/drucker_prager_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/tresca_yield_surface.h"

// Plastic potentials
#include "custom_constitutive/auxiliary_files/plastic_potentials/generic_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/von_mises_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/tresca_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/modified_mohr_coulomb_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/drucker_prager_plastic_potential.h"

// Constitutive law
#include "custom_constitutive/composites/rule_of_mixtures_law.h"
#include "custom_constitutive/small_strains/plastic_damage/generic_small_strain_plastic_damage_model.h"
#include "includes/model_part.h"
#include "geometries/tetrahedra_3d_4.h"

namespace Kratos::Testing
{
namespace
{

using PlasticDamageLaw = GenericSmallStrainPlasticDamageModel<
    GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>,
    GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>;

constexpr char plastic_damage_law_name[] = "SmallStrainPlasticDamageModel3DVonMisesVonMisesVonMises";
constexpr double elastic_strain = 0.005;
constexpr double initial_strain = 0.006;
constexpr double continuation_strain = 0.007;
constexpr double advancing_strain = 0.008;
constexpr double damage_dissipation_saturation = 0.99999;
constexpr double stress_tolerance = 1.0e-6;

void set_plastic_damage_properties(Properties& rProperties)
{
    rProperties.SetValue(DENSITY, 7850.0);
    rProperties.SetValue(YOUNG_MODULUS, 2.01e11);
    rProperties.SetValue(POISSON_RATIO, 0.3);
    rProperties.SetValue(FRACTURE_ENERGY, 50e5);
    rProperties.SetValue(FRACTURE_ENERGY_DAMAGE_PROCESS, 50e5);
    rProperties.SetValue(YIELD_STRESS_COMPRESSION, 8.389e8);
    rProperties.SetValue(YIELD_STRESS_TENSION, 8.389e8);
    rProperties.SetValue(HARDENING_CURVE, 0);
    rProperties.SetValue(SOFTENING_TYPE, 0);
}

class PlasticDamageTestSetup
{
public:
    PlasticDamageTestSetup()
        : mGeometry(
            Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0),
            Kratos::make_intrusive<Node>(2, 1.0, 0.0, 0.0),
            Kratos::make_intrusive<Node>(3, 0.0, 1.0, 0.0),
            Kratos::make_intrusive<Node>(4, 0.0, 0.0, 1.0))
    {
    }

    void Initialize(ConstitutiveLaw& rLaw, const Properties& rProperties) const
    {
        rLaw.InitializeMaterial(rProperties, mGeometry, ScalarVector(4, 0.25));
        KRATOS_EXPECT_EQ(rLaw.Check(rProperties, mGeometry, mProcessInfo), 0);
    }

    Vector CalculateAndFinalize(
        ConstitutiveLaw& rLaw,
        const Properties& rProperties,
        const double Strain) const
    {
        Vector strain = ZeroVector(6);
        strain[0] = Strain;
        Vector stress = ZeroVector(6);
        Matrix tangent = ZeroMatrix(6, 6);
        ConstitutiveLaw::Parameters values(mGeometry, rProperties, mProcessInfo);
        values.SetStrainVector(strain);
        values.SetStressVector(stress);
        values.SetConstitutiveMatrix(tangent);
        values.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        values.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
        rLaw.InitializeMaterialResponseCauchy(values);
        rLaw.CalculateMaterialResponseCauchy(values);
        const Vector response = stress;
        rLaw.FinalizeMaterialResponseCauchy(values);
        return response;
    }

private:
    Tetrahedra3D4<Node> mGeometry;
    ProcessInfo mProcessInfo;
};

double get_damage_dissipation(const ConstitutiveLaw& rLaw)
{
    // The public state variables omit this accepted history. Read just its named trace field.
    StreamSerializer serializer(Serializer::SERIALIZER_TRACE_ERROR);
    serializer.pGetBuffer()->precision(std::numeric_limits<double>::max_digits10);
    serializer.save("Law", rLaw);
    const std::string trace = serializer.GetStringRepresentation();
    const std::string tag = "\"DamageDissipation\"\n";
    const auto tag_position = trace.find(tag);
    KRATOS_ERROR_IF(tag_position == std::string::npos) << "Missing DamageDissipation in the law's trace.";
    return std::stod(trace.substr(tag_position + tag.size()));
}

} // namespace

/**
    * Check the correct calculation of the integrated stress with the CL's
    */

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawPlasticDamageInternalVariables, KratosConstitutiveLawsFastSuite)
{
    //
    // Test: check correct behavior of internal and calculated variables
    //

    typedef GenericSmallStrainPlasticDamageModel<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> PD_CL;

    Model current_model;
    ModelPart& r_test_model_part = current_model.CreateModelPart("Main");
    PD_CL cl = PD_CL();

    KRATOS_EXPECT_TRUE(cl.Has(DAMAGE));  // = True
    KRATOS_EXPECT_TRUE(cl.Has(UNIAXIAL_STRESS));  // = True
    KRATOS_EXPECT_TRUE(cl.Has(PLASTIC_DISSIPATION));  // = True
    KRATOS_EXPECT_TRUE(cl.Has(PLASTIC_STRAIN_VECTOR));  // = True
    KRATOS_EXPECT_TRUE(cl.Has(INTERNAL_VARIABLES));  // = True

    Vector internal_variables_w(9);
    internal_variables_w[0] = 0.0;
    internal_variables_w[1] = 0.1;
    internal_variables_w[2] = 0.2;
    internal_variables_w[3] = 0.3;
    internal_variables_w[4] = 0.4;
    internal_variables_w[5] = 0.5;
    internal_variables_w[6] = 0.6;
    internal_variables_w[7] = 0.7;
    internal_variables_w[8] = 0.8;
    cl.SetValue(INTERNAL_VARIABLES, internal_variables_w, r_test_model_part.GetProcessInfo());
    Vector internal_variables_r;  // CL should internally resize it to 6
    cl.GetValue(INTERNAL_VARIABLES, internal_variables_r);

    KRATOS_EXPECT_NEAR(internal_variables_r.size(), 9., 1.e-5);  // = True
    KRATOS_EXPECT_NEAR(internal_variables_r[0], 0.0, 1.e-5);  // = True
    KRATOS_EXPECT_NEAR(internal_variables_r[1], 0.1, 1.e-5);  // = True
    KRATOS_EXPECT_NEAR(internal_variables_r[2], 0.2, 1.e-5);  // = True
    KRATOS_EXPECT_NEAR(internal_variables_r[3], 0.3, 1.e-5);  // = True
    KRATOS_EXPECT_NEAR(internal_variables_r[4], 0.4, 1.e-5);  // = True
    KRATOS_EXPECT_NEAR(internal_variables_r[5], 0.5, 1.e-5);  // = True
    KRATOS_EXPECT_NEAR(internal_variables_r[6], 0.6, 1.e-5);  // = True
    KRATOS_EXPECT_NEAR(internal_variables_r[7], 0.7, 1.e-5);  // = True
    KRATOS_EXPECT_NEAR(internal_variables_r[8], 0.8, 1.e-5);  // = True
}

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawPlasticDamageExponential, KratosConstitutiveLawsFastSuite)
{
    typedef GenericSmallStrainPlasticDamageModel<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>, GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> PD_CL;

    ConstitutiveLaw::Parameters cl_parameters;
    Properties material_properties;
    Vector stress_vector, strain_vector;

    Model current_model;

    ModelPart& test_model_part = current_model.CreateModelPart("Main");

    Node::Pointer p_node_1 = test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    Node::Pointer p_node_2 = test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    Node::Pointer p_node_3 = test_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    Node::Pointer p_node_4 = test_model_part.CreateNewNode(4, 0.0, 0.0, 1.0);

    Tetrahedra3D4<Node> Geom = Tetrahedra3D4<Node>(p_node_1, p_node_2, p_node_3, p_node_4);

    stress_vector = ZeroVector(6);

    strain_vector = ZeroVector(6);
    strain_vector[0] = 0.01;
    strain_vector[1] = 0.0;
    strain_vector[2] = 0.0;
    strain_vector[3] = 0.0;
    strain_vector[4] = 0.0;
    strain_vector[5] = 0.0;

    material_properties.SetValue(DENSITY, 7850.0);
    material_properties.SetValue(YOUNG_MODULUS, 2.01e11);
    material_properties.SetValue(POISSON_RATIO, 0.3);
    material_properties.SetValue(FRACTURE_ENERGY, 50e5);
    material_properties.SetValue(FRACTURE_ENERGY_DAMAGE_PROCESS, 50e5);
    material_properties.SetValue(YIELD_STRESS_COMPRESSION, 8.389e8);
    material_properties.SetValue(YIELD_STRESS_TENSION, 8.389e8);
    material_properties.SetValue(HARDENING_CURVE, 0);
    material_properties.SetValue(SOFTENING_TYPE, 0);

    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions = cl_parameters.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

    cl_parameters.SetElementGeometry(Geom);
    cl_parameters.SetProcessInfo(test_model_part.GetProcessInfo());
    cl_parameters.SetMaterialProperties(material_properties);
    cl_parameters.SetStrainVector(strain_vector);
    cl_parameters.SetStressVector(stress_vector);
    Matrix const_matrix;
    cl_parameters.SetConstitutiveMatrix(const_matrix);

    // Create the CL's
    PD_CL PlasticDamage_CL = PD_CL();

    std::vector<double> expected_value;
    expected_value = {1.60206e+09,8.81364e+08,8.81364e+08,0,0,0};

    Vector TestStress;
    Vector dummy;
    PlasticDamage_CL.InitializeMaterial(material_properties, Geom, dummy);
    PlasticDamage_CL.CalculateMaterialResponseCauchy(cl_parameters);
    PlasticDamage_CL.FinalizeMaterialResponseCauchy(cl_parameters);
    TestStress = cl_parameters.GetStressVector();

    //Check the results
    for (int comp = 0; comp < 6; comp++) {
        KRATOS_EXPECT_NEAR(expected_value[comp], TestStress[comp], 0.0001e+08);
    }
}

KRATOS_TEST_CASE_IN_SUITE(SmallStrainPlasticDamage3DClonePreservesDamageDissipation, KratosConstitutiveLawsFastSuite)
{
    const PlasticDamageTestSetup setup;
    Properties properties;
    set_plastic_damage_properties(properties);
    const auto& r_registered_law = KratosComponents<ConstitutiveLaw>::Get(plastic_damage_law_name);
    const auto p_original = r_registered_law.Clone();
    setup.Initialize(*p_original, properties);
    setup.CalculateAndFinalize(*p_original, properties, elastic_strain);
    KRATOS_EXPECT_EQ(get_damage_dissipation(*p_original), 0.0);
    setup.CalculateAndFinalize(*p_original, properties, initial_strain);
    const double original_dissipation = get_damage_dissipation(*p_original);
    KRATOS_EXPECT_GT(original_dissipation, 0.0);
    KRATOS_EXPECT_LT(original_dissipation, damage_dissipation_saturation);

    const auto p_clone = p_original->Clone();
    KRATOS_EXPECT_NE(p_clone.get(), p_original.get());
    KRATOS_EXPECT_NE(dynamic_cast<const PlasticDamageLaw*>(p_clone.get()), nullptr);
    KRATOS_EXPECT_EQ(get_damage_dissipation(*p_clone), original_dissipation);

    const Vector original_stress = setup.CalculateAndFinalize(*p_original, properties, continuation_strain);
    const Vector cloned_stress = setup.CalculateAndFinalize(*p_clone, properties, continuation_strain);
    KRATOS_EXPECT_VECTOR_NEAR(cloned_stress, original_stress, stress_tolerance);
    const double continued_dissipation = get_damage_dissipation(*p_original);
    KRATOS_EXPECT_GT(continued_dissipation, original_dissipation);
    KRATOS_EXPECT_NEAR(get_damage_dissipation(*p_clone), continued_dissipation, 1.0e-12);

    setup.CalculateAndFinalize(*p_clone, properties, advancing_strain);
    KRATOS_EXPECT_GT(get_damage_dissipation(*p_clone), continued_dissipation);
    KRATOS_EXPECT_EQ(get_damage_dissipation(*p_original), continued_dissipation);
}

KRATOS_TEST_CASE_IN_SUITE(ParallelRuleOfMixturesClonePreservesPlasticDamageDissipation, KratosConstitutiveLawsFastSuite)
{
    const auto& r_plastic_damage_law = KratosComponents<ConstitutiveLaw>::Get(plastic_damage_law_name);
    const auto& r_isotropic_damage_law = KratosComponents<ConstitutiveLaw>::Get("SmallStrainIsotropicDamage3DVonMises");
    const auto set_layer_properties = [&](Properties& rMixtureProperties) {
        auto p_plastic_damage_properties = Kratos::make_shared<Properties>(11);
        set_plastic_damage_properties(*p_plastic_damage_properties);
        p_plastic_damage_properties->SetValue(CONSTITUTIVE_LAW, r_plastic_damage_law.Clone());
        rMixtureProperties.AddSubProperties(p_plastic_damage_properties);

        auto p_isotropic_damage_properties = Kratos::make_shared<Properties>(*p_plastic_damage_properties);
        p_isotropic_damage_properties->SetId(12);
        p_isotropic_damage_properties->SetValue(CONSTITUTIVE_LAW, r_isotropic_damage_law.Clone());
        rMixtureProperties.AddSubProperties(p_isotropic_damage_properties);
    };

    Properties source_properties(1);
    Properties clone_properties(1);
    set_layer_properties(source_properties);
    set_layer_properties(clone_properties);
    const PlasticDamageTestSetup source_setup;
    const PlasticDamageTestSetup clone_setup;

    const auto& r_mixture_law = KratosComponents<ConstitutiveLaw>::Get("ParallelRuleOfMixturesLaw3D");
    const auto p_source = r_mixture_law.Create(Parameters(R"({"combination_factors": [1.0, 3.0]})"));
    source_setup.Initialize(*p_source, source_properties);
    source_setup.CalculateAndFinalize(*p_source, source_properties, elastic_strain);
    source_setup.CalculateAndFinalize(*p_source, source_properties, initial_strain);

    auto* p_source_mixture = dynamic_cast<ParallelRuleOfMixturesLaw<3>*>(p_source.get());
    ASSERT_NE(p_source_mixture, nullptr);
    const auto& r_source_layers = p_source_mixture->GetConstitutiveLaws();
    ASSERT_EQ(r_source_layers.size(), 2);
    ASSERT_NE(r_source_layers[0], nullptr);
    ASSERT_NE(r_source_layers[1], nullptr);
    const auto& r_plastic_damage_layer = *r_source_layers[0];
    const auto& r_isotropic_damage_layer = *r_source_layers[1];
    ASSERT_TRUE(typeid(r_plastic_damage_layer) == typeid(r_plastic_damage_law));
    ASSERT_TRUE(typeid(r_isotropic_damage_layer) == typeid(r_isotropic_damage_law));
    const double accepted_dissipation = get_damage_dissipation(r_plastic_damage_layer);
    ASSERT_GT(accepted_dissipation, 0.0);
    ASSERT_LT(accepted_dissipation, damage_dissipation_saturation);

    const auto p_clone = p_source->Clone();
    auto* p_clone_mixture = dynamic_cast<ParallelRuleOfMixturesLaw<3>*>(p_clone.get());
    ASSERT_NE(p_clone_mixture, nullptr);
    const auto& r_cloned_layers = p_clone_mixture->GetConstitutiveLaws();
    ASSERT_EQ(r_cloned_layers.size(), r_source_layers.size());
    for (std::size_t i = 0; i < r_source_layers.size(); ++i) {
        ASSERT_NE(r_cloned_layers[i], nullptr);
        KRATOS_EXPECT_NE(r_cloned_layers[i].get(), r_source_layers[i].get());
        const auto& r_cloned_layer = *r_cloned_layers[i];
        const auto& r_source_layer = *r_source_layers[i];
        KRATOS_EXPECT_TRUE(typeid(r_cloned_layer) == typeid(r_source_layer));
    }
    const std::vector<double> expected_factors{0.25, 0.75};
    KRATOS_EXPECT_EQ(p_source_mixture->GetCombinationFactors(), expected_factors);
    KRATOS_EXPECT_EQ(p_clone_mixture->GetCombinationFactors(), expected_factors);
    const double cloned_dissipation = get_damage_dissipation(*r_cloned_layers[0]);
    KRATOS_EXPECT_DOUBLE_EQ(cloned_dissipation, accepted_dissipation);

    const Vector source_stress = source_setup.CalculateAndFinalize(*p_source, source_properties, continuation_strain);
    const Vector clone_stress = clone_setup.CalculateAndFinalize(*p_clone, clone_properties, continuation_strain);
    KRATOS_EXPECT_VECTOR_NEAR(clone_stress, source_stress, stress_tolerance);
    const double continued_source_dissipation = get_damage_dissipation(*r_source_layers[0]);
    const double continued_clone_dissipation = get_damage_dissipation(*r_cloned_layers[0]);
    KRATOS_EXPECT_DOUBLE_EQ(continued_clone_dissipation, continued_source_dissipation);

    clone_setup.CalculateAndFinalize(*p_clone, clone_properties, advancing_strain);
    const double advanced_clone_dissipation = get_damage_dissipation(*r_cloned_layers[0]);
    const double unchanged_source_dissipation = get_damage_dissipation(*r_source_layers[0]);
    KRATOS_EXPECT_GT(advanced_clone_dissipation, continued_clone_dissipation);
    KRATOS_EXPECT_DOUBLE_EQ(unchanged_source_dissipation, continued_source_dissipation);
}

} // namespace Kratos::Testing
