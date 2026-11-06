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
//  Main authors:    Sergio Jimenez/Alejandro Cornejo
//

// System includes
#include <cmath>
#include <limits>
#include <string>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/kratos_components.h"
#include "includes/stream_serializer.h"

// Application includes
#include "constitutive_laws_application_variables.h"
#include "tests/cpp_tests/constitutive_laws_fast_suite.h"

// Integrator
#include "custom_constitutive/auxiliary_files/cl_integrators/generic_cl_integrator_damage.h"

// Yield surfaces
#include "custom_constitutive/auxiliary_files/yield_surfaces/von_mises_yield_surface.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/von_mises_plastic_potential.h"

// Constitutive law
#include "custom_constitutive/small_strains/fatigue/generic_small_strain_high_cycle_fatigue_law.h"
#include "geometries/tetrahedra_3d_4.h"

namespace Kratos::Testing
{

typedef GenericSmallStrainHighCycleFatigueLaw<GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> FatigueLawType;

/**
    * Check the correct calculation of the integrated stress with the CL's
    */
KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawHighCycleFatigueExponential, KratosConstitutiveLawsFastSuite)
{
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
    strain_vector[0] = 5.0e-03;
    strain_vector[1] = 0.0;
    strain_vector[2] = 0.0;
    strain_vector[3] = 0.0;
    strain_vector[4] = 0.0;
    strain_vector[5] = 0.0;

    material_properties.SetValue(DENSITY, 7850.0);
    material_properties.SetValue(YOUNG_MODULUS, 2.01e11);
    material_properties.SetValue(POISSON_RATIO, 0.3);
    material_properties.SetValue(FRACTURE_ENERGY, 20e5);
    material_properties.SetValue(YIELD_STRESS, 838.9e6);
    material_properties.SetValue(SOFTENING_TYPE, 1);
    Vector high_cycle_fatigue_coefficients(7);
    high_cycle_fatigue_coefficients[0]=0.5;
    high_cycle_fatigue_coefficients[1]=0.7;
    high_cycle_fatigue_coefficients[2]=0.5;
    high_cycle_fatigue_coefficients[3]=0.0068;
    high_cycle_fatigue_coefficients[4]=3.35;
    high_cycle_fatigue_coefficients[5]=0.0133;
    high_cycle_fatigue_coefficients[6]=0.0068;
    material_properties.SetValue(HIGH_CYCLE_FATIGUE_COEFFICIENTS, high_cycle_fatigue_coefficients);

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
    FatigueLawType HCF_CL = FatigueLawType(8.8914e-01, 7.6380e+08, 6.9760e+08, 7.6380e+08, -2.2914e+08, 984, 5.3239E-07);

    std::vector<double> expected_value;
    expected_value = {1.07636e+09,4.61299e+08,4.61299e+08,0,0,0};

    Vector TestStress;
    Vector dummy;
    HCF_CL.InitializeMaterial(material_properties, Geom, dummy);
    HCF_CL.CalculateMaterialResponseCauchy(cl_parameters);
    HCF_CL.FinalizeMaterialResponseCauchy(cl_parameters);
    TestStress = cl_parameters.GetStressVector();

    //Check the results
    KRATOS_EXPECT_VECTOR_NEAR(expected_value, TestStress, 0.0001e+08);
}

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawHighCycleFatigueClonePreservesHistory, KratosConstitutiveLawsFastSuite)
{
    const Tetrahedra3D4<Node> geometry(
        Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0),
        Kratos::make_intrusive<Node>(2, 1.0, 0.0, 0.0),
        Kratos::make_intrusive<Node>(3, 0.0, 1.0, 0.0),
        Kratos::make_intrusive<Node>(4, 0.0, 0.0, 1.0));
    Properties material_properties;
    material_properties.SetValue(DENSITY, 7850.0);
    material_properties.SetValue(YOUNG_MODULUS, 2.01e11);
    material_properties.SetValue(POISSON_RATIO, 0.3);
    material_properties.SetValue(FRACTURE_ENERGY, 20.0e6);
    material_properties.SetValue(YIELD_STRESS, 838.9e6);
    material_properties.SetValue(SOFTENING_TYPE, 1);
    constexpr double c_factor = 0.5;
    const array_1d<double, 8> fatigue_coefficients{0.5, 0.7, 0.5, 0.0068, 3.35, 0.0133, 0.0068, c_factor};
    material_properties.SetValue(HIGH_CYCLE_FATIGUE_COEFFICIENTS, Vector(fatigue_coefficients));

    ProcessInfo process_info;
    process_info[CURRENT_LOAD_TYPE] = false;
    process_info[NEW_MODEL_PART] = false;
    process_info[ADVANCE_STRATEGY_APPLIED] = false;
    Vector strain = ZeroVector(6);
    Vector stress = ZeroVector(6);
    Matrix tangent;
    ConstitutiveLaw::Parameters values(geometry, material_properties, process_info);
    values.SetStrainVector(strain);
    values.SetStressVector(stress);
    values.SetConstitutiveMatrix(tangent);
    values.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    values.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

    const auto advance = [&](ConstitutiveLaw& rLaw, double AxialStrain) {
        strain[0] = AxialStrain;
        rLaw.InitializeMaterialResponseCauchy(values);
        rLaw.CalculateMaterialResponseCauchy(values);
        const Vector response = stress;
        rLaw.FinalizeMaterialResponseCauchy(values);
        return response;
    };
    const auto read_scalar = [](ConstitutiveLaw& rLaw, const Variable<double>& rVariable) {
        double result = 0.0;
        return rLaw.GetValue(rVariable, result);
    };
    const auto read_cycles = [](ConstitutiveLaw& rLaw) {
        int result = 0;
        return rLaw.GetValue(NUMBER_OF_CYCLES, result);
    };
    const auto read_private_history = [](const ConstitutiveLaw& rLaw, const std::string& rName) {
        StreamSerializer serializer(Serializer::SERIALIZER_TRACE_ERROR);
        serializer.pGetBuffer()->precision(std::numeric_limits<double>::max_digits10);
        serializer.save("Law", rLaw);
        const std::string trace = serializer.GetStringRepresentation();
        const std::string tag = "\"" + rName + "\"\n";
        const auto position = trace.find(tag);
        KRATOS_ERROR_IF(position == std::string::npos) << "Missing fatigue history in the trace: " << rName;
        return std::stod(trace.substr(position + tag.size()));
    };

    const auto p_source = KratosComponents<ConstitutiveLaw>::Get("SmallStrainHighCycleFatigue3DLawVonMises").Clone();
    p_source->InitializeMaterial(material_properties, geometry, ScalarVector(4, 0.25));
    KRATOS_EXPECT_EQ(p_source->Check(material_properties, geometry, process_info), 0);
    advance(*p_source, 0.006);
    advance(*p_source, 0.002);
    const double damage = read_scalar(*p_source, DAMAGE);
    ASSERT_GT(damage, 0.0);
    ASSERT_LT(damage, 1.0);

    // Record the overload damage as the reference state for a new cyclic load block.
    constexpr double valley = -0.0045;
    constexpr double peak = 0.005;
    process_info[CURRENT_LOAD_TYPE] = true;
    process_info[NEW_MODEL_PART] = true;
    advance(*p_source, valley);
    process_info[NEW_MODEL_PART] = false;
    advance(*p_source, peak);
    advance(*p_source, valley);
    advance(*p_source, peak);

    const double cycles_to_failure = read_scalar(*p_source, CYCLES_TO_FAILURE);
    ASSERT_TRUE(std::isfinite(cycles_to_failure));
    ASSERT_GT(cycles_to_failure, 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(read_private_history(*p_source, "ReferenceDamage"), damage);
    KRATOS_EXPECT_DOUBLE_EQ(read_private_history(*p_source, "CFactor"), c_factor);

    // The time-advance process writes this state through SetValue.
    constexpr double previous_cycle_time = 2.0;
    constexpr double cycle_period = 0.25;
    p_source->SetValue(PREVIOUS_CYCLE, previous_cycle_time, process_info);
    p_source->SetValue(CYCLE_PERIOD, cycle_period, process_info);
    p_source->SetValue(PREVIOUS_CYCLE_DAMAGE, damage, process_info);

    const auto p_clone = p_source->Clone();
    KRATOS_EXPECT_NE(p_clone.get(), p_source.get());
    KRATOS_EXPECT_DOUBLE_EQ(read_private_history(*p_clone, "ReferenceDamage"), damage);
    KRATOS_EXPECT_DOUBLE_EQ(read_private_history(*p_clone, "CFactor"), c_factor);
    KRATOS_EXPECT_DOUBLE_EQ(read_scalar(*p_clone, CYCLES_TO_FAILURE), cycles_to_failure);
    KRATOS_EXPECT_DOUBLE_EQ(read_scalar(*p_clone, THRESHOLD_STRESS), read_scalar(*p_source, THRESHOLD_STRESS));
    KRATOS_EXPECT_DOUBLE_EQ(read_scalar(*p_clone, PREVIOUS_CYCLE), previous_cycle_time);
    KRATOS_EXPECT_DOUBLE_EQ(read_scalar(*p_clone, CYCLE_PERIOD), cycle_period);
    KRATOS_EXPECT_DOUBLE_EQ(read_scalar(*p_clone, PREVIOUS_CYCLE_DAMAGE), damage);

    for (const double axial_strain : {valley, peak}) {
        const Vector source_response = advance(*p_source, axial_strain);
        const Vector clone_response = advance(*p_clone, axial_strain);
        KRATOS_EXPECT_VECTOR_NEAR(clone_response, source_response, 1.0e-6);
    }
    const double continued_cycles_to_failure = read_scalar(*p_source, CYCLES_TO_FAILURE);
    KRATOS_EXPECT_DOUBLE_EQ(read_scalar(*p_clone, CYCLES_TO_FAILURE), continued_cycles_to_failure);
    const int continued_cycles = read_cycles(*p_source);
    KRATOS_EXPECT_EQ(read_cycles(*p_clone), continued_cycles);

    advance(*p_clone, valley);
    advance(*p_clone, peak);
    KRATOS_EXPECT_EQ(read_cycles(*p_clone), continued_cycles + 1);
    KRATOS_EXPECT_EQ(read_cycles(*p_source), continued_cycles);
    KRATOS_EXPECT_DOUBLE_EQ(read_scalar(*p_source, CYCLES_TO_FAILURE), continued_cycles_to_failure);
}

} // namespace Kratos::Testing
