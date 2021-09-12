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
//  Main authors:    Marcelo Raschi
//                   

// System includes

// External includes

// Project includes
#include "includes/process_info.h"
#include "testing/testing.h"
#include "containers/model.h"

// Application includes
#include "structural_mechanics_application_variables.h"
#include "constitutive_laws_application_variables.h"

// Constitutive law
#include "custom_constitutive/small_strain_isotropic_damage_traction_only_implex_3d.h"
#include "includes/model_part.h"
#include "geometries/tetrahedra_3d_4.h"

namespace Kratos
{
namespace Testing
{
typedef Node<3> NodeType;

KRATOS_TEST_CASE_IN_SUITE(_ConstitutiveLaw_SmallStrainIsotropicDamageTractionOnlyImplex3D_0, KratosConstitutiveLawsFastSuite)
{

    //
    //  CREATE LAW
    //

    ConstitutiveLaw::Parameters cl_parameters;
    Properties material_properties;
    Vector stress_vector(6), strain_vector(6);
    Matrix const_matrix;
    // Create gauss point
    Model current_model;
    ModelPart& test_model_part = current_model.CreateModelPart("Main");
    test_model_part.SetBufferSize(2);
    NodeType::Pointer p_node_1 = test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    NodeType::Pointer p_node_2 = test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    NodeType::Pointer p_node_3 = test_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    NodeType::Pointer p_node_4 = test_model_part.CreateNewNode(4, 0.0, 0.0, 1.0);
    Tetrahedra3D4<NodeType> geometry = Tetrahedra3D4<NodeType>(p_node_1, p_node_2, p_node_3, p_node_4);
    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions=cl_parameters.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    // Set required constitutive law parameters:
    cl_parameters.SetElementGeometry(geometry);
    cl_parameters.SetProcessInfo(test_model_part.GetProcessInfo());
    cl_parameters.SetMaterialProperties(material_properties);
    cl_parameters.SetStrainVector(strain_vector);
    cl_parameters.SetStressVector(stress_vector);
    cl_parameters.SetConstitutiveMatrix(const_matrix);
    // Create the CL
    SmallStrainIsotropicDamageTractionOnlyImplex3D cl = SmallStrainIsotropicDamageTractionOnlyImplex3D();

    //
    //  TESTS
    //

    //
    // Test 0: check correct behavior of internal and calculated variables
    //

    KRATOS_CHECK_IS_FALSE(cl.Has(SCALE_FACTOR));
    KRATOS_CHECK_IS_FALSE(cl.Has(STRAIN_ENERGY));
    KRATOS_CHECK_IS_FALSE(cl.Has(DAMAGE_VARIABLE));
    KRATOS_CHECK_IS_FALSE(cl.Has(STRAIN));
    KRATOS_CHECK(cl.Has(INTERNAL_VARIABLES));
    Vector internal_variables_w(2);
    internal_variables_w[0] = 0.123;
    internal_variables_w[1] = 0.456;
    cl.SetValue(INTERNAL_VARIABLES, internal_variables_w, test_model_part.GetProcessInfo());
    Vector internal_variables_r;
    cl.GetValue(INTERNAL_VARIABLES, internal_variables_r);  // CL resizes it
    KRATOS_CHECK_NEAR(internal_variables_r.size(), 2., 1.e-5);
    KRATOS_CHECK_NEAR(internal_variables_r[0], 0.123, 1.e-5);
    KRATOS_CHECK_NEAR(internal_variables_r[1], 0.456, 1.e-5);

    // Check CalculateValue(INITIAL_STRAIN_VECTOR)
    Vector initial_strain_w(6);
    initial_strain_w(0) = 0.000000;
    initial_strain_w(1) = 0.000001;
    initial_strain_w(2) = 0.000002;
    initial_strain_w(3) = 0.000003;
    initial_strain_w(4) = 0.000004;
    initial_strain_w(5) = 0.000005;
    InitialState::Pointer p_initial_state = Kratos::make_intrusive<InitialState>
	    (initial_strain_w, InitialState::InitialImposingType::STRAIN_ONLY);
    cl.SetInitialState(p_initial_state);
    Vector initial_strain_r;  // Should be accordingly resized in the CL
    cl.CalculateValue(cl_parameters, INITIAL_STRAIN_VECTOR, initial_strain_r);
    KRATOS_CHECK_VECTOR_NEAR(initial_strain_r, initial_strain_w, 1e-8);
    p_initial_state->SetInitialStrainVector(ZeroVector(6));  // reset values

    // Check CalculateValue(STRAIN)
    strain_vector(0) = 0.0000;
    strain_vector(1) = 0.0001;
    strain_vector(2) = 0.0002;
    strain_vector(3) = 0.0003;
    strain_vector(4) = 0.0004;
    strain_vector(5) = 0.0005;
    Vector value_vector = ZeroVector(0);
    cl.CalculateValue(cl_parameters, STRAIN, value_vector);  // from Elastic CL
    KRATOS_CHECK_VECTOR_NEAR(value_vector, strain_vector, 1e-8);

    // Check CalculateValue(STRESSES)
    stress_vector(0) = 0.0000;
    stress_vector(1) = 0.0001;
    stress_vector(2) = 0.0002;
    stress_vector(3) = 0.0003;
    stress_vector(4) = 0.0004;
    stress_vector(5) = 0.0005;
    value_vector = ZeroVector(0);
    cl.CalculateValue(cl_parameters, STRESSES, value_vector);  // from Elastic CL
    KRATOS_CHECK_VECTOR_NEAR(value_vector, stress_vector, 1e-8);
    cl.CalculateValue(cl_parameters, CAUCHY_STRESS_VECTOR, value_vector);  // from Elastic CL
    KRATOS_CHECK_VECTOR_NEAR(value_vector, stress_vector, 1e-8);

    //
    // Test 1: exponential hardening model, load in traction
    //

    // Set material properties for this test
    material_properties.SetValue(YOUNG_MODULUS, 6);
    material_properties.SetValue(POISSON_RATIO, 0.3);
    material_properties.SetValue(HARDENING_CURVE, 0);
    Vector stress_limits(2);
    stress_limits(0) = 1.5;
    stress_limits(1) = 3.0;
    material_properties.SetValue(STRESS_LIMITS, stress_limits);
    Vector hardening_params(1);
    hardening_params(0) = 0.5;
    material_properties.SetValue(HARDENING_PARAMETERS, hardening_params);

    // Simulate the call sequence of the element
    cl.Check(material_properties, geometry, test_model_part.GetProcessInfo());
    Vector dummy_vector;
    cl.InitializeMaterial(material_properties, geometry, dummy_vector);
    Vector imposed_strain = ZeroVector(6);
    imposed_strain(0) = -0.0759;
    imposed_strain(1) =  0.7483;
    imposed_strain(2) =  0.1879;
    imposed_strain(3) =  0.5391;
    imposed_strain(4) =  0.0063;
    imposed_strain(5) = -0.3292;

    double time = 0;
    for (int t = 0; t < 10; ++t){

        time += 1. / 9.;
        test_model_part.CloneTimeStep(time);

        for (std::size_t comp = 0; comp < 6; ++comp)
            strain_vector(comp) = t / 9 * imposed_strain(comp);

        if (cl.RequiresInitializeMaterialResponse()){
            cl.InitializeMaterialResponseCauchy(cl_parameters);
        }

        cl.CalculateMaterialResponseCauchy(cl_parameters);

        if (cl.RequiresFinalizeMaterialResponse()){
            cl.FinalizeMaterialResponseCauchy(cl_parameters);
        }
    }

    // Checks
    const double tolerance = 1e-4;
    double value_double;

    cl.CalculateValue(cl_parameters, SCALE_FACTOR, value_double);
    KRATOS_CHECK_NEAR(0.24946, value_double, tolerance);

    cl.CalculateValue(cl_parameters, DAMAGE_VARIABLE, value_double);  // from Damage CL
    KRATOS_CHECK_NEAR(0.56273, value_double, tolerance);

    cl.CalculateValue(cl_parameters, STRAIN_ENERGY, value_double);  // from Damage CL
    KRATOS_CHECK_NEAR(1.36795, value_double, tolerance);

    cl.GetValue(INTERNAL_VARIABLES, internal_variables_r);
    KRATOS_CHECK_NEAR(internal_variables_r.size(), 2., tolerance);
    KRATOS_CHECK_NEAR(internal_variables_r[0], 2.50135, tolerance);  // r_{t}
    KRATOS_CHECK_NEAR(internal_variables_r[1], 0.61237, tolerance);  // r_{t-1}

    cl.CalculateValue(cl_parameters, INITIAL_STRAIN_VECTOR, value_vector);  // from Damage CL
    KRATOS_CHECK_VECTOR_NEAR(value_vector, ZeroVector(6), tolerance);

    Vector ref_strain = imposed_strain;
    KRATOS_CHECK_VECTOR_NEAR(strain_vector, ref_strain, tolerance); // check strain from Cl properties
    cl.CalculateValue(cl_parameters, STRAIN, value_vector);    // from Damage CL
    KRATOS_CHECK_VECTOR_NEAR(value_vector, ref_strain, tolerance);  // check strain from CL CalculateValue

    Vector ref_stress = ZeroVector(6);
    ref_stress(0) =  2.62765;
    ref_stress(1) =  6.43165;
    ref_stress(2) =  3.84519;
    ref_stress(3) =  1.24408;
    ref_stress(4) =  0.01454;
    ref_stress(5) = -0.75969;
    KRATOS_CHECK_VECTOR_NEAR(stress_vector, ref_stress, tolerance);
    //cl.CalculateValue(cl_parameters, CAUCHY_STRESS_VECTOR, value_vector);  // from Elastic CL
    //KRATOS_CHECK_VECTOR_NEAR(value_vector, ref_stress, tolerance);  // check strain from CL CalculateValue

    Matrix ref_C = ZeroMatrix(6, 6);
    ref_C(0,0) = 8.07692; ref_C(0,1) = 3.46154; ref_C(0,2) = 3.46154; ref_C(0,3) = 0      ; ref_C(0,4) = 0      ; ref_C(0,5) = 0      ;
    ref_C(1,0) = 3.46154; ref_C(1,1) = 8.07692; ref_C(1,2) = 3.46154; ref_C(1,3) = 0      ; ref_C(1,4) = 0      ; ref_C(1,5) = 0      ;
    ref_C(2,0) = 3.46154; ref_C(2,1) = 3.46154; ref_C(2,2) = 8.07692; ref_C(2,3) = 0      ; ref_C(2,4) = 0      ; ref_C(2,5) = 0      ;
    ref_C(3,0) = 0      ; ref_C(3,1) = 0      ; ref_C(3,2) = 0      ; ref_C(3,3) = 2.30769; ref_C(3,4) = 0      ; ref_C(3,5) = 0      ;
    ref_C(4,0) = 0      ; ref_C(4,1) = 0      ; ref_C(4,2) = 0      ; ref_C(4,3) = 0      ; ref_C(4,4) = 2.30769; ref_C(4,5) = 0      ;
    ref_C(5,0) = 0      ; ref_C(5,1) = 0      ; ref_C(5,2) = 0      ; ref_C(5,3) = 0      ; ref_C(5,4) = 0      ; ref_C(5,5) = 2.30769;
    KRATOS_CHECK_MATRIX_NEAR(const_matrix, ref_C, tolerance);

    //
    // Test 2: trilinear hardening model, load in traction
    //

    // Set material properties for this test
    material_properties.SetValue(YOUNG_MODULUS, 6);
    material_properties.SetValue(POISSON_RATIO, 0.3);
    material_properties.SetValue(HARDENING_CURVE, 1);
    stress_limits.resize(3);
    stress_limits(0) = 1.5;
    stress_limits(1) = 2.0;
    stress_limits(2) = 3.0;
    material_properties.SetValue(STRESS_LIMITS, stress_limits);
    hardening_params.resize(3);
    hardening_params(0) = 0.6;
    hardening_params(1) = 0.4;
    hardening_params(2) = 0.0;
    material_properties.SetValue(HARDENING_PARAMETERS, hardening_params);

    // Simulate the call sequence of the element
    cl.Check(material_properties, geometry, test_model_part.GetProcessInfo());
    cl.InitializeMaterial(material_properties, geometry, dummy_vector);

    time = 0;
    for (int t = 0; t < 10; ++t){

        time += 1. / 9.;
        test_model_part.CloneTimeStep(time);

        for (std::size_t comp = 0; comp < 6; ++comp)
            strain_vector(comp) = t / 9 * imposed_strain(comp);

        if (cl.RequiresInitializeMaterialResponse()){
            cl.InitializeMaterialResponseCauchy(cl_parameters);
        }

        cl.CalculateMaterialResponseCauchy(cl_parameters);

        if (cl.RequiresFinalizeMaterialResponse()){
            cl.FinalizeMaterialResponseCauchy(cl_parameters);
        }
    }

    // Checks
    cl.CalculateValue(cl_parameters, DAMAGE_VARIABLE, value_double);
    KRATOS_CHECK_NEAR(0.51037, value_double, tolerance);

    cl.CalculateValue(cl_parameters, STRAIN_ENERGY, value_double);
    KRATOS_CHECK_NEAR(1.53176, value_double, tolerance);

    cl.CalculateValue(cl_parameters, SCALE_FACTOR, value_double);
    KRATOS_CHECK_NEAR(0.36976, value_double, tolerance);

    cl.GetValue(INTERNAL_VARIABLES, internal_variables_r);
    KRATOS_CHECK_NEAR(internal_variables_r.size(), 2., tolerance);
    KRATOS_CHECK_NEAR(internal_variables_r[0], 2.50135, tolerance);  // r_{t}
    KRATOS_CHECK_NEAR(internal_variables_r[1], 0.61237, tolerance);  // r_{t-1}

    //cl.CalculateValue(cl_parameters, INITIAL_STRAIN_VECTOR, value_vector);    // from Damage CL
    //KRATOS_CHECK_VECTOR_NEAR(value_vector, ZeroVector(6), tolerance);

    cl.CalculateValue(cl_parameters, STRAIN, value_vector);    // from Damage CL
    ref_strain = imposed_strain;
    KRATOS_CHECK_VECTOR_NEAR(value_vector, ref_strain, tolerance);  // check strain from CL CalculateValue
    KRATOS_CHECK_VECTOR_NEAR(strain_vector, ref_strain, tolerance); // check strain from Cl properties

    ref_stress(0) =  2.62765;
    ref_stress(1) =  6.43165;
    ref_stress(2) =  3.84519;
    ref_stress(3) =  1.24408;
    ref_stress(4) =  0.01454;
    ref_stress(5) = -0.75969;
    KRATOS_CHECK_VECTOR_NEAR(stress_vector, ref_stress, tolerance);  // check stress from CL properties
    //cl.CalculateValue(cl_parameters, CAUCHY_STRESS_VECTOR, value_vector);
    //KRATOS_CHECK_VECTOR_NEAR(value_vector, ref_stress, tolerance);   // check stress from elements

    ref_C(0,0) = 8.07692; ref_C(0,1) = 3.46154; ref_C(0,2) = 3.46154; ref_C(0,3) = 0      ; ref_C(0,4) = 0      ; ref_C(0,5) = 0      ;
    ref_C(1,0) = 3.46154; ref_C(1,1) = 8.07692; ref_C(1,2) = 3.46154; ref_C(1,3) = 0      ; ref_C(1,4) = 0      ; ref_C(1,5) = 0      ;
    ref_C(2,0) = 3.46154; ref_C(2,1) = 3.46154; ref_C(2,2) = 8.07692; ref_C(2,3) = 0      ; ref_C(2,4) = 0      ; ref_C(2,5) = 0      ;
    ref_C(3,0) = 0      ; ref_C(3,1) = 0      ; ref_C(3,2) = 0      ; ref_C(3,3) = 2.30769; ref_C(3,4) = 0      ; ref_C(3,5) = 0      ;
    ref_C(4,0) = 0      ; ref_C(4,1) = 0      ; ref_C(4,2) = 0      ; ref_C(4,3) = 0      ; ref_C(4,4) = 2.30769; ref_C(4,5) = 0      ;
    ref_C(5,0) = 0      ; ref_C(5,1) = 0      ; ref_C(5,2) = 0      ; ref_C(5,3) = 0      ; ref_C(5,4) = 0      ; ref_C(5,5) = 2.30769;
    KRATOS_CHECK_MATRIX_NEAR(const_matrix, ref_C, tolerance);
}


KRATOS_TEST_CASE_IN_SUITE(_ConstitutiveLaw_SmallStrainIsotropicDamageTractionOnlyImplex3D_1, KratosConstitutiveLawsFastSuite)
{

    //
    //  CREATE LAW
    //

    ConstitutiveLaw::Parameters cl_parameters;
    Properties material_properties;
    Vector stress_vector(6), strain_vector(6);
    Matrix const_matrix;
    // Create gauss point
    Model current_model;
    ModelPart& test_model_part = current_model.CreateModelPart("Main");
    test_model_part.SetBufferSize(2);
    NodeType::Pointer p_node_1 = test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    NodeType::Pointer p_node_2 = test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    NodeType::Pointer p_node_3 = test_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    NodeType::Pointer p_node_4 = test_model_part.CreateNewNode(4, 0.0, 0.0, 1.0);
    Tetrahedra3D4<NodeType> geometry = Tetrahedra3D4<NodeType>(p_node_1, p_node_2, p_node_3, p_node_4);
    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions=cl_parameters.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    // Set required constitutive law parameters:
    cl_parameters.SetElementGeometry(geometry);
    cl_parameters.SetProcessInfo(test_model_part.GetProcessInfo());
    cl_parameters.SetMaterialProperties(material_properties);
    cl_parameters.SetStrainVector(strain_vector);
    cl_parameters.SetStressVector(stress_vector);
    cl_parameters.SetConstitutiveMatrix(const_matrix);
    // Create the CL
    SmallStrainIsotropicDamageTractionOnlyImplex3D cl = SmallStrainIsotropicDamageTractionOnlyImplex3D();

    //
    //  TESTS
    //

    //
    // Test 1: exponential hardening model, load in traction
    //

    // Set material properties for this test
    material_properties.SetValue(YOUNG_MODULUS, 6);
    material_properties.SetValue(POISSON_RATIO, 0.3);
    material_properties.SetValue(HARDENING_CURVE, 0);
    Vector stress_limits(2);
    stress_limits(0) = 1.5;
    stress_limits(1) = 3.0;
    material_properties.SetValue(STRESS_LIMITS, stress_limits);
    Vector hardening_params(1);
    hardening_params(0) = 0.5;
    material_properties.SetValue(HARDENING_PARAMETERS, hardening_params);

    // Simulate the call sequence of the element
    cl.Check(material_properties, geometry, test_model_part.GetProcessInfo());
    Vector dummy_vector;
    cl.InitializeMaterial(material_properties, geometry, dummy_vector);
    Vector imposed_strain = ZeroVector(6);
    imposed_strain(0) = -0.0759;
    imposed_strain(1) =  0.7483;
    imposed_strain(2) =  0.1879;
    imposed_strain(3) =  0.5391;
    imposed_strain(4) =  0.0063;
    imposed_strain(5) = -0.3292;

    double time = 0;
    for (int t = 0; t < 10; ++t){

        time += 1. / 9.;
        test_model_part.CloneTimeStep(time);

        for (std::size_t comp = 0; comp < 6; ++comp)
            strain_vector(comp) = t / 9 * imposed_strain(comp);

        if (cl.RequiresInitializeMaterialResponse()){
            cl.InitializeMaterialResponseCauchy(cl_parameters);
        }

        cl.CalculateMaterialResponseCauchy(cl_parameters);

        if (cl.RequiresFinalizeMaterialResponse()){
            cl.FinalizeMaterialResponseCauchy(cl_parameters);
        }
    }

    for (int t = 8; t > -14; --t){

        time += 1. / 9.;
        test_model_part.CloneTimeStep(time);

        for (std::size_t comp = 0; comp < 6; ++comp)
            strain_vector(comp) = t / 9 * imposed_strain(comp);

        if (cl.RequiresInitializeMaterialResponse()){
            cl.InitializeMaterialResponseCauchy(cl_parameters);
        }

        cl.CalculateMaterialResponseCauchy(cl_parameters);

        if (cl.RequiresFinalizeMaterialResponse()){
            cl.FinalizeMaterialResponseCauchy(cl_parameters);
        }
    }

    // Checks
    const double tolerance = 1e-4;
    double value_double;
    Vector value_vector;

    cl.CalculateValue(cl_parameters, DAMAGE_VARIABLE, value_double);
    KRATOS_CHECK_NEAR(0.562727, value_double, tolerance);

    cl.CalculateValue(cl_parameters, STRAIN_ENERGY, value_double);
    KRATOS_CHECK_NEAR(1.36795, value_double, tolerance);

    Vector internal_variables_r = ZeroVector(0);
    cl.GetValue(INTERNAL_VARIABLES, internal_variables_r);
    KRATOS_CHECK_NEAR(internal_variables_r.size(), 2., tolerance);
    KRATOS_CHECK_NEAR(internal_variables_r[0], 2.50135, tolerance);  // r_{t}
    KRATOS_CHECK_NEAR(internal_variables_r[1], 2.50135, tolerance);  // r_{t-1}

    cl.CalculateValue(cl_parameters, STRAIN, value_vector);    // from Damage CL
    Vector ref_strain = ZeroVector(6);
    ref_strain(0) =  0.0759;
    ref_strain(1) = -0.7483;
    ref_strain(2) = -0.1879;
    ref_strain(3) = -0.5391;
    ref_strain(4) = -0.0063;
    ref_strain(5) =  0.3292;
    KRATOS_CHECK_VECTOR_NEAR(value_vector, ref_strain, tolerance);
    KRATOS_CHECK_VECTOR_NEAR(strain_vector, ref_strain, tolerance);

    Vector ref_stress = ZeroVector(6);
    ref_stress(0) = -1.14900;
    ref_stress(1) = -2.81239;
    ref_stress(2) = -1.68140;
    ref_stress(3) = -0.54400;
    ref_stress(4) = -0.00636;
    ref_stress(5) =  0.33219;
    KRATOS_CHECK_VECTOR_NEAR(stress_vector, ref_stress, tolerance);


    //cl.CalculateValue(cl_parameters, CAUCHY_STRESS_VECTOR, value_vector);    // from Element
    //ref_stress(0) = -2.62765;
    //ref_stress(1) = -6.43165;
    //ref_stress(2) = -3.84519;
    //ref_stress(3) = -1.24408;
    //ref_stress(4) = -0.01454;
    //ref_stress(5) =  0.75969;
    //KRATOS_CHECK_VECTOR_NEAR(value_vector, ref_stress, tolerance);
    //KRATOS_CHECK_VECTOR_NEAR(stress_vector, ref_stress, tolerance);
    //
    //Matrix ref_C = ZeroMatrix(6, 6);
    //ref_C(0,0) = 8.07692; ref_C(0,1) = 3.46154; ref_C(0,2) = 3.46154; ref_C(0,3) = 0      ; ref_C(0,4) = 0      ; ref_C(0,5) = 0      ;
    //ref_C(1,0) = 3.46154; ref_C(1,1) = 8.07692; ref_C(1,2) = 3.46154; ref_C(1,3) = 0      ; ref_C(1,4) = 0      ; ref_C(1,5) = 0      ;
    //ref_C(2,0) = 3.46154; ref_C(2,1) = 3.46154; ref_C(2,2) = 8.07692; ref_C(2,3) = 0      ; ref_C(2,4) = 0      ; ref_C(2,5) = 0      ;
    //ref_C(3,0) = 0      ; ref_C(3,1) = 0      ; ref_C(3,2) = 0      ; ref_C(3,3) = 2.30769; ref_C(3,4) = 0      ; ref_C(3,5) = 0      ;
    //ref_C(4,0) = 0      ; ref_C(4,1) = 0      ; ref_C(4,2) = 0      ; ref_C(4,3) = 0      ; ref_C(4,4) = 2.30769; ref_C(4,5) = 0      ;
    //ref_C(5,0) = 0      ; ref_C(5,1) = 0      ; ref_C(5,2) = 0      ; ref_C(5,3) = 0      ; ref_C(5,4) = 0      ; ref_C(5,5) = 2.30769;
    //KRATOS_CHECK_MATRIX_NEAR(const_matrix, ref_C, tolerance);

    Matrix ref_C = ZeroMatrix(6, 6);
    ref_C(0,0) = 3.53182; ref_C(0,1) = 1.51364; ref_C(0,2) = 1.51364; ref_C(0,3) = 0.0    ; ref_C(0,4) = 0.0    ; ref_C(0,5) = 0.0    ;
    ref_C(1,0) = 1.51364; ref_C(1,1) = 3.53182; ref_C(1,2) = 1.51364; ref_C(1,3) = 0.0    ; ref_C(1,4) = 0.0    ; ref_C(1,5) = 0.0    ;
    ref_C(2,0) = 1.51364; ref_C(2,1) = 1.51364; ref_C(2,2) = 3.53182; ref_C(2,3) = 0.0    ; ref_C(2,4) = 0.0    ; ref_C(2,5) = 0.0    ;
    ref_C(3,0) = 0.0    ; ref_C(3,1) = 0.0    ; ref_C(3,2) = 0.0    ; ref_C(3,3) = 1.00909; ref_C(3,4) = 0.0    ; ref_C(3,5) = 0.0    ;
    ref_C(4,0) = 0.0    ; ref_C(4,1) = 0.0    ; ref_C(4,2) = 0.0    ; ref_C(4,3) = 0.0    ; ref_C(4,4) = 1.00909; ref_C(4,5) = 0.0    ;
    ref_C(5,0) = 0.0    ; ref_C(5,1) = 0.0    ; ref_C(5,2) = 0.0    ; ref_C(5,3) = 0.0    ; ref_C(5,4) = 0.0    ; ref_C(5,5) = 1.00909;
    KRATOS_CHECK_MATRIX_NEAR(const_matrix, ref_C, tolerance);

    //
    // Test 2: trilinear hardening model, load in traction
    //

    // Update properties for this test
    material_properties.SetValue(HARDENING_CURVE, 1);
    stress_limits.resize(3);
    stress_limits(0) = 1.5;
    stress_limits(1) = 2.0;
    stress_limits(2) = 3.0;
    material_properties.SetValue(STRESS_LIMITS, stress_limits);
    hardening_params.resize(3);
    hardening_params(0) = 0.6;
    hardening_params(1) = 0.4;
    hardening_params(2) = 0.0;
    material_properties.SetValue(HARDENING_PARAMETERS, hardening_params);

    // Simulate the call sequence of the element
    cl.Check(material_properties, geometry, test_model_part.GetProcessInfo());
    cl.InitializeMaterial(material_properties, geometry, dummy_vector);

    time = 0.;
    for (int t = 0; t < 10; ++t){

        time += 1. / 9.;
        test_model_part.CloneTimeStep(time);

        for (std::size_t comp = 0; comp < 6; ++comp)
            strain_vector(comp) = t / 9 * imposed_strain(comp);

        if (cl.RequiresInitializeMaterialResponse()){
            cl.InitializeMaterialResponseCauchy(cl_parameters);
        }

        cl.CalculateMaterialResponseCauchy(cl_parameters);

        if (cl.RequiresFinalizeMaterialResponse()){
            cl.FinalizeMaterialResponseCauchy(cl_parameters);
        }
    }

    for (int t = 8; t > -14; --t){

        time += 1. / 9.;
        test_model_part.CloneTimeStep(time);

        for (std::size_t comp = 0; comp < 6; ++comp)
            strain_vector(comp) = t / 9 * imposed_strain(comp);

        if (cl.RequiresInitializeMaterialResponse()){
            cl.InitializeMaterialResponseCauchy(cl_parameters);
        }

        cl.CalculateMaterialResponseCauchy(cl_parameters);

        if (cl.RequiresFinalizeMaterialResponse()){
            cl.FinalizeMaterialResponseCauchy(cl_parameters);
        }
    }

    // Checks
    cl.CalculateValue(cl_parameters, DAMAGE_VARIABLE, value_double);
    KRATOS_CHECK_NEAR(0.510366, value_double, tolerance);

    cl.CalculateValue(cl_parameters, STRAIN_ENERGY, value_double);
    KRATOS_CHECK_NEAR(1.53176, value_double, tolerance);

    internal_variables_r = ZeroVector(0);
    cl.GetValue(INTERNAL_VARIABLES, internal_variables_r);
    KRATOS_CHECK_NEAR(internal_variables_r.size(), 2., tolerance);
    KRATOS_CHECK_NEAR(internal_variables_r[0], 2.50135, tolerance);  // r_{t}
    KRATOS_CHECK_NEAR(internal_variables_r[1], 2.50135, tolerance);  // r_{t-1}

    cl.CalculateValue(cl_parameters, STRAIN, value_vector);    // from Damage CL
    ref_strain(0) =  0.0759;
    ref_strain(1) = -0.7483;
    ref_strain(2) = -0.1879;
    ref_strain(3) = -0.5391;
    ref_strain(4) = -0.0063;
    ref_strain(5) =  0.3292;
    KRATOS_CHECK_VECTOR_NEAR(value_vector, ref_strain, tolerance);  // check strain from CL CalculateValue
    KRATOS_CHECK_VECTOR_NEAR(strain_vector, ref_strain, tolerance); // check strain from Cl properties

    ref_stress(0) = -1.28659;
    ref_stress(1) = -3.14916;
    ref_stress(2) = -1.88274;
    ref_stress(3) = -0.60914;
    ref_stress(4) = -0.00712;
    ref_stress(5) =  0.37197;
    KRATOS_CHECK_VECTOR_NEAR(stress_vector, ref_stress, tolerance);

    //cl.CalculateValue(cl_parameters, CAUCHY_STRESS_VECTOR, value_vector);    // from Element
    //ref_stress(0) = -2.62765;
    //ref_stress(1) = -6.43165;
    //ref_stress(2) = -3.84519;
    //ref_stress(3) = -1.24408;
    //ref_stress(4) = -0.01454;
    //ref_stress(5) =  0.75969;
    //KRATOS_CHECK_VECTOR_NEAR(value_vector, ref_stress, tolerance);
    //KRATOS_CHECK_VECTOR_NEAR(stress_vector, ref_stress, tolerance);
    //
    //ref_C(0,0) = 8.07692; ref_C(0,1) = 3.46154; ref_C(0,2) = 3.46154; ref_C(0,3) = 0      ; ref_C(0,4) = 0      ; ref_C(0,5) = 0      ;
    //ref_C(1,0) = 3.46154; ref_C(1,1) = 8.07692; ref_C(1,2) = 3.46154; ref_C(1,3) = 0      ; ref_C(1,4) = 0      ; ref_C(1,5) = 0      ;
    //ref_C(2,0) = 3.46154; ref_C(2,1) = 3.46154; ref_C(2,2) = 8.07692; ref_C(2,3) = 0      ; ref_C(2,4) = 0      ; ref_C(2,5) = 0      ;
    //ref_C(3,0) = 0      ; ref_C(3,1) = 0      ; ref_C(3,2) = 0      ; ref_C(3,3) = 2.30769; ref_C(3,4) = 0      ; ref_C(3,5) = 0      ;
    //ref_C(4,0) = 0      ; ref_C(4,1) = 0      ; ref_C(4,2) = 0      ; ref_C(4,3) = 0      ; ref_C(4,4) = 2.30769; ref_C(4,5) = 0      ;
    //ref_C(5,0) = 0      ; ref_C(5,1) = 0      ; ref_C(5,2) = 0      ; ref_C(5,3) = 0      ; ref_C(5,4) = 0      ; ref_C(5,5) = 2.30769;
    //KRATOS_CHECK_MATRIX_NEAR(const_matrix, ref_C, tolerance);

    ref_C(0,0) = 3.95474; ref_C(0,1) = 1.69489; ref_C(0,2) = 1.69489; ref_C(0,3) = 0.0    ; ref_C(0,4) = 0.0    ; ref_C(0,5) = 0.0    ;
    ref_C(1,0) = 1.69489; ref_C(1,1) = 3.95474; ref_C(1,2) = 1.69489; ref_C(1,3) = 0.0    ; ref_C(1,4) = 0.0    ; ref_C(1,5) = 0.0    ;
    ref_C(2,0) = 1.69489; ref_C(2,1) = 1.69489; ref_C(2,2) = 3.95474; ref_C(2,3) = 0.0    ; ref_C(2,4) = 0.0    ; ref_C(2,5) = 0.0    ;
    ref_C(3,0) = 0.0    ; ref_C(3,1) = 0.0    ; ref_C(3,2) = 0.0    ; ref_C(3,3) = 1.12992; ref_C(3,4) = 0.0    ; ref_C(3,5) = 0.0    ;
    ref_C(4,0) = 0.0    ; ref_C(4,1) = 0.0    ; ref_C(4,2) = 0.0    ; ref_C(4,3) = 0.0    ; ref_C(4,4) = 1.12992; ref_C(4,5) = 0.0    ;
    ref_C(5,0) = 0.0    ; ref_C(5,1) = 0.0    ; ref_C(5,2) = 0.0    ; ref_C(5,3) = 0.0    ; ref_C(5,4) = 0.0    ; ref_C(5,5) = 1.12992;
    KRATOS_CHECK_MATRIX_NEAR(const_matrix, ref_C, tolerance);
}
} // namespace Testing
} // namespace Kratos
