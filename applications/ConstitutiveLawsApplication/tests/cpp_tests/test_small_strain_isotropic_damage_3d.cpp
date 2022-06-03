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
//  Main authors:    Alejandro Cornejo
//                   Marcelo Raschi
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
#include "custom_constitutive/small_strain_isotropic_damage_3d.h"
#include "includes/model_part.h"
#include "geometries/tetrahedra_3d_4.h"

namespace Kratos
{
namespace Testing
{
typedef Node<3> NodeType;

KRATOS_TEST_CASE_IN_SUITE(_ConstitutiveLaw_SmallStrainIsotropicDamage3D, KratosConstitutiveLawsFastSuite)
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
    NodeType::Pointer p_node_1 = test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    NodeType::Pointer p_node_2 = test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    NodeType::Pointer p_node_3 = test_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    NodeType::Pointer p_node_4 = test_model_part.CreateNewNode(4, 0.0, 0.0, 1.0);
    Tetrahedra3D4<NodeType> geometry = Tetrahedra3D4<NodeType>(p_node_1, p_node_2, p_node_3, p_node_4);
    // Set material properties
    material_properties.SetValue(YOUNG_MODULUS, 6);
    material_properties.SetValue(POISSON_RATIO, 0.3);
    material_properties.SetValue(HARDENING_CURVE, 0);
    Vector stress_limits(2);
    stress_limits(0) = 1.5;
    stress_limits(1) = 3.0;
    material_properties.SetValue(STRESS_LIMITS, stress_limits);
    //Vector hardening_params(3);
    Vector hardening_params(1);
    hardening_params(0) = 0.5;
    material_properties.SetValue(HARDENING_PARAMETERS, hardening_params);
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
    SmallStrainIsotropicDamage3D cl = SmallStrainIsotropicDamage3D();

    //
    //  TESTS
    //

    //
    // Test: check correct behavior of internal and calculated variables
    //
    KRATOS_CHECK_IS_FALSE(cl.Check(material_properties, geometry, test_model_part.GetProcessInfo()));

    KRATOS_CHECK_IS_FALSE(cl.Has(STRAIN_ENERGY));  // = False, in order to use CalculateValue())
    KRATOS_CHECK_IS_FALSE(cl.Has(DAMAGE_VARIABLE));  // = False, in order to use CalculateValue())
    KRATOS_CHECK(cl.Has(INTERNAL_VARIABLES));  // = True
    Vector internal_variables_w(1);
    internal_variables_w[0] = 0.123;
    cl.SetValue(INTERNAL_VARIABLES, internal_variables_w, test_model_part.GetProcessInfo());
    Vector internal_variables_r;  // CL should internally resize it to 1
    cl.GetValue(INTERNAL_VARIABLES, internal_variables_r);
    KRATOS_CHECK_NEAR(internal_variables_r.size(), 1., 1.e-5);  // = True
    KRATOS_CHECK_NEAR(internal_variables_r[0], 0.123, 1.e-5);  // = True

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

    // Check CalculateValue(STRAIN)
    strain_vector(0) = 0.0000;
    strain_vector(1) = 0.0001;
    strain_vector(2) = 0.0002;
    strain_vector(3) = 0.0003;
    strain_vector(4) = 0.0004;
    strain_vector(5) = 0.0005;
    Vector strain_ref = strain_vector - initial_strain_w;
    Vector strain_r;
    cl.CalculateValue(cl_parameters, STRAIN, strain_r);
    KRATOS_CHECK_VECTOR_NEAR(strain_r, strain_ref, 1e-8);

    //Check CalculateVector(STRESSES)
    Vector stress_r;
    cl.CalculateValue(cl_parameters, STRESSES, stress_r);
    KRATOS_CHECK_VECTOR_NEAR(stress_r, stress_vector, 1e-8);

    // reset values
    p_initial_state->SetInitialStrainVector(ZeroVector(6));


    //
    // Test: exponential hardening model, load in traction
    //
    const double tolerance = 1e-4;
    double ref_damage_variable = 0.562727;
    double ref_strain_energy = 1.36795;
    Vector imposed_strain = ZeroVector(6);
    imposed_strain(0) = -0.0759;
    imposed_strain(1) =  0.7483;
    imposed_strain(2) =  0.1879;
    imposed_strain(3) =  0.5391;
    imposed_strain(4) =  0.0063;
    imposed_strain(5) = -0.3292;
    Vector ref_stress = ZeroVector(6);
    ref_stress(0) =  1.1490;
    ref_stress(1) =  2.8123;
    ref_stress(2) =  1.6813;
    ref_stress(3) =  0.5440;
    ref_stress(4) =  0.0063;
    ref_stress(5) = -0.3321;
    Matrix ref_C = ZeroMatrix(6, 6);
    ref_C(0,0) =  3.1672e+00; ref_C(0,1) =  6.2137e-01; ref_C(0,2) =  9.8019e-01; ref_C(0,3) = -1.7259e-01; ref_C(0,4) = -2.0169e-03; ref_C(0,5) =  1.0539e-01;
    ref_C(1,0) =  6.2137e-01; ref_C(1,1) =  1.3478e+00; ref_C(1,2) =  2.0793e-01; ref_C(1,3) = -4.2244e-01; ref_C(1,4) = -4.9367e-03; ref_C(1,5) =  2.5796e-01;
    ref_C(2,0) =  9.8019e-01; ref_C(2,1) =  2.0793e-01; ref_C(2,2) =  2.7511e+00; ref_C(2,3) = -2.5256e-01; ref_C(2,4) = -2.9514e-03; ref_C(2,5) =  1.5422e-01;
    ref_C(3,0) = -1.7259e-01; ref_C(3,1) = -4.2244e-01; ref_C(3,2) = -2.5256e-01; ref_C(3,3) =  9.2737e-01; ref_C(3,4) = -9.5492e-04; ref_C(3,5) =  4.9898e-02;
    ref_C(4,0) = -2.0169e-03; ref_C(4,1) = -4.9367e-03; ref_C(4,2) = -2.9514e-03; ref_C(4,3) = -9.5492e-04; ref_C(4,4) =  1.0090e+00; ref_C(4,5) =  5.8312e-04;
    ref_C(5,0) =  1.0539e-01; ref_C(5,1) =  2.5796e-01; ref_C(5,2) =  1.5422e-01; ref_C(5,3) =  4.9898e-02; ref_C(5,4) =  5.8312e-04; ref_C(5,5) =  9.7862e-01;

    // Simulate the call sequence of the element
    Vector dummy_vector;
    cl.InitializeMaterial(material_properties, geometry, dummy_vector);

    for (std::size_t t = 0; t < 10; ++t){
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

    // Check damage variable and strain energy
    double value;
    cl.CalculateValue(cl_parameters, DAMAGE_VARIABLE, value);
    KRATOS_CHECK_NEAR(ref_damage_variable, value, tolerance);
    cl.CalculateValue(cl_parameters, STRAIN_ENERGY, value);
    KRATOS_CHECK_NEAR(ref_strain_energy, value, tolerance);

    // Check stress
    for (std::size_t comp = 0; comp < 6; ++comp)
        KRATOS_CHECK_NEAR(stress_vector(comp), ref_stress(comp), tolerance);

    // Check constitutive tensor
    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t j = 0; j < 6; ++j)
            KRATOS_CHECK_NEAR(const_matrix(i, j), ref_C(j, i), tolerance);

    //
    // Test: trilinear hardening model, load in traction
    //
    ref_damage_variable = 0.510366;
    ref_strain_energy = 1.53176;
    imposed_strain(0) = -0.0759;
    imposed_strain(1) =  0.7483;
    imposed_strain(2) =  0.1879;
    imposed_strain(3) =  0.5391;
    imposed_strain(4) =  0.0063;
    imposed_strain(5) = -0.3292;
    ref_stress(0) =  1.28659;
    ref_stress(1) =  3.14916;
    ref_stress(2) =  1.88274;
    ref_stress(3) =  0.60914;
    ref_stress(4) =  0.00712;
    ref_stress(5) = -0.37197;
    ref_C(0,0) =  3.41440e+00; ref_C(0,1) =  3.72331e-01; ref_C(0,2) =  9.04190e-01; ref_C(0,3) = -2.55822e-01; ref_C(0,4) = -2.98958e-03; ref_C(0,5) =  1.56217e-01;
    ref_C(1,0) =  3.72331e-01; ref_C(1,1) =  7.17543e-01; ref_C(1,2) = -2.40482e-01; ref_C(1,3) = -6.26171e-01; ref_C(1,4) = -7.31753e-03; ref_C(1,5) =  3.82370e-01;
    ref_C(2,0) =  9.04191e-01; ref_C(2,1) = -2.40483e-01; ref_C(2,2) =  2.79766e+00; ref_C(2,3) = -3.74359e-01; ref_C(2,4) = -4.37481e-03; ref_C(2,5) =  2.28601e-01;
    ref_C(3,0) = -2.55822e-01; ref_C(3,1) = -6.26171e-01; ref_C(3,2) = -3.74359e-01; ref_C(3,3) =  1.00880e+00; ref_C(3,4) = -1.41543e-03; ref_C(3,5) =  7.39619e-02;
    ref_C(4,0) = -2.98958e-03; ref_C(4,1) = -7.31753e-03; ref_C(4,2) = -4.37481e-03; ref_C(4,3) = -1.41543e-03; ref_C(4,4) =  1.12991e+00; ref_C(4,5) =  8.64329e-04;
    ref_C(5,0) =  1.56217e-01; ref_C(5,1) =  3.82370e-01; ref_C(5,2) =  2.28601e-01; ref_C(5,3) =  7.39619e-02; ref_C(5,4) =  8.64329e-04; ref_C(5,5) =  1.08476e+00;

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

    for (std::size_t t = 0; t < 10; ++t){
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

    // Check damage variable and strain energy
    cl.CalculateValue(cl_parameters, DAMAGE_VARIABLE, value);
    KRATOS_CHECK_NEAR(ref_damage_variable, value, tolerance);
    cl.CalculateValue(cl_parameters, STRAIN_ENERGY, value);
    KRATOS_CHECK_NEAR(ref_strain_energy, value, tolerance);

    // Check stress
    for (std::size_t comp = 0; comp < 6; ++comp)
        KRATOS_CHECK_NEAR(stress_vector(comp), ref_stress(comp), tolerance);

    // Check constitutive tensor
    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t j = 0; j < 6; ++j)
            KRATOS_CHECK_NEAR(const_matrix(i, j), ref_C(j, i), tolerance);


}
} // namespace Testing
} // namespace Kratos
