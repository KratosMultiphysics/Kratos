// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                     license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//

// System includes

// External includes

// Project includes
#include "includes/process_info.h"
#include "testing/testing.h"
#include "containers/model.h"

// Application includes
#include "structural_mechanics_application_variables.h"

// Integrator
#include "custom_advanced_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_damage.h"

// Yield surfaces
#include "custom_advanced_constitutive/yield_surfaces/generic_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/von_mises_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/modified_mohr_coulomb_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/rankine_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/simo_ju_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/drucker_prager_yield_surface.h"
#include "custom_advanced_constitutive/yield_surfaces/tresca_yield_surface.h"

// Plastic potentials
#include "custom_advanced_constitutive/plastic_potentials/generic_plastic_potential.h"
#include "custom_advanced_constitutive/plastic_potentials/von_mises_plastic_potential.h"
#include "custom_advanced_constitutive/plastic_potentials/tresca_plastic_potential.h"
#include "custom_advanced_constitutive/plastic_potentials/modified_mohr_coulomb_plastic_potential.h"
#include "custom_advanced_constitutive/plastic_potentials/drucker_prager_plastic_potential.h"

// Constitutive law
#include "custom_advanced_constitutive/generic_small_strain_isotropic_damage.h"
#include "includes/model_part.h"
#include "geometries/tetrahedra_3d_4.h"

namespace Kratos
{
namespace Testing
{
// We test the associated damage Constitutive laws...
typedef Node<3> NodeType;

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawIntegrateStressDamageInternalVariables,
                          KratosStructuralMechanicsFastSuite)
{
    //
    // Test: check correct behavior of internal and calculated variables
    //

    typedef GenericSmallStrainIsotropicDamage<GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> MC;
    Model current_model;
    ModelPart &test_model_part = current_model.CreateModelPart("Main");
    MC cl = MC();

    KRATOS_CHECK(cl.Has(DAMAGE));  // = True
    KRATOS_CHECK(cl.Has(THRESHOLD));  // = True
    KRATOS_CHECK(cl.Has(UNIAXIAL_STRESS));  // = True
    KRATOS_CHECK(cl.Has(INTERNAL_VARIABLES));  // = True

    Vector internal_variables_w(3);
    internal_variables_w[0] = 0.0;
    internal_variables_w[1] = 0.1;
    internal_variables_w[2] = 0.2;
    cl.SetValue(INTERNAL_VARIABLES, internal_variables_w, test_model_part.GetProcessInfo());
    Vector internal_variables_r;
    cl.GetValue(INTERNAL_VARIABLES, internal_variables_r);

    KRATOS_CHECK_NEAR(internal_variables_r.size(), 3., 1.e-5);  // = True
    KRATOS_CHECK_NEAR(internal_variables_r[0], 0.0, 1.e-5);  // = True
    KRATOS_CHECK_NEAR(internal_variables_r[1], 0.1, 1.e-5);  // = True
    KRATOS_CHECK_NEAR(internal_variables_r[2], 0.2, 1.e-5);  // = True
}


KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawIntegrateStressDamageLinear, KratosStructuralMechanicsFastSuite)
{
    typedef GenericSmallStrainIsotropicDamage<GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> MC;
    typedef GenericSmallStrainIsotropicDamage<GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> VM;
    typedef GenericSmallStrainIsotropicDamage<GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> DP;
    typedef GenericSmallStrainIsotropicDamage<GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> T;
    typedef GenericSmallStrainIsotropicDamage<GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> SJ;
    typedef GenericSmallStrainIsotropicDamage<GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> R;

    ConstitutiveLaw::Parameters cl_parameters;
    Properties material_properties;
    Vector stress_vector, strain_vector;

    Model current_model;

    ModelPart& test_model_part = current_model.CreateModelPart("Main");

    NodeType::Pointer p_node_1 = test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    NodeType::Pointer p_node_2 = test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    NodeType::Pointer p_node_3 = test_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    NodeType::Pointer p_node_4 = test_model_part.CreateNewNode(4, 0.0, 0.0, 1.0);

    Tetrahedra3D4<NodeType> Geom = Tetrahedra3D4<NodeType>(p_node_1, p_node_2, p_node_3, p_node_4);

    stress_vector = ZeroVector(6);
    stress_vector[0] = 5.40984e+06;
    stress_vector[1] = 5.40984e+06;
    stress_vector[2] = 1.91803e+07;
    stress_vector[3] = 0.0;
    stress_vector[4] = 0.0;
    stress_vector[5] = 1.45804e-10;

    strain_vector = ZeroVector(6);
    strain_vector[0] = 0.0;
    strain_vector[1] = 0.0;
    strain_vector[2] = 8.0e-5;
    strain_vector[3] = 0.0;
    strain_vector[4] = 0.0;
    strain_vector[5] = 1.6941e-21;

    material_properties.SetValue(YOUNG_MODULUS, 210e9);
    material_properties.SetValue(POISSON_RATIO, 0.22);
    material_properties.SetValue(YIELD_STRESS_COMPRESSION, 3.0e6);
    material_properties.SetValue(YIELD_STRESS_TENSION, 3.0e6);
    material_properties.SetValue(FRICTION_ANGLE, 32.0);
    material_properties.SetValue(DILATANCY_ANGLE, 16.0);
    material_properties.SetValue(FRACTURE_ENERGY, 1000.0);
    material_properties.SetValue(SOFTENING_TYPE, 0);

    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions=cl_parameters.GetOptions();
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
    MC MohrCoulombCL = MC();
    VM VonMisesCL = VM();
    DP DruckerPragerCL = DP();
    T TrescaCL = T();
    R RankineCL = R();
    SJ SimoJuCL = SJ();

    std::vector<double> MCres, VMres, DPres, Tres, Rres, SJres;
    MCres = {1.10203e+06,1.10203e+06,3.9072e+06,0,0,2.97015e-11};
    VMres = {1.10203e+06,1.10203e+06,3.9072e+06,0,0,2.97015e-11};
    DPres = {788415,788415,2.79529e+06,0,0,2.12491e-11};
    Tres  = {1.10203e+06,1.10203e+06,3.9072e+06,0,0,2.97015e-11};
    Rres  = {763601,763601,2.70731e+06,0,0,2.05803e-11};
    SJres = {5.40984e+06, 5.40984e+06, 1.91803e+07, 0, 0, 1.45804e-10};

    Vector TestMC, TestVM, TestDP, TestT, TestR, TestSJ;
    MohrCoulombCL.CalculateMaterialResponseCauchy(cl_parameters);
    TestMC = cl_parameters.GetStressVector();

    VonMisesCL.CalculateMaterialResponseCauchy(cl_parameters);
    TestVM = cl_parameters.GetStressVector();

    DruckerPragerCL.CalculateMaterialResponseCauchy(cl_parameters);
    TestDP = cl_parameters.GetStressVector();

    TrescaCL.CalculateMaterialResponseCauchy(cl_parameters);
    TestT = cl_parameters.GetStressVector();

    RankineCL.CalculateMaterialResponseCauchy(cl_parameters);
    TestR = cl_parameters.GetStressVector();

    material_properties.SetValue(FRACTURE_ENERGY, 1.0e5);
    cl_parameters.SetMaterialProperties(material_properties);
    SimoJuCL.CalculateMaterialResponseCauchy(cl_parameters);
    TestSJ = cl_parameters.GetStressVector();

    //Check the results
    KRATOS_CHECK_VECTOR_NEAR(MCres, TestMC, 0.0001e+06);
    KRATOS_CHECK_VECTOR_NEAR(VMres, TestVM, 0.0001e+06);
    KRATOS_CHECK_VECTOR_NEAR(DPres, TestDP, 0.001e+06);
    KRATOS_CHECK_VECTOR_NEAR(Tres, TestT, 0.0001e+06);
    KRATOS_CHECK_VECTOR_NEAR(Rres, TestR, 0.001e+06);
    KRATOS_CHECK_VECTOR_NEAR(SJres, TestSJ, 0.0001e+06);
}

KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawIntegrateStressDamageExponential, KratosStructuralMechanicsFastSuite)
{
    typedef GenericSmallStrainIsotropicDamage<GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> MC;
    typedef GenericSmallStrainIsotropicDamage<GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> VM;
    typedef GenericSmallStrainIsotropicDamage<GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> DP;
    typedef GenericSmallStrainIsotropicDamage<GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> T;
    typedef GenericSmallStrainIsotropicDamage<GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> SJ;
    typedef GenericSmallStrainIsotropicDamage<GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> R;

    ConstitutiveLaw::Parameters cl_parameters;
    Properties material_properties;
    Vector stress_vector, strain_vector;

    Model current_model;

    ModelPart& test_model_part = current_model.CreateModelPart("Main");

    NodeType::Pointer p_node_1 = test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    NodeType::Pointer p_node_2 = test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    NodeType::Pointer p_node_3 = test_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    NodeType::Pointer p_node_4 = test_model_part.CreateNewNode(4, 0.0, 0.0, 1.0);

    Tetrahedra3D4<NodeType> Geom = Tetrahedra3D4<NodeType>(p_node_1, p_node_2, p_node_3, p_node_4);

    stress_vector = ZeroVector(6);
    stress_vector[0] = 5.40984e+06;
    stress_vector[1] = 5.40984e+06;
    stress_vector[2] = 1.91803e+07;
    stress_vector[3] = 0.0;
    stress_vector[4] = 0.0;
    stress_vector[5] = 1.45804e-10;

    strain_vector = ZeroVector(6);
    strain_vector[0] = 0.0;
    strain_vector[1] = 0.0;
    strain_vector[2] = 8.0e-5;
    strain_vector[3] = 0.0;
    strain_vector[4] = 0.0;
    strain_vector[5] = 1.6941e-21;

    material_properties.SetValue(YOUNG_MODULUS, 210e9);
    material_properties.SetValue(POISSON_RATIO, 0.22);
    material_properties.SetValue(YIELD_STRESS_COMPRESSION, 3.0e6);
    material_properties.SetValue(YIELD_STRESS_TENSION, 3.0e6);
    material_properties.SetValue(FRICTION_ANGLE, 32.0);
    material_properties.SetValue(DILATANCY_ANGLE, 16.0);
    material_properties.SetValue(FRACTURE_ENERGY, 1.0e5);
    material_properties.SetValue(SOFTENING_TYPE, 1);

    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions=cl_parameters.GetOptions();
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
    MC MohrCoulombCL = MC();
    VM VonMisesCL = VM();
    DP DruckerPragerCL = DP();
    T TrescaCL = T();
    R RankineCL = R();
    SJ SimoJuCL = SJ();

    std::vector<double> MCres, VMres, DPres, Tres, Rres, SJres;
    MCres = {1.17707e+06,1.17707e+06,4.17324e+06,0,0,3.17239e-11};
    VMres = {1.17707e+06,1.17707e+06,4.17324e+06,0,0,3.17239e-11};
    DPres = {868915,868915,3.0807e+06,0,0,2.34186e-11};
    Tres  = {1.17707e+06,1.17707e+06,4.17324e+06,0,0,3.17239e-11};
    Rres  = {844533,844533,2.99425e+06,0,0,2.27615e-11};
    SJres = {870986,870986,3.08804e+06,0,0,2.34744e-11};

    Vector TestMC, TestVM, TestDP, TestT, TestR, TestSJ;
    MohrCoulombCL.CalculateMaterialResponseCauchy(cl_parameters);
    TestMC = cl_parameters.GetStressVector();

    VonMisesCL.CalculateMaterialResponseCauchy(cl_parameters);
    TestVM = cl_parameters.GetStressVector();

    DruckerPragerCL.CalculateMaterialResponseCauchy(cl_parameters);
    TestDP = cl_parameters.GetStressVector();

    TrescaCL.CalculateMaterialResponseCauchy(cl_parameters);
    TestT = cl_parameters.GetStressVector();

    RankineCL.CalculateMaterialResponseCauchy(cl_parameters);
    TestR = cl_parameters.GetStressVector();

    material_properties.SetValue(FRACTURE_ENERGY, 1.0e15);
    cl_parameters.SetMaterialProperties(material_properties);
    SimoJuCL.CalculateMaterialResponseCauchy(cl_parameters);
    TestSJ = cl_parameters.GetStressVector();

    //Check the results
    KRATOS_CHECK_VECTOR_NEAR(MCres, TestMC, 0.00001e+06);
    KRATOS_CHECK_VECTOR_NEAR(VMres, TestVM, 0.00001e+06);
    KRATOS_CHECK_VECTOR_NEAR(DPres, TestDP, 0.00001e+06);
    KRATOS_CHECK_VECTOR_NEAR(Tres, TestT, 0.00001e+06);
    KRATOS_CHECK_VECTOR_NEAR(Rres, TestR, 0.00001e+06);
    KRATOS_CHECK_VECTOR_NEAR(SJres, TestSJ, 0.00001e+06);
}
} // namespace Testing
} // namespace Kratos
