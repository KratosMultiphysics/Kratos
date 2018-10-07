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

// Application includes

// Hyperelastic behaviours
#include "custom_constitutive/hyper_elastic_isotropic_kirchhoff_3d.h"
#include "custom_constitutive/hyper_elastic_isotropic_neo_hookean_3d.h"

// Integrator
#include "custom_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_plasticity.h"
#include "custom_constitutive/constitutive_laws_integrators/generic_finite_strain_constitutive_law_integrator_plasticity.h"

/* Small strain*/
// Yield surfaces
#include "custom_constitutive/yield_surfaces/generic_yield_surface.h"
#include "custom_constitutive/yield_surfaces/von_mises_yield_surface.h"
#include "custom_constitutive/yield_surfaces/modified_mohr_coulomb_yield_surface.h"
#include "custom_constitutive/yield_surfaces/rankine_yield_surface.h"
#include "custom_constitutive/yield_surfaces/simo_ju_yield_surface.h"
#include "custom_constitutive/yield_surfaces/drucker_prager_yield_surface.h"
#include "custom_constitutive/yield_surfaces/tresca_yield_surface.h"

// Plastic potentials
#include "custom_constitutive/plastic_potentials/generic_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/von_mises_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/tresca_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/modified_mohr_coulomb_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/drucker_prager_plastic_potential.h"

/* Finite strain*/
// Yield surfaces
#include "custom_constitutive/yield_surfaces/finite_strain/generic_yield_surface.h"
#include "custom_constitutive/yield_surfaces/finite_strain/von_mises_yield_surface.h"
#include "custom_constitutive/yield_surfaces/finite_strain/modified_mohr_coulomb_yield_surface.h"
#include "custom_constitutive/yield_surfaces/finite_strain/rankine_yield_surface.h"
#include "custom_constitutive/yield_surfaces/finite_strain/simo_ju_yield_surface.h"
#include "custom_constitutive/yield_surfaces/finite_strain/drucker_prager_yield_surface.h"
#include "custom_constitutive/yield_surfaces/finite_strain/tresca_yield_surface.h"

// Plastic potentials
#include "custom_constitutive/plastic_potentials/finite_strain/generic_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/finite_strain/von_mises_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/finite_strain/tresca_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/finite_strain/modified_mohr_coulomb_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/finite_strain/drucker_prager_plastic_potential.h"

// Constitutive law
#include "custom_constitutive/generic_small_strain_isotropic_plasticity.h"
#include "custom_constitutive/generic_finite_strain_isotropic_plasticity.h"
#include "includes/model_part.h"
#include "geometries/tetrahedra_3d_4.h"

namespace Kratos
{
namespace Testing
{
// We test the associated plasticity Constitutive laws...
typedef Node<3> NodeType;

/**
* Check the correct calculation of the integrated stress with the CL's in small strain
*/
KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawIntegrateStressPlasticitySmallStrain, KratosStructuralMechanicsFastSuite)
{
    typedef GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> MC;
    typedef GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> VM;
    typedef GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>> DP;
    typedef GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>> T;

    ConstitutiveLaw::Parameters cl_parameters;
    Properties material_properties;
    Vector stress_vector, strain_vector;

    ModelPart test_model_part("Main");

    NodeType::Pointer p_node_1 = test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    NodeType::Pointer p_node_2 = test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    NodeType::Pointer p_node_3 = test_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    NodeType::Pointer p_node_4 = test_model_part.CreateNewNode(4, 0.0, 0.0, 1.0);

    Tetrahedra3D4<NodeType> Geom = Tetrahedra3D4<NodeType>(p_node_1, p_node_2, p_node_3, p_node_4);

    stress_vector = ZeroVector(6);
    strain_vector = ZeroVector(6);
    strain_vector[0] = 0.0;
    strain_vector[1] = 0.0;
    strain_vector[2] = -8.0e-5;
    strain_vector[3] = 0.0;
    strain_vector[4] = 0.0;
    strain_vector[5] = -1.6941e-21;

    material_properties.SetValue(YOUNG_MODULUS, 210e9);
    material_properties.SetValue(POISSON_RATIO, 0.22);
    material_properties.SetValue(YIELD_STRESS_COMPRESSION, 3.0e6);
    material_properties.SetValue(YIELD_STRESS_TENSION, 3.0e6);
    material_properties.SetValue(FRICTION_ANGLE, 32.0);
    material_properties.SetValue(DILATANCY_ANGLE, 16.0);
    material_properties.SetValue(SOFTENING_TYPE, 1);
    material_properties.SetValue(FRACTURE_ENERGY, 1000.0);
    material_properties.SetValue(HARDENING_CURVE, 0);

    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions=cl_parameters.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

    cl_parameters.SetElementGeometry(Geom);
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

    std::vector<double> MCres, VMres, DPres, Tres;
    MCres = {-9.07094e+06, -9.07094e+06, -1.18581e+07, 0, 0, -2.95111e-11};
    VMres = {-9.09506e+06, -9.09506e+06, -1.18099e+07, 0, 0, -2.87449e-11};
    DPres = {-5.40984e+06, -5.40984e+06, -1.91803e+07, 0, 0, -1.45804e-10};
    Tres = {-9.09506e+06, -9.09506e+06, -1.18099e+07, 0, 0, -2.87449e-11};

    Vector TestMC, TestVM, TestDP, TestT;
    MohrCoulombCL.CalculateMaterialResponseCauchy(cl_parameters);
    TestMC = cl_parameters.GetStressVector();

    VonMisesCL.CalculateMaterialResponseCauchy(cl_parameters);
    TestVM = cl_parameters.GetStressVector();

    DruckerPragerCL.CalculateMaterialResponseCauchy(cl_parameters);
    TestDP = cl_parameters.GetStressVector();

    TrescaCL.CalculateMaterialResponseCauchy(cl_parameters);
    TestT = cl_parameters.GetStressVector();

    //Check the results
    for (int comp = 0; comp < 6; comp++)
    {
        KRATOS_CHECK_NEAR(MCres[comp], TestMC[comp], 0.0001e6);
        KRATOS_CHECK_NEAR(VMres[comp], TestVM[comp], 0.0001e6);
        KRATOS_CHECK_NEAR(DPres[comp], TestDP[comp], 0.0001e6);
        KRATOS_CHECK_NEAR(Tres[comp], TestT[comp], 0.0001e6);
    }
}

/**
* Check the correct calculation of the integrated stress with the CL's in finite strain
*/
KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawIntegrateStressPlasticityFiniteStrainKirchoff, KratosStructuralMechanicsFastSuite)
{
    typedef GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<FiniteStrainModifiedMohrCoulombYieldSurface<FiniteStrainModifiedMohrCoulombPlasticPotential<6>>>> MC;
    typedef GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<FiniteStrainVonMisesYieldSurface<FiniteStrainVonMisesPlasticPotential<6>>>> VM;
    typedef GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<FiniteStrainDruckerPragerYieldSurface<FiniteStrainDruckerPragerPlasticPotential<6>>>> DP;
    typedef GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<FiniteStrainTrescaYieldSurface<FiniteStrainTrescaPlasticPotential<6>>>> T;

    ConstitutiveLaw::Parameters cl_parameters;
    Properties material_properties;
    Vector stress_vector, strain_vector;

    ModelPart test_model_part("Main");

    NodeType::Pointer p_node_1 = test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    NodeType::Pointer p_node_2 = test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    NodeType::Pointer p_node_3 = test_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    NodeType::Pointer p_node_4 = test_model_part.CreateNewNode(4, 0.0, 0.0, 1.0);

    Tetrahedra3D4<NodeType> Geom = Tetrahedra3D4<NodeType>(p_node_1, p_node_2, p_node_3, p_node_4);

    stress_vector = ZeroVector(6);
    strain_vector = ZeroVector(6);
//     strain_vector[0] = 0.0;
//     strain_vector[1] = 0.0;
//     strain_vector[2] = -8.0e-5;
//     strain_vector[3] = 0.0;
//     strain_vector[4] = 0.0;
//     strain_vector[5] = -1.6941e-21;
    Matrix deformation_gradient = ZeroMatrix(3, 3);
    deformation_gradient(0,0) = 1.0;
    deformation_gradient(1,1) = 1.0;
    deformation_gradient(2,2) = (1.0 - 8.0e-5);
    deformation_gradient(2,0) = -1.6941e-21;
    deformation_gradient(0,2) = -1.6941e-21;

    material_properties.SetValue(YOUNG_MODULUS, 210e9);
    material_properties.SetValue(POISSON_RATIO, 0.22);
    material_properties.SetValue(YIELD_STRESS_COMPRESSION, 3.0e6);
    material_properties.SetValue(YIELD_STRESS_TENSION, 3.0e6);
    material_properties.SetValue(FRICTION_ANGLE, 32.0);
    material_properties.SetValue(DILATANCY_ANGLE, 16.0);
    material_properties.SetValue(SOFTENING_TYPE, 1);
    material_properties.SetValue(FRACTURE_ENERGY, 1000.0);
    material_properties.SetValue(HARDENING_CURVE, 0);

    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions=cl_parameters.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    cl_parameters.SetElementGeometry(Geom);
    cl_parameters.SetProcessInfo(test_model_part.GetProcessInfo());
    cl_parameters.SetMaterialProperties(material_properties);
    cl_parameters.SetStrainVector(strain_vector);
    cl_parameters.SetStressVector(stress_vector);
    cl_parameters.SetDeformationGradientF(deformation_gradient);
    cl_parameters.SetDeterminantF(MathUtils<double>::DetMat(deformation_gradient));
    Matrix const_matrix = ZeroMatrix(6, 6);
    cl_parameters.SetConstitutiveMatrix(const_matrix);

    // Create the CL's
    MC MohrCoulombCL = MC();
    VM VonMisesCL = VM();
    DP DruckerPragerCL = DP();
    T TrescaCL = T();

    std::vector<double> MCres, VMres, DPres, Tres;
    MCres = {-9.07094e+06, -9.07094e+06, -1.18581e+07, 0, 0, -2.95111e-11};
    VMres = {-9.09506e+06, -9.09506e+06, -1.18099e+07, 0, 0, -2.87449e-11};
    DPres = {-5.40984e+06, -5.40984e+06, -1.91803e+07, 0, 0, -1.45804e-10};
    Tres = {-9.09506e+06, -9.09506e+06, -1.18099e+07, 0, 0, -2.87449e-11};

    Vector TestMC, TestVM, TestDP, TestT;
    ConstitutiveLaw::Parameters cl_parameters_MC(cl_parameters);
    MohrCoulombCL.CalculateMaterialResponsePK2(cl_parameters_MC);
    TestMC = cl_parameters_MC.GetStressVector();

    ConstitutiveLaw::Parameters cl_parameters_CL(cl_parameters);
    VonMisesCL.CalculateMaterialResponsePK2(cl_parameters_CL);
    TestVM = cl_parameters_CL.GetStressVector();

    ConstitutiveLaw::Parameters cl_parameters_DP(cl_parameters);
    DruckerPragerCL.CalculateMaterialResponsePK2(cl_parameters_DP);
    TestDP = cl_parameters_DP.GetStressVector();

    ConstitutiveLaw::Parameters cl_parameters_T(cl_parameters);
    TrescaCL.CalculateMaterialResponsePK2(cl_parameters_T);
    TestT = cl_parameters_T.GetStressVector();

    //Check the results
    for (int comp = 0; comp < 3; comp++)
    {
        KRATOS_CHECK_NEAR(MCres[comp], TestMC[comp], 0.0001e6);
        KRATOS_CHECK_NEAR(VMres[comp], TestVM[comp], 0.0001e6);
        KRATOS_CHECK_NEAR(DPres[comp], TestDP[comp], 0.0001e6);
        KRATOS_CHECK_NEAR(Tres[comp], TestT[comp], 0.0001e6);
    }
}
} // namespace Testing
} // namespace Kratos
