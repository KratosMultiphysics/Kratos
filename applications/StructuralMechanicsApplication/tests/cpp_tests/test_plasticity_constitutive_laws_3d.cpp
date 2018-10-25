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

// Hyperelastic behaviours
#include "custom_constitutive/hyper_elastic_isotropic_kirchhoff_3d.h"
#include "custom_constitutive/hyper_elastic_isotropic_neo_hookean_3d.h"

// Integrator
#include "custom_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_plasticity.h"
#include "custom_constitutive/constitutive_laws_integrators/generic_finite_strain_constitutive_law_integrator_plasticity.h"

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

    Model current_model;
    ModelPart& test_model_part = current_model.CreateModelPart("Main");

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

    std::vector<double> MCres, VMres, DPres, Tres;
    MCres = {-9.07094e+06, -9.07094e+06, -1.18581e+07, 0, 0, -2.95111e-11};
    VMres = {-9.09506e+06, -9.09506e+06, -1.18099e+07, 0, 0, -2.87449e-11};
    DPres = {-5.40984e+06, -5.40984e+06, -1.91803e+07, 0, 0, -1.45804e-10};
    Tres = {-9.09506e+06, -9.09506e+06, -1.18099e+07, 0, 0, -2.87449e-11};

    double plastic_dissipation;
    Vector TestMC, TestVM, TestDP, TestT;
    MohrCoulombCL.CalculateMaterialResponseCauchy(cl_parameters);
    TestMC = cl_parameters.GetStressVector();
    MohrCoulombCL.GetValue(PLASTIC_DISSIPATION, plastic_dissipation);
    KRATOS_WARNING_IF("TestPlasticity", plastic_dissipation < 1.0e-12) << "MohrCoulomb:: This test is not in plastic range" << std::endl;

    VonMisesCL.CalculateMaterialResponseCauchy(cl_parameters);
    TestVM = cl_parameters.GetStressVector();
    VonMisesCL.GetValue(PLASTIC_DISSIPATION, plastic_dissipation);
    KRATOS_WARNING_IF("TestPlasticity", plastic_dissipation < 1.0e-12) << "VonMises:: This test is not in plastic range" << std::endl;

    DruckerPragerCL.CalculateMaterialResponseCauchy(cl_parameters);
    TestDP = cl_parameters.GetStressVector();
    DruckerPragerCL.GetValue(PLASTIC_DISSIPATION, plastic_dissipation);
    KRATOS_WARNING_IF("TestPlasticity", plastic_dissipation < 1.0e-12) << "DruckerPrager:: This test is not in plastic range" << std::endl;

    TrescaCL.CalculateMaterialResponseCauchy(cl_parameters);
    TestT = cl_parameters.GetStressVector();
    TrescaCL.GetValue(PLASTIC_DISSIPATION, plastic_dissipation);
    KRATOS_WARNING_IF("TestPlasticity", plastic_dissipation < 1.0e-12) << "Tresca:: This test is not in plastic range" << std::endl;

    // Check the results
    const double tolerance = 1.0e-4;
    for (std::size_t comp = 0; comp < 6; ++comp){
        KRATOS_CHECK(!std::isnan(TestMC[comp]));
        KRATOS_CHECK_LESS_EQUAL(std::abs((MCres[comp] - TestMC[comp])/MCres[comp]), tolerance);
        KRATOS_CHECK(!std::isnan(VMres[comp]));
        KRATOS_CHECK_LESS_EQUAL(std::abs((VMres[comp] - TestVM[comp])/VMres[comp]), tolerance);
        KRATOS_CHECK(!std::isnan(DPres[comp]));
        KRATOS_CHECK_LESS_EQUAL(std::abs((DPres[comp] - TestDP[comp])/DPres[comp]), tolerance);
        KRATOS_CHECK(!std::isnan(TestT[comp]));
        KRATOS_CHECK_LESS_EQUAL(std::abs((Tres[comp] - TestT[comp])/Tres[comp]), tolerance);
    }
}

/**
* Check the correct calculation of the integrated stress with the CL's in small strain
*/
KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawCTensorPlasticitySmallStrain, KratosStructuralMechanicsFastSuite)
{
    typedef GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> MC;
    typedef GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> VM;
    typedef GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>> DP;
    typedef GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>> T;

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
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

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

    Matrix MCres = ZeroMatrix(6, 6);
    MCres(0, 1) = 6.7623e+10;
    MCres(0, 2) = 6.7623e+10;
    MCres(1, 0) = 6.7623e+10;
    MCres(1, 1) = 2.39754e+11;
    MCres(1, 2) = 6.7623e+10;
    MCres(2, 0) = 6.7623e+10;
    MCres(2, 1) = 6.7623e+10;
    MCres(2, 2) = 2.39754e+11;
    MCres(3, 3) = 8.60656e+10;
    MCres(4, 4) = 8.60656e+10;

    Matrix VMres = ZeroMatrix(6, 6);
    VMres(0, 0) = 2.39754e+11;
    VMres(0, 1) = 6.7623e+10;
    VMres(0, 2) = 6.7623e+10;
    VMres(1, 0) = 6.7623e+10;
    VMres(1, 1) = 2.39754e+11;
    VMres(1, 2) = 6.7623e+10;
    VMres(2, 0) = 6.7623e+10;
    VMres(2, 1) = 6.7623e+10;
    VMres(2, 2) = 2.39754e+11;
    VMres(3, 3) = 8.60656e+10;
    VMres(4, 4) = 8.60656e+10;

    Matrix DPres = ZeroMatrix(6, 6);
    DPres(0, 0) = 2.39754e+11;
    DPres(0, 1) = 6.7623e+10;
    DPres(0, 2) = 6.7623e+10;
    DPres(1, 0) = 6.7623e+10;
    DPres(1, 1) = 2.39754e+11;
    DPres(1, 2) = 6.7623e+10;
    DPres(2, 0) = 6.7623e+10;
    DPres(2, 1) = 6.7623e+10;
    DPres(2, 2) = 2.39754e+11;
    DPres(3, 3) = 8.60656e+10;
    DPres(4, 4) = 8.60656e+10;
    DPres(5, 5) = 8.60656e+10;

    Matrix Tres = ZeroMatrix(6, 6);
    Tres(0, 0) = 2.39754e+11;
    Tres(0, 1) = 6.7623e+10;
    Tres(0, 2) = 6.7623e+10;
    Tres(1, 0) = 6.7623e+10;
    Tres(1, 1) = 2.39754e+11;
    Tres(1, 2) = 6.7623e+10;
    Tres(2, 0) = 6.7623e+10;
    Tres(2, 1) = 6.7623e+10;
    Tres(2, 2) = 2.39754e+11;
    Tres(3, 3) = 8.60656e+10;
    Tres(4, 4) = 8.60656e+10;

    double plastic_dissipation;
    Matrix TestMC, TestVM, TestDP, TestT;

    // Mohr Coulomb
    ConstitutiveLaw::Parameters cl_parameters_MC(cl_parameters);
    MohrCoulombCL.CalculateMaterialResponseCauchy(cl_parameters_MC);
    TestMC = cl_parameters_MC.GetConstitutiveMatrix();
    MohrCoulombCL.GetValue(PLASTIC_DISSIPATION, plastic_dissipation);
    KRATOS_WARNING_IF("TestPlasticity", plastic_dissipation < 1.0e-12) << "MohrCoulomb:: This test is not in plastic range" << std::endl;


    // Von Mises
    ConstitutiveLaw::Parameters cl_parameters_VM(cl_parameters);
    VonMisesCL.CalculateMaterialResponseCauchy(cl_parameters_VM);
    TestVM = cl_parameters_VM.GetConstitutiveMatrix();
    VonMisesCL.GetValue(PLASTIC_DISSIPATION, plastic_dissipation);
    KRATOS_WARNING_IF("TestPlasticity", plastic_dissipation < 1.0e-12) << "VonMises:: This test is not in plastic range" << std::endl;


    // Drucker Pragger
    ConstitutiveLaw::Parameters cl_parameters_DP(cl_parameters);
    DruckerPragerCL.CalculateMaterialResponseCauchy(cl_parameters_DP);
    TestDP = cl_parameters_DP.GetConstitutiveMatrix();
    DruckerPragerCL.GetValue(PLASTIC_DISSIPATION, plastic_dissipation);
    KRATOS_WARNING_IF("TestPlasticity", plastic_dissipation < 1.0e-12) << "DruckerPrager:: This test is not in plastic range" << std::endl;


    // Tresca
    ConstitutiveLaw::Parameters cl_parameters_T(cl_parameters);
    TrescaCL.CalculateMaterialResponseCauchy(cl_parameters_T);
    TestT = cl_parameters_T.GetConstitutiveMatrix();
    TrescaCL.GetValue(PLASTIC_DISSIPATION, plastic_dissipation);
    KRATOS_WARNING_IF("TestPlasticity", plastic_dissipation < 1.0e-12) << "Tresca:: This test is not in plastic range" << std::endl;


    // Check the results
    const double tolerance = 1.0e-4;
    for (std::size_t i = 0; i < 6 ; i++) {
        for (std::size_t j = 0; j < 6 ; j++) {
            KRATOS_CHECK(!std::isnan(TestMC(i, j)));
            if (std::abs(MCres(i, j)) > 0.0)
                KRATOS_CHECK_LESS_EQUAL(std::abs((MCres(i, j) - TestMC(i, j))/MCres(i, j)), tolerance);
            KRATOS_CHECK(!std::isnan(VMres(i, j)));
            if (std::abs(VMres(i, j)) > 0.0)
                KRATOS_CHECK_LESS_EQUAL(std::abs((VMres(i, j) - TestVM(i, j))/VMres(i, j)), tolerance);
            KRATOS_CHECK(!std::isnan(DPres(i, j)));
            if (std::abs(DPres(i, j)) > 0.0)
                KRATOS_CHECK_LESS_EQUAL(std::abs((DPres(i, j) - TestDP(i, j))/DPres(i, j)), tolerance);
            KRATOS_CHECK(!std::isnan(TestT(i, j)));
            if (std::abs(Tres(i, j)) > 0.0)
                KRATOS_CHECK_LESS_EQUAL(std::abs((Tres(i, j) - TestT(i, j))/Tres(i, j)), tolerance);
        }
    }
}

/**
* Check the correct calculation of the integrated stress with the CL's in finite strain
*/
KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawIntegrateStressPlasticityFiniteStrainKirchoff, KratosStructuralMechanicsFastSuite)
{
    typedef GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> MC;
    typedef GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> VM;
    typedef GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>> DP;
    typedef GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>> T;

    typedef GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> MCFS;
    typedef GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> VMFS;
    typedef GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>> DPFS;
    typedef GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>> TFS;

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
    strain_vector = ZeroVector(6);
    strain_vector[0] = 0.0;
    strain_vector[1] = 0.0;
    strain_vector[2] = -8.0e-5;
    strain_vector[3] = 0.0;
    strain_vector[4] = 0.0;
    strain_vector[5] = -1.6941e-21;
    Matrix deformation_gradient = ZeroMatrix(3, 3);
    deformation_gradient(0,0) = 1.0;
    deformation_gradient(1,1) = 1.0;
    deformation_gradient(2,2) = (1.0 - 8.0e-5);
    deformation_gradient(2,0) = 0.5 * -1.6941e-21;
    deformation_gradient(0,2) = 0.5 * -1.6941e-21;

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
    cl_parameters.SetProcessInfo(test_model_part.GetProcessInfo());
    cl_parameters.SetMaterialProperties(material_properties);
    cl_parameters.SetStrainVector(strain_vector);
    cl_parameters.SetStressVector(stress_vector);
    cl_parameters.SetDeformationGradientF(deformation_gradient);
    cl_parameters.SetDeterminantF(MathUtils<double>::DetMat(deformation_gradient));
    Matrix const_matrix = ZeroMatrix(6, 6);
    cl_parameters.SetConstitutiveMatrix(const_matrix);

    // Create the CL's (Infinitesimal strain) (We will check if for the same results in finite strain)
    MC MohrCoulombCL = MC();
    VM VonMisesCL = VM();
    DP DruckerPragerCL = DP();
    T TrescaCL = T();

    Vector MCres, VMres, DPres, Tres;

    // Mohr Coulomb
    ConstitutiveLaw::Parameters cl_parameters_MC_res(cl_parameters);
    MohrCoulombCL.CalculateMaterialResponseCauchy(cl_parameters_MC_res);
    MCres = cl_parameters_MC_res.GetStressVector();

    // Von Mises
    ConstitutiveLaw::Parameters cl_parameters_VM_res(cl_parameters);
    VonMisesCL.CalculateMaterialResponseCauchy(cl_parameters_VM_res);
    VMres = cl_parameters_MC_res.GetStressVector();

    // Drucker Pragger
    ConstitutiveLaw::Parameters cl_parameters_DP_res(cl_parameters);
    DruckerPragerCL.CalculateMaterialResponseCauchy(cl_parameters_DP_res);
    DPres = cl_parameters_DP_res.GetStressVector();

    // Tresca
    ConstitutiveLaw::Parameters cl_parameters_T_res(cl_parameters);
    TrescaCL.CalculateMaterialResponseCauchy(cl_parameters_T_res);
    Tres = cl_parameters_T_res.GetStressVector();

    // Create the CL's (Finite strain)
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
    MCFS MohrCoulombCLFS = MCFS();
    VMFS VonMisesCLFS = VMFS();
    DPFS DruckerPragerCLFS = DPFS();
    TFS TrescaCLFS = TFS();

    double plastic_dissipation;
    Vector TestMC, TestVM, TestDP, TestT;

    // Mohr Coulomb
    ConstitutiveLaw::Parameters cl_parameters_MC(cl_parameters);
    MohrCoulombCLFS.CalculateMaterialResponsePK2(cl_parameters_MC);
    TestMC = cl_parameters_MC.GetStressVector();
    MohrCoulombCLFS.GetValue(PLASTIC_DISSIPATION, plastic_dissipation);
    KRATOS_WARNING_IF("TestPlasticity", plastic_dissipation < 1.0e-12) << "MohrCoulomb:: This test is not in plastic range" << std::endl;

    // Von Mises
    ConstitutiveLaw::Parameters cl_parameters_VM(cl_parameters);
    VonMisesCLFS.CalculateMaterialResponsePK2(cl_parameters_VM);
    TestVM = cl_parameters_VM.GetStressVector();
    VonMisesCLFS.GetValue(PLASTIC_DISSIPATION, plastic_dissipation);
    KRATOS_WARNING_IF("TestPlasticity", plastic_dissipation < 1.0e-12) << "VonMises:: This test is not in plastic range" << std::endl;

    // Drucker Pragger
    ConstitutiveLaw::Parameters cl_parameters_DP(cl_parameters);
    DruckerPragerCLFS.CalculateMaterialResponsePK2(cl_parameters_DP);
    TestDP = cl_parameters_DP.GetStressVector();
    DruckerPragerCLFS.GetValue(PLASTIC_DISSIPATION, plastic_dissipation);
    KRATOS_WARNING_IF("TestPlasticity", plastic_dissipation < 1.0e-12) << "DruckerPrager:: This test is not in plastic range" << std::endl;

    // Tresca
    ConstitutiveLaw::Parameters cl_parameters_T(cl_parameters);
    TrescaCLFS.CalculateMaterialResponsePK2(cl_parameters_T);
    TestT = cl_parameters_T.GetStressVector();
    TrescaCLFS.GetValue(PLASTIC_DISSIPATION, plastic_dissipation);
    KRATOS_WARNING_IF("TestPlasticity", plastic_dissipation < 1.0e-12) << "Tresca:: This test is not in plastic range" << std::endl;

    // Check the results
    const double tolerance = 5.0e-3;
    for (int comp = 0; comp < 6; ++comp) {
        KRATOS_CHECK(!std::isnan(TestMC[comp]));
        if (TestMC[comp] > 1.0e-3) KRATOS_CHECK_LESS_EQUAL(std::abs((MCres[comp] - TestMC[comp])/MCres[comp]), tolerance);
        KRATOS_CHECK(!std::isnan(VMres[comp]));
        if (TestVM[comp] > 1.0e-3) KRATOS_CHECK_LESS_EQUAL(std::abs((VMres[comp] - TestVM[comp])/VMres[comp]), tolerance);
        KRATOS_CHECK(!std::isnan(DPres[comp]));
        if (TestDP[comp] > 1.0e-3) KRATOS_CHECK_LESS_EQUAL(std::abs((DPres[comp] - TestDP[comp])/DPres[comp]), tolerance);
        KRATOS_CHECK(!std::isnan(TestT[comp]));
        if (TestT[comp] > 1.0e-3) KRATOS_CHECK_LESS_EQUAL(std::abs((Tres[comp] - TestT[comp])/Tres[comp]), tolerance);
    }
}

/**
* Check the correct calculation of the integrated stress with the CL's in finite strain
*/
KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawCTensorPlasticityFiniteStrainKirchoff, KratosStructuralMechanicsFastSuite)
{
    typedef GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> MC;
    typedef GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> VM;
    typedef GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>> DP;
    typedef GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>> T;

    typedef GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> MCFS;
    typedef GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> VMFS;
    typedef GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>> DPFS;
    typedef GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicKirchhoff3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>> TFS;

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
    strain_vector = ZeroVector(6);
    strain_vector[0] = 0.0;
    strain_vector[1] = 0.0;
    strain_vector[2] = -8.0e-5;
    strain_vector[3] = 0.0;
    strain_vector[4] = 0.0;
    strain_vector[5] = -1.6941e-21;
    Matrix deformation_gradient = ZeroMatrix(3, 3);
    deformation_gradient(0,0) = 1.0;
    deformation_gradient(1,1) = 1.0;
    deformation_gradient(2,2) = (1.0 - 8.0e-5);
    deformation_gradient(2,0) = 0.5 * -1.6941e-21;
    deformation_gradient(0,2) = 0.5 * -1.6941e-21;

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
    cl_parameters.SetProcessInfo(test_model_part.GetProcessInfo());
    cl_parameters.SetMaterialProperties(material_properties);
    cl_parameters.SetStrainVector(strain_vector);
    cl_parameters.SetStressVector(stress_vector);
    cl_parameters.SetDeformationGradientF(deformation_gradient);
    cl_parameters.SetDeterminantF(MathUtils<double>::DetMat(deformation_gradient));
    Matrix const_matrix = ZeroMatrix(6, 6);
    cl_parameters.SetConstitutiveMatrix(const_matrix);

    // Create the CL's (Infinitesimal strain) (We will check if for the same results in finite strain)
    MC MohrCoulombCL = MC();
    VM VonMisesCL = VM();
    DP DruckerPragerCL = DP();
    T TrescaCL = T();

    Matrix MCres, VMres, DPres, Tres;

    // Mohr Coulomb
    ConstitutiveLaw::Parameters cl_parameters_MC_res(cl_parameters);
    MohrCoulombCL.CalculateMaterialResponseCauchy(cl_parameters_MC_res);
    MCres = cl_parameters_MC_res.GetConstitutiveMatrix();

    // Von Mises
    ConstitutiveLaw::Parameters cl_parameters_VM_res(cl_parameters);
    VonMisesCL.CalculateMaterialResponseCauchy(cl_parameters_VM_res);
    VMres = cl_parameters_MC_res.GetConstitutiveMatrix();

    // Drucker Pragger
    ConstitutiveLaw::Parameters cl_parameters_DP_res(cl_parameters);
    DruckerPragerCL.CalculateMaterialResponseCauchy(cl_parameters_DP_res);
    DPres = cl_parameters_DP_res.GetConstitutiveMatrix();

    // Tresca
    ConstitutiveLaw::Parameters cl_parameters_T_res(cl_parameters);
    TrescaCL.CalculateMaterialResponseCauchy(cl_parameters_T_res);
    Tres = cl_parameters_T_res.GetConstitutiveMatrix();

    // Create the CL's (Finite strain)
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
    MCFS MohrCoulombCLFS = MCFS();
    VMFS VonMisesCLFS = VMFS();
    DPFS DruckerPragerCLFS = DPFS();
    TFS TrescaCLFS = TFS();

    double plastic_dissipation;
    Matrix TestMC, TestVM, TestDP, TestT;

    // Mohr Coulomb
    ConstitutiveLaw::Parameters cl_parameters_MC(cl_parameters);
    MohrCoulombCLFS.CalculateMaterialResponsePK2(cl_parameters_MC);
    TestMC = cl_parameters_MC.GetConstitutiveMatrix();
    MohrCoulombCLFS.GetValue(PLASTIC_DISSIPATION, plastic_dissipation);
    KRATOS_WARNING_IF("TestPlasticity", plastic_dissipation < 1.0e-12) << "MohrCoulomb:: This test is not in plastic range" << std::endl;

    // Von Mises
    ConstitutiveLaw::Parameters cl_parameters_VM(cl_parameters);
    VonMisesCLFS.CalculateMaterialResponsePK2(cl_parameters_VM);
    TestVM = cl_parameters_VM.GetConstitutiveMatrix();
    VonMisesCLFS.GetValue(PLASTIC_DISSIPATION, plastic_dissipation);
    KRATOS_WARNING_IF("TestPlasticity", plastic_dissipation < 1.0e-12) << "VonMises:: This test is not in plastic range" << std::endl;

    // Drucker Pragger
    ConstitutiveLaw::Parameters cl_parameters_DP(cl_parameters);
    DruckerPragerCLFS.CalculateMaterialResponsePK2(cl_parameters_DP);
    TestDP = cl_parameters_DP.GetConstitutiveMatrix();
    DruckerPragerCLFS.GetValue(PLASTIC_DISSIPATION, plastic_dissipation);
    KRATOS_WARNING_IF("TestPlasticity", plastic_dissipation < 1.0e-12) << "DruckerPrager:: This test is not in plastic range" << std::endl;

    // Tresca
    ConstitutiveLaw::Parameters cl_parameters_T(cl_parameters);
    TrescaCLFS.CalculateMaterialResponsePK2(cl_parameters_T);
    TestT = cl_parameters_T.GetConstitutiveMatrix();
    TrescaCLFS.GetValue(PLASTIC_DISSIPATION, plastic_dissipation);
    KRATOS_WARNING_IF("TestPlasticity", plastic_dissipation < 1.0e-12) << "Tresca:: This test is not in plastic range" << std::endl;

    // Check the results
    const double tolerance = 5.0e-3;
    for (std::size_t i = 0; i < 6 ; i++) {
        for (std::size_t j = 0; j < 6 ; j++) {
            KRATOS_CHECK(!std::isnan(TestMC(i, j)));
            if (std::abs(MCres(i, j)) > 0.0)
                KRATOS_CHECK_LESS_EQUAL(std::abs((MCres(i, j) - TestMC(i, j))/MCres(i, j)), tolerance);
            KRATOS_CHECK(!std::isnan(VMres(i, j)));
            if (std::abs(VMres(i, j)) > 0.0)
                KRATOS_CHECK_LESS_EQUAL(std::abs((VMres(i, j) - TestVM(i, j))/VMres(i, j)), tolerance);
            KRATOS_CHECK(!std::isnan(DPres(i, j)));
            if (std::abs(DPres(i, j)) > 0.0)
                KRATOS_CHECK_LESS_EQUAL(std::abs((DPres(i, j) - TestDP(i, j))/DPres(i, j)), tolerance);
            KRATOS_CHECK(!std::isnan(TestT(i, j)));
            if (std::abs(Tres(i, j)) > 0.0)
                KRATOS_CHECK_LESS_EQUAL(std::abs((Tres(i, j) - TestT(i, j))/Tres(i, j)), tolerance);
        }
    }
}

/**
* Check the correct calculation of the integrated stress with the CL's in finite strain
*/
KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawIntegrateStressPlasticityFiniteStrainNeoHookean, KratosStructuralMechanicsFastSuite)
{
    typedef GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> MC;
    typedef GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> VM;
    typedef GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>> DP;
    typedef GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>> T;

    typedef GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> MCFS;
    typedef GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> VMFS;
    typedef GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>> DPFS;
    typedef GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>> TFS;

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
    strain_vector = ZeroVector(6);
    strain_vector[0] = 0.0;
    strain_vector[1] = 0.0;
    strain_vector[2] = -8.0e-5;
    strain_vector[3] = 0.0;
    strain_vector[4] = 0.0;
    strain_vector[5] = -1.6941e-21;
    Matrix deformation_gradient = ZeroMatrix(3, 3);
    deformation_gradient(0,0) = 1.0;
    deformation_gradient(1,1) = 1.0;
    deformation_gradient(2,2) = (1.0 - 8.0e-5);
    deformation_gradient(2,0) = 0.5 * -1.6941e-21;
    deformation_gradient(0,2) = 0.5 * -1.6941e-21;

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
    cl_parameters.SetProcessInfo(test_model_part.GetProcessInfo());
    cl_parameters.SetMaterialProperties(material_properties);
    cl_parameters.SetStrainVector(strain_vector);
    cl_parameters.SetStressVector(stress_vector);
    cl_parameters.SetDeformationGradientF(deformation_gradient);
    cl_parameters.SetDeterminantF(MathUtils<double>::DetMat(deformation_gradient));
    Matrix const_matrix = ZeroMatrix(6, 6);
    cl_parameters.SetConstitutiveMatrix(const_matrix);

    // Create the CL's (Infinitesimal strain) (We will check if for the same results in finite strain)
    MC MohrCoulombCL = MC();
    VM VonMisesCL = VM();
    DP DruckerPragerCL = DP();
    T TrescaCL = T();

    Vector MCres, VMres, DPres, Tres;

    // Mohr Coulomb
    ConstitutiveLaw::Parameters cl_parameters_MC_res(cl_parameters);
    MohrCoulombCL.CalculateMaterialResponseCauchy(cl_parameters_MC_res);
    MCres = cl_parameters_MC_res.GetStressVector();

    // Von Mises
    ConstitutiveLaw::Parameters cl_parameters_VM_res(cl_parameters);
    VonMisesCL.CalculateMaterialResponseCauchy(cl_parameters_VM_res);
    VMres = cl_parameters_MC_res.GetStressVector();

    // Drucker Pragger
    ConstitutiveLaw::Parameters cl_parameters_DP_res(cl_parameters);
    DruckerPragerCL.CalculateMaterialResponseCauchy(cl_parameters_DP_res);
    DPres = cl_parameters_DP_res.GetStressVector();

    // Tresca
    ConstitutiveLaw::Parameters cl_parameters_T_res(cl_parameters);
    TrescaCL.CalculateMaterialResponseCauchy(cl_parameters_T_res);
    Tres = cl_parameters_T_res.GetStressVector();

    // Create the CL's (Finite strain)
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
    MCFS MohrCoulombCLFS = MCFS();
    VMFS VonMisesCLFS = VMFS();
    DPFS DruckerPragerCLFS = DPFS();
    TFS TrescaCLFS = TFS();

    double plastic_dissipation;
    Vector TestMC, TestVM, TestDP, TestT;

    // Mohr Coulomb
    ConstitutiveLaw::Parameters cl_parameters_MC(cl_parameters);
    MohrCoulombCLFS.CalculateMaterialResponsePK2(cl_parameters_MC);
    TestMC = cl_parameters_MC.GetStressVector();
    MohrCoulombCLFS.GetValue(PLASTIC_DISSIPATION, plastic_dissipation);
    KRATOS_WARNING_IF("TestPlasticity", plastic_dissipation < 1.0e-12) << "MohrCoulomb:: This test is not in plastic range" << std::endl;

    // Von Mises
    ConstitutiveLaw::Parameters cl_parameters_VM(cl_parameters);
    VonMisesCLFS.CalculateMaterialResponsePK2(cl_parameters_VM);
    TestVM = cl_parameters_VM.GetStressVector();
    VonMisesCLFS.GetValue(PLASTIC_DISSIPATION, plastic_dissipation);
    KRATOS_WARNING_IF("TestPlasticity", plastic_dissipation < 1.0e-12) << "VonMises:: This test is not in plastic range" << std::endl;

    // Drucker Pragger
    ConstitutiveLaw::Parameters cl_parameters_DP(cl_parameters);
    DruckerPragerCLFS.CalculateMaterialResponsePK2(cl_parameters_DP);
    TestDP = cl_parameters_DP.GetStressVector();
    DruckerPragerCLFS.GetValue(PLASTIC_DISSIPATION, plastic_dissipation);
    KRATOS_WARNING_IF("TestPlasticity", plastic_dissipation < 1.0e-12) << "DruckerPrager:: This test is not in plastic range" << std::endl;

    // Tresca
    ConstitutiveLaw::Parameters cl_parameters_T(cl_parameters);
    TrescaCLFS.CalculateMaterialResponsePK2(cl_parameters_T);
    TestT = cl_parameters_T.GetStressVector();
    TrescaCLFS.GetValue(PLASTIC_DISSIPATION, plastic_dissipation);
    KRATOS_WARNING_IF("TestPlasticity", plastic_dissipation < 1.0e-12) << "Tresca:: This test is not in plastic range" << std::endl;

    // Check the results
    const double tolerance = 5.0e-3;
    for (int comp = 0; comp < 6; ++comp) {
        KRATOS_CHECK(!std::isnan(TestMC[comp]));
        if (TestMC[comp] > 1.0e-3) KRATOS_CHECK_LESS_EQUAL(std::abs((MCres[comp] - TestMC[comp])/MCres[comp]), tolerance);
        KRATOS_CHECK(!std::isnan(VMres[comp]));
        if (TestVM[comp] > 1.0e-3) KRATOS_CHECK_LESS_EQUAL(std::abs((VMres[comp] - TestVM[comp])/VMres[comp]), tolerance);
        KRATOS_CHECK(!std::isnan(DPres[comp]));
        if (TestDP[comp] > 1.0e-3) KRATOS_CHECK_LESS_EQUAL(std::abs((DPres[comp] - TestDP[comp])/DPres[comp]), tolerance);
        KRATOS_CHECK(!std::isnan(TestT[comp]));
        if (TestT[comp] > 1.0e-3) KRATOS_CHECK_LESS_EQUAL(std::abs((Tres[comp] - TestT[comp])/Tres[comp]), tolerance);
    }
}

/**
* Check the correct calculation of the integrated stress with the CL's in finite strain
*/
KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawCTensorPlasticityFiniteStrainNeoHookean, KratosStructuralMechanicsFastSuite)
{
    typedef GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> MC;
    typedef GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> VM;
    typedef GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>> DP;
    typedef GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>> T;

    typedef GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>> MCFS;
    typedef GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> VMFS;
    typedef GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>> DPFS;
    typedef GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>> TFS;

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
    strain_vector = ZeroVector(6);
    strain_vector[0] = 0.0;
    strain_vector[1] = 0.0;
    strain_vector[2] = -8.0e-5;
    strain_vector[3] = 0.0;
    strain_vector[4] = 0.0;
    strain_vector[5] = -1.6941e-21;
    Matrix deformation_gradient = ZeroMatrix(3, 3);
    deformation_gradient(0,0) = 1.0;
    deformation_gradient(1,1) = 1.0;
    deformation_gradient(2,2) = (1.0 - 8.0e-5);
    deformation_gradient(2,0) = 0.5 * -1.6941e-21;
    deformation_gradient(0,2) = 0.5 * -1.6941e-21;

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
    cl_parameters.SetProcessInfo(test_model_part.GetProcessInfo());
    cl_parameters.SetMaterialProperties(material_properties);
    cl_parameters.SetStrainVector(strain_vector);
    cl_parameters.SetStressVector(stress_vector);
    cl_parameters.SetDeformationGradientF(deformation_gradient);
    cl_parameters.SetDeterminantF(MathUtils<double>::DetMat(deformation_gradient));
    Matrix const_matrix = ZeroMatrix(6, 6);
    cl_parameters.SetConstitutiveMatrix(const_matrix);

    // Create the CL's (Infinitesimal strain) (We will check if for the same results in finite strain)
    MC MohrCoulombCL = MC();
    VM VonMisesCL = VM();
    DP DruckerPragerCL = DP();
    T TrescaCL = T();

    Matrix MCres, VMres, DPres, Tres;

    // Mohr Coulomb
    ConstitutiveLaw::Parameters cl_parameters_MC_res(cl_parameters);
    MohrCoulombCL.CalculateMaterialResponseCauchy(cl_parameters_MC_res);
    MCres = cl_parameters_MC_res.GetConstitutiveMatrix();

    // Von Mises
    ConstitutiveLaw::Parameters cl_parameters_VM_res(cl_parameters);
    VonMisesCL.CalculateMaterialResponseCauchy(cl_parameters_VM_res);
    VMres = cl_parameters_MC_res.GetConstitutiveMatrix();

    // Drucker Pragger
    ConstitutiveLaw::Parameters cl_parameters_DP_res(cl_parameters);
    DruckerPragerCL.CalculateMaterialResponseCauchy(cl_parameters_DP_res);
    DPres = cl_parameters_DP_res.GetConstitutiveMatrix();

    // Tresca
    ConstitutiveLaw::Parameters cl_parameters_T_res(cl_parameters);
    TrescaCL.CalculateMaterialResponseCauchy(cl_parameters_T_res);
    Tres = cl_parameters_T_res.GetConstitutiveMatrix();

    // Create the CL's (Finite strain)
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
    MCFS MohrCoulombCLFS = MCFS();
    VMFS VonMisesCLFS = VMFS();
    DPFS DruckerPragerCLFS = DPFS();
    TFS TrescaCLFS = TFS();

    double plastic_dissipation;
    Matrix TestMC, TestVM, TestDP, TestT;

    // Mohr Coulomb
    ConstitutiveLaw::Parameters cl_parameters_MC(cl_parameters);
    MohrCoulombCLFS.CalculateMaterialResponsePK2(cl_parameters_MC);
    TestMC = cl_parameters_MC.GetConstitutiveMatrix();
    MohrCoulombCLFS.GetValue(PLASTIC_DISSIPATION, plastic_dissipation);
    KRATOS_WARNING_IF("TestPlasticity", plastic_dissipation < 1.0e-12) << "MohrCoulomb:: This test is not in plastic range" << std::endl;

    // Von Mises
    ConstitutiveLaw::Parameters cl_parameters_VM(cl_parameters);
    VonMisesCLFS.CalculateMaterialResponsePK2(cl_parameters_VM);
    TestVM = cl_parameters_VM.GetConstitutiveMatrix();
    VonMisesCLFS.GetValue(PLASTIC_DISSIPATION, plastic_dissipation);
    KRATOS_WARNING_IF("TestPlasticity", plastic_dissipation < 1.0e-12) << "VonMises:: This test is not in plastic range" << std::endl;

    // Drucker Pragger
    ConstitutiveLaw::Parameters cl_parameters_DP(cl_parameters);
    DruckerPragerCLFS.CalculateMaterialResponsePK2(cl_parameters_DP);
    TestDP = cl_parameters_DP.GetConstitutiveMatrix();
    DruckerPragerCLFS.GetValue(PLASTIC_DISSIPATION, plastic_dissipation);
    KRATOS_WARNING_IF("TestPlasticity", plastic_dissipation < 1.0e-12) << "DruckerPrager:: This test is not in plastic range" << std::endl;

    // Tresca
    ConstitutiveLaw::Parameters cl_parameters_T(cl_parameters);
    TrescaCLFS.CalculateMaterialResponsePK2(cl_parameters_T);
    TestT = cl_parameters_T.GetConstitutiveMatrix();
    TrescaCLFS.GetValue(PLASTIC_DISSIPATION, plastic_dissipation);
    KRATOS_WARNING_IF("TestPlasticity", plastic_dissipation < 1.0e-12) << "Tresca:: This test is not in plastic range" << std::endl;

    // Check the results
    const double tolerance = 5.0e-3;
    for (std::size_t i = 0; i < 6 ; i++) {
        for (std::size_t j = 0; j < 6 ; j++) {
            KRATOS_CHECK(!std::isnan(TestMC(i, j)));
            if (std::abs(MCres(i, j)) > 0.0)
                KRATOS_CHECK_LESS_EQUAL(std::abs((MCres(i, j) - TestMC(i, j))/MCres(i, j)), tolerance);
            KRATOS_CHECK(!std::isnan(VMres(i, j)));
            if (std::abs(VMres(i, j)) > 0.0)
                KRATOS_CHECK_LESS_EQUAL(std::abs((VMres(i, j) - TestVM(i, j))/VMres(i, j)), tolerance);
            KRATOS_CHECK(!std::isnan(DPres(i, j)));
            if (std::abs(DPres(i, j)) > 0.0)
                KRATOS_CHECK_LESS_EQUAL(std::abs((DPres(i, j) - TestDP(i, j))/DPres(i, j)), tolerance);
            KRATOS_CHECK(!std::isnan(TestT(i, j)));
            if (std::abs(Tres(i, j)) > 0.0)
                KRATOS_CHECK_LESS_EQUAL(std::abs((Tres(i, j) - TestT(i, j))/Tres(i, j)), tolerance);
        }
    }
}
} // namespace Testing
} // namespace Kratos
