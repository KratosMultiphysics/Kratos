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

// Integrator
#include "custom_constitutive/constitutive_laws_integrators/generic_constitutive_law_integrator_kinematic_plasticity.h"

// Yield surfaces
#include "custom_constitutive/yield_surfaces/generic_yield_surface.h"
#include "custom_constitutive/yield_surfaces/von_mises_yield_surface.h"

// Plastic potentials
#include "custom_constitutive/plastic_potentials/generic_plastic_potential.h"
#include "custom_constitutive/plastic_potentials/von_mises_plastic_potential.h"

// Constitutive law
#include "custom_constitutive/generic_small_strain_kinematic_plasticity.h"
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
KRATOS_TEST_CASE_IN_SUITE(ConstitutiveLawIntegrateStressPlasticitySmallStrainKinematic, KratosStructuralMechanicsFastSuite)
{
    typedef GenericSmallStrainKinematicPlasticity <GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>> VM;

    ConstitutiveLaw::Parameters cl_parameters;
    Properties material_properties;
    Vector stress_vector, strain_vector;

    Model current_model;
    ModelPart& r_test_model_part = current_model.CreateModelPart("Main");

    NodeType::Pointer p_node_1 = r_test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    NodeType::Pointer p_node_2 = r_test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    NodeType::Pointer p_node_3 = r_test_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
    NodeType::Pointer p_node_4 = r_test_model_part.CreateNewNode(4, 0.0, 0.0, 1.0);

    Tetrahedra3D4<NodeType> Geom = Tetrahedra3D4<NodeType>(p_node_1, p_node_2, p_node_3, p_node_4);

    stress_vector = ZeroVector(6);
    strain_vector = ZeroVector(6);
    strain_vector[0] = 0.0;
    strain_vector[1] = 0.0;
    strain_vector[2] = -1.1e-04;
    strain_vector[3] = 0.0;
    strain_vector[4] = 0.0;
    strain_vector[5] = 0.0;

    material_properties.SetValue(YOUNG_MODULUS, 206900000000.0);
    material_properties.SetValue(POISSON_RATIO, 0.29);
    material_properties.SetValue(YIELD_STRESS, 1.5e6);
    material_properties.SetValue(FRICTION_ANGLE, 32.0);
    material_properties.SetValue(DILATANCY_ANGLE, 16.0);
    material_properties.SetValue(SOFTENING_TYPE, 1);
    material_properties.SetValue(FRACTURE_ENERGY, 1.5e2);
    material_properties.SetValue(HARDENING_CURVE, 3);
    material_properties.SetValue(KINEMATIC_HARDENING_TYPE, 1);
    Vector kin_hard_parameters = ZeroVector(3);
    kin_hard_parameters[0] = 15.0e9;
    material_properties.SetValue(KINEMATIC_PLASTICITY_PARAMETERS, kin_hard_parameters);

    // Set constitutive law flags:
    Flags& r_constitutive_law_options=cl_parameters.GetOptions();
    r_constitutive_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

    cl_parameters.SetElementGeometry(Geom);
    cl_parameters.SetProcessInfo(r_test_model_part.GetProcessInfo());
    cl_parameters.SetMaterialProperties(material_properties);
    cl_parameters.SetStrainVector(strain_vector);
    cl_parameters.SetStressVector(stress_vector);
    Matrix const_matrix;
    cl_parameters.SetConstitutiveMatrix(const_matrix);

    // Create the CL
    VM VonMisesCL = VM();
    std::vector<double> VMres = {-1.81882e+07,-1.81882e+07,-1.78117e+07,0,0,0};

    double plastic_dissipation;
    Vector TestMC, TestVM, TestDP, TestT;

    VonMisesCL.CalculateMaterialResponseCauchy(cl_parameters);
    VonMisesCL.FinalizeMaterialResponseCauchy(cl_parameters);
    TestVM = cl_parameters.GetStressVector();
    VonMisesCL.GetValue(PLASTIC_DISSIPATION, plastic_dissipation);
    KRATOS_WARNING_IF("TestPlasticity", plastic_dissipation < 1.0e-12) << "VonMises:: This test is not in plastic range" << std::endl;
	

    // Check the results
    const double tolerance = 1.0e-4;
    for (std::size_t comp = 0; comp < 6; ++comp){
        KRATOS_CHECK_LESS_EQUAL(std::abs((VMres[comp] - TestVM[comp]) / VMres[comp]), tolerance);
    }
}

} // namespace Testing
} // namespace Kratos
