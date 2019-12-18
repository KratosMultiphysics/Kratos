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

// Iga utility
#include "custom_utilities/geometry_utilities/iga_curve_on_surface_utilities.h"


// Constitutive law

#include "geometries/geometry.h"

namespace Kratos
{
namespace Testing
{
// We test the associated plasticity Constitutive laws...
typedef Node<3> NodeType;

/**
* Check the correct calculation of the integrated stress with the CL's in small strain
*/
KRATOS_TEST_CASE_IN_SUITE(CalculateVariations, KratosIgaFastSuite)
{
	//Vector N;
	//Matrix DN_De;
	//array_1d<double, 2> tangents;

	//std::vector<NodeType::Pointer> node_vector;
 //   NodeType::Pointer p_node_1 = test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
 //   NodeType::Pointer p_node_2 = test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
 //   NodeType::Pointer p_node_3 = test_model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
 //   NodeType::Pointer p_node_4 = test_model_part.CreateNewNode(4, 0.0, 0.0, 1.0);

	//node_vector.push_back(p_node_1);
	//node_vector.push_back(p_node_2);
	//node_vector.push_back(p_node_3);
	//node_vector.push_back(p_node_4);

	//Geometry<NodeType> Geom = Geometry<NodeType>(node_vector);

	//IgaCurveOnSurfaceUtilities::CalculateVariations(Geom, )

 //   stress_vector = ZeroVector(6);
 //   strain_vector = ZeroVector(6);
 //   strain_vector[0] = 0.0;
 //   strain_vector[1] = 0.0;
 //   strain_vector[2] = 8.0e-5;
 //   strain_vector[3] = 0.0;
 //   strain_vector[4] = 0.0;
 //   strain_vector[5] = 1.6941e-21;

 //   material_properties.SetValue(YOUNG_MODULUS, 210e9);
 //   material_properties.SetValue(POISSON_RATIO, 0.22);
 //   material_properties.SetValue(YIELD_STRESS_COMPRESSION, 1.0e5);
 //   material_properties.SetValue(YIELD_STRESS_TENSION, 1.0e5);
 //   material_properties.SetValue(FRICTION_ANGLE, 32.0);
 //   material_properties.SetValue(DILATANCY_ANGLE, 16.0);
 //   material_properties.SetValue(SOFTENING_TYPE, 1);
 //   material_properties.SetValue(FRACTURE_ENERGY, 10.0);
 //   material_properties.SetValue(HARDENING_CURVE, 0);

 //   // Set constitutive law flags:
 //   Flags& ConstitutiveLawOptions=cl_parameters.GetOptions();
 //   ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
 //   ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
 //   ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

 //   cl_parameters.SetElementGeometry(Geom);
 //   cl_parameters.SetProcessInfo(test_model_part.GetProcessInfo());
 //   cl_parameters.SetMaterialProperties(material_properties);
 //   cl_parameters.SetStrainVector(strain_vector);
 //   cl_parameters.SetStressVector(stress_vector);
 //   Matrix const_matrix;
 //   cl_parameters.SetConstitutiveMatrix(const_matrix);

 //   // Create the CL's
 //   MC MohrCoulombCL = MC();
 //   VM VonMisesCL = VM();
 //   DP DruckerPragerCL = DP();
 //   T TrescaCL = T();

 //   std::vector<double> MCres, VMres, DPres, Tres;
 //   MCres = {9.96862e+06,9.96862e+06,1.00628e+07,0,0,9.96661e-13};
 //   VMres = {9.96935e+06,9.96935e+06,1.00613e+07,0,0,9.73593e-13};
 //   DPres = {102202,102202,102199,0,0,-3.67488e-17};
 //   Tres = {9.96935e+06,9.96935e+06,1.00613e+07,0,0,9.73593e-13};

 //   double plastic_dissipation;
 //   Vector TestMC, TestVM, TestDP, TestT;
 //   MohrCoulombCL.CalculateMaterialResponseCauchy(cl_parameters);
 //   MohrCoulombCL.FinalizeMaterialResponseCauchy(cl_parameters);
 //   TestMC = cl_parameters.GetStressVector();
 //   MohrCoulombCL.GetValue(PLASTIC_DISSIPATION, plastic_dissipation);
 //   KRATOS_WARNING_IF("TestPlasticity", plastic_dissipation < 1.0e-12) << "MohrCoulomb:: This test is not in plastic range" << std::endl;

 //   VonMisesCL.CalculateMaterialResponseCauchy(cl_parameters);
 //   VonMisesCL.FinalizeMaterialResponseCauchy(cl_parameters);
 //   TestVM = cl_parameters.GetStressVector();
 //   VonMisesCL.GetValue(PLASTIC_DISSIPATION, plastic_dissipation);
 //   KRATOS_WARNING_IF("TestPlasticity", plastic_dissipation < 1.0e-12) << "VonMises:: This test is not in plastic range" << std::endl;

 //   DruckerPragerCL.CalculateMaterialResponseCauchy(cl_parameters);
 //   DruckerPragerCL.FinalizeMaterialResponseCauchy(cl_parameters);
 //   TestDP = cl_parameters.GetStressVector();
 //   DruckerPragerCL.GetValue(PLASTIC_DISSIPATION, plastic_dissipation);
 //   KRATOS_WARNING_IF("TestPlasticity", plastic_dissipation < 1.0e-12) << "DruckerPrager:: This test is not in plastic range" << std::endl;

 //   TrescaCL.CalculateMaterialResponseCauchy(cl_parameters);
 //   TrescaCL.FinalizeMaterialResponseCauchy(cl_parameters);
 //   TestT = cl_parameters.GetStressVector();
 //   TrescaCL.GetValue(PLASTIC_DISSIPATION, plastic_dissipation);
 //   KRATOS_WARNING_IF("TestPlasticity", plastic_dissipation < 1.0e-12) << "Tresca:: This test is not in plastic range" << std::endl;

 //   // Check the results
 //   const double tolerance = 1.0e-4;
 //   for (std::size_t comp = 0; comp < 6; ++comp){
 //       KRATOS_CHECK(!std::isnan(TestMC[comp]));
 //       KRATOS_CHECK_LESS_EQUAL(std::abs((MCres[comp] - TestMC[comp])/MCres[comp]), tolerance);
 //       KRATOS_CHECK(!std::isnan(VMres[comp]));
 //       KRATOS_CHECK_LESS_EQUAL(std::abs((VMres[comp] - TestVM[comp])/VMres[comp]), tolerance);
 //       KRATOS_CHECK(!std::isnan(DPres[comp]));
 //       KRATOS_CHECK_LESS_EQUAL(std::abs((DPres[comp] - TestDP[comp])/DPres[comp]), tolerance);
 //       KRATOS_CHECK(!std::isnan(TestT[comp]));
 //       KRATOS_CHECK_LESS_EQUAL(std::abs((Tres[comp] - TestT[comp])/Tres[comp]), tolerance);
    //}
}
} // namespace Testing
} // namespace Kratos
