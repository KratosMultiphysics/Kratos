//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Peter Wilson
//


// System includes

// External includes

// Project includes
#include "includes/process_info.h"
#include "testing/testing.h"
#include "containers/model.h"

// Application includes
#include "mpm_application_variables.h"

// Material law
#include "custom_constitutive/johnson_cook_thermal_plastic_plane_strain_2D_law.hpp"
#include "includes/model_part.h"
#include "geometries/quadrilateral_2d_4.h"


namespace Kratos
{
namespace Testing
{
    typedef Node NodeType;

    KRATOS_TEST_CASE_IN_SUITE(MPMConstitutiveLawJohnsonCookWithThermalSoftening, KratosMPMFastSuite)
    {
        ConstitutiveLaw::Parameters cl_parameters;
        Properties material_properties;
        Vector stress_vector = ZeroVector(3);
        Vector strain_vector(3);

        // Create gauss point
        Model current_model;
        ModelPart& test_model_part = current_model.CreateModelPart("Main");
        NodeType::Pointer p_node_1 = test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        NodeType::Pointer p_node_2 = test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
        NodeType::Pointer p_node_3 = test_model_part.CreateNewNode(3, 1.0, 1.0, 0.0);
        NodeType::Pointer p_node_4 = test_model_part.CreateNewNode(4, 0.0, 1.0, 0.0);
        Quadrilateral2D4<NodeType> geometry = Quadrilateral2D4<NodeType>(p_node_1, p_node_2, p_node_3, p_node_4);

        // Set material properties
        material_properties.SetValue(DENSITY, 7830);
        material_properties.SetValue(YOUNG_MODULUS, 200e9);
        material_properties.SetValue(POISSON_RATIO, 0.29);
        material_properties.SetValue(TEMPERATURE, 294.0);
        material_properties.SetValue(JC_PARAMETER_A, 792e6);
        material_properties.SetValue(JC_PARAMETER_B, 510e6);
        material_properties.SetValue(JC_PARAMETER_C, 0.014);
        material_properties.SetValue(JC_PARAMETER_m, 1.03);
        material_properties.SetValue(JC_PARAMETER_n, 0.26);
        material_properties.SetValue(REFERENCE_STRAIN_RATE, 1.0);
        material_properties.SetValue(REFERENCE_TEMPERATURE, 294.0);
        material_properties.SetValue(MELD_TEMPERATURE, 1793.0);
        material_properties.SetValue(SPECIFIC_HEAT, 477.0);
        material_properties.SetValue(TAYLOR_QUINNEY_COEFFICIENT, 0.9);

        test_model_part.GetProcessInfo().SetValue(DELTA_TIME, 0.001);
        test_model_part.GetProcessInfo().SetValue(IS_EXPLICIT, true);

        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions = cl_parameters.GetOptions();
        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

        const ProcessInfo& r_current_process_info = test_model_part.GetProcessInfo();

        // Set required constitutive law parameters:
        cl_parameters.SetElementGeometry(geometry);
        cl_parameters.SetProcessInfo(r_current_process_info);
        cl_parameters.SetMaterialProperties(material_properties);
        cl_parameters.SetStrainVector(strain_vector);
        cl_parameters.SetStressVector(stress_vector);

        Matrix F = IdentityMatrix(2);
        double detF = 1.0;
        cl_parameters.SetDeformationGradientF(F);
        cl_parameters.SetDeterminantF(detF);
        Matrix dummy_matrix;
        Vector dummy_vector;
        cl_parameters.SetConstitutiveMatrix(dummy_matrix);
        cl_parameters.SetShapeFunctionsValues(dummy_vector);
        cl_parameters.SetShapeFunctionsDerivatives(dummy_matrix);

        // Create the CL
        JohnsonCookThermalPlastic2DPlaneStrainLaw cl = JohnsonCookThermalPlastic2DPlaneStrainLaw();

        // Set variables for the test
        const double tolerance = 1.0e-4;
        double ref_temperature = 294.22441567308590;
        double ref_plastic_strain = 0.0011138185110505848;
        double ref_plastic_strain_rate = 1.1138185110505847;
        double ref_equivalent_stress = 880269994.07343519;
        strain_vector[0] = 0.004;
        strain_vector[1] = 0.002;
        strain_vector[2] = 0.008;

        // Run test
        Vector dummy;
        double value;
        cl.InitializeMaterial(material_properties, geometry, dummy);
        cl.CalculateMaterialResponseKirchhoff(cl_parameters);

        // Check MP_TEMPERATURE
        cl.GetValue(MP_TEMPERATURE, value);
        KRATOS_EXPECT_NEAR(ref_temperature, value, tolerance);

        // Check MP_EQUIVALENT_PLASTIC_STRAIN
        cl.GetValue(MP_EQUIVALENT_PLASTIC_STRAIN, value);
        KRATOS_EXPECT_NEAR(ref_plastic_strain, value, tolerance*tolerance);

        // Check MP_EQUIVALENT_PLASTIC_STRAIN_RATE
        cl.GetValue(MP_EQUIVALENT_PLASTIC_STRAIN_RATE, value);
        KRATOS_EXPECT_NEAR(ref_plastic_strain_rate, value, tolerance*tolerance);

        // Check MP_EQUIVALENT_STRESS
        cl.GetValue(MP_EQUIVALENT_STRESS, value);
        KRATOS_EXPECT_NEAR(ref_equivalent_stress, value, tolerance);
    }

    KRATOS_TEST_CASE_IN_SUITE(MPMConstitutiveLawJohnsonCookWithoutThermalSoftening, KratosMPMFastSuite)
    {
        ConstitutiveLaw::Parameters cl_parameters;
        Properties material_properties;
        Vector stress_vector = ZeroVector(3);
        Vector strain_vector(3);

        // Create gauss point
        Model current_model;
        ModelPart& test_model_part = current_model.CreateModelPart("Main");
        NodeType::Pointer p_node_1 = test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        NodeType::Pointer p_node_2 = test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
        NodeType::Pointer p_node_3 = test_model_part.CreateNewNode(3, 1.0, 1.0, 0.0);
        NodeType::Pointer p_node_4 = test_model_part.CreateNewNode(4, 0.0, 1.0, 0.0);
        Quadrilateral2D4<NodeType> geometry = Quadrilateral2D4<NodeType>(p_node_1, p_node_2, p_node_3, p_node_4);

        // Set material properties
        material_properties.SetValue(DENSITY, 7830);
        material_properties.SetValue(YOUNG_MODULUS, 200e9);
        material_properties.SetValue(POISSON_RATIO, 0.29);
        material_properties.SetValue(TEMPERATURE, 294.0);
        material_properties.SetValue(JC_PARAMETER_A, 792e6);
        material_properties.SetValue(JC_PARAMETER_B, 510e6);
        material_properties.SetValue(JC_PARAMETER_C, 0.014);
        material_properties.SetValue(JC_PARAMETER_m, 1.03);
        material_properties.SetValue(JC_PARAMETER_n, 0.26);
        material_properties.SetValue(REFERENCE_STRAIN_RATE, 1.0);
        material_properties.SetValue(REFERENCE_TEMPERATURE, 294.0);
        material_properties.SetValue(MELD_TEMPERATURE, 1793.0);
        material_properties.SetValue(SPECIFIC_HEAT, 477.0);
        material_properties.SetValue(TAYLOR_QUINNEY_COEFFICIENT, 0.0);

        test_model_part.GetProcessInfo().SetValue(DELTA_TIME, 0.001);
        test_model_part.GetProcessInfo().SetValue(IS_EXPLICIT, true);

        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions = cl_parameters.GetOptions();
        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

        const ProcessInfo& r_current_process_info = test_model_part.GetProcessInfo();

        // Set required constitutive law parameters:
        cl_parameters.SetElementGeometry(geometry);
        cl_parameters.SetProcessInfo(r_current_process_info);
        cl_parameters.SetMaterialProperties(material_properties);
        cl_parameters.SetStrainVector(strain_vector);
        cl_parameters.SetStressVector(stress_vector);

        Matrix F = IdentityMatrix(2);
        double detF = 1.0;
        cl_parameters.SetDeformationGradientF(F);
        cl_parameters.SetDeterminantF(detF);
        Matrix dummy_matrix;
        Vector dummy_vector;
        cl_parameters.SetConstitutiveMatrix(dummy_matrix);
        cl_parameters.SetShapeFunctionsValues(dummy_vector);
        cl_parameters.SetShapeFunctionsDerivatives(dummy_matrix);

        // Create the CL
        JohnsonCookThermalPlastic2DPlaneStrainLaw cl = JohnsonCookThermalPlastic2DPlaneStrainLaw();

        // Set variables for the test
        const double tolerance = 1.0e-4;
        double ref_temperature = 294.0;
        double ref_plastic_strain = 0.0011134350919909683;
        double ref_plastic_strain_rate = 1.1134350919909684;
        double ref_equivalent_stress = 880359161.29660177;
        strain_vector[0] = 0.004;
        strain_vector[1] = 0.002;
        strain_vector[2] = 0.008;

        // Run test
        Vector dummy;
        double value;
        cl.InitializeMaterial(material_properties, geometry, dummy);
        cl.CalculateMaterialResponseKirchhoff(cl_parameters);

        // Check MP_TEMPERATURE
        cl.GetValue(MP_TEMPERATURE, value);
        KRATOS_EXPECT_NEAR(ref_temperature, value, tolerance);

        // Check MP_EQUIVALENT_PLASTIC_STRAIN
        cl.GetValue(MP_EQUIVALENT_PLASTIC_STRAIN, value);
        KRATOS_EXPECT_NEAR(ref_plastic_strain, value, tolerance * tolerance);

        // Check MP_EQUIVALENT_PLASTIC_STRAIN_RATE
        cl.GetValue(MP_EQUIVALENT_PLASTIC_STRAIN_RATE, value);
        KRATOS_EXPECT_NEAR(ref_plastic_strain_rate, value, tolerance * tolerance);

        // Check MP_EQUIVALENT_STRESS
        cl.GetValue(MP_EQUIVALENT_STRESS, value);
        KRATOS_EXPECT_NEAR(ref_equivalent_stress, value, tolerance);
    }

} // namespace Testing
} // namespace Kratos
