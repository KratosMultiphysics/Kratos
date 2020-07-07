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
#include "particle_mechanics_application_variables.h"

// Material law
#include "custom_constitutive/rht_concrete_3D_law.hpp"
#include "includes/model_part.h"
#include "geometries/hexahedra_3d_8.h"


namespace Kratos
{
namespace Testing
{
    typedef Node<3> NodeType;

    KRATOS_TEST_CASE_IN_SUITE(ParticleConstitutiveLawRHTConcrete3D, KratosParticleMechanicsFastSuite)
    {
        const bool print_results = true;

        ConstitutiveLaw::Parameters cl_parameters;
        Properties material_properties;
        Vector stress_vector = ZeroVector(6);
        Vector strain_vector = ZeroVector(6);

        // Create gauss point
        Model current_model;
        ModelPart& test_model_part = current_model.CreateModelPart("Main");
        NodeType::Pointer p_node_1 = test_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        NodeType::Pointer p_node_2 = test_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
        NodeType::Pointer p_node_3 = test_model_part.CreateNewNode(3, 1.0, 1.0, 0.0);
        NodeType::Pointer p_node_4 = test_model_part.CreateNewNode(4, 0.0, 1.0, 0.0);
        NodeType::Pointer p_node_5 = test_model_part.CreateNewNode(5, 0.0, 0.0, 1.0);
        NodeType::Pointer p_node_6 = test_model_part.CreateNewNode(6, 1.0, 0.0, 1.0);
        NodeType::Pointer p_node_7 = test_model_part.CreateNewNode(7, 1.0, 1.0, 1.0);
        NodeType::Pointer p_node_8 = test_model_part.CreateNewNode(8, 0.0, 1.0, 1.0);
        Hexahedra3D8<NodeType> geometry = Hexahedra3D8<NodeType>(p_node_1, p_node_2, p_node_3, p_node_4,
            p_node_5, p_node_6, p_node_7, p_node_8);
        const double mp_volume = 0.000027; //30mm char length
        geometry.SetValue(MP_VOLUME, mp_volume);

        // Set material properties
        const double density_ref = 2314.0;
        material_properties.SetValue(DENSITY, density_ref);
        material_properties.SetValue(SHEAR_MODULUS, 16700e6);
        material_properties.SetValue(RHT_A, 1.6);
        material_properties.SetValue(RHT_N, 0.61);
        material_properties.SetValue(RHT_COMPRESSIVE_STRENGTH, 35e6);
        material_properties.SetValue(RHT_RELATIVE_SHEAR_STRENGTH, 0.18);
        material_properties.SetValue(RHT_RELATIVE_TENSILE_STRENGTH, 0.1);
        material_properties.SetValue(RHT_Q0, 0.6805);
        material_properties.SetValue(RHT_B, 0.0105);
        material_properties.SetValue(REFERENCE_TENSION_STRAIN_RATE, 3e-6);
        material_properties.SetValue(REFERENCE_COMPRESSION_STRAIN_RATE, 3e-5);
        material_properties.SetValue(RHT_GC_STAR, 0.53);
        material_properties.SetValue(RHT_GT_STAR, 0.7);
        material_properties.SetValue(RHT_SHEAR_MOD_REDUCTION_FACTOR, 0.5);
        material_properties.SetValue(RHT_D1, 0.04);
        material_properties.SetValue(RHT_D2, 1.0);
        material_properties.SetValue(RHT_MIN_DAMAGED_RESIDUAL_STRAIN, 0.01);
        material_properties.SetValue(RHT_AF, 1.6);
        material_properties.SetValue(RHT_NF, 0.61);
        material_properties.SetValue(RHT_EOS_A1, 3.527e10);
        material_properties.SetValue(RHT_EOS_A2, 3.958e10);
        material_properties.SetValue(RHT_EOS_A3, 9.04e9);
        material_properties.SetValue(RHT_EOS_B0, 1.22);
        material_properties.SetValue(RHT_EOS_B1, 1.22);
        material_properties.SetValue(RHT_EOS_T1, 3.527e10);
        material_properties.SetValue(RHT_EOS_T2, 0.0);
        material_properties.SetValue(RHT_EOS_ALPHA0, 1.1884);
        material_properties.SetValue(RHT_EOS_NP, 3.0);
        material_properties.SetValue(RHT_CRUSH_PRESSURE, 33e6);
        material_properties.SetValue(RHT_COMPACTION_PRESSURE, 6000e6);
        material_properties.SetValue(FRACTURE_ENERGY, 120.0);


        test_model_part.GetProcessInfo().SetValue(DELTA_TIME, 1e12); // disable strain rate effect
        test_model_part.GetProcessInfo().SetValue(IS_EXPLICIT, true);

        // Set constitutive law flags:
        Flags& ConstitutiveLawOptions = cl_parameters.GetOptions();
        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

        // Set required constitutive law parameters:
        cl_parameters.SetElementGeometry(geometry);
        cl_parameters.SetProcessInfo(test_model_part.GetProcessInfo());
        cl_parameters.SetMaterialProperties(material_properties);
        cl_parameters.SetStrainVector(strain_vector);
        cl_parameters.SetStressVector(stress_vector);

        // Create the CL
        RHTConcrete3DLaw cl = RHTConcrete3DLaw();


        // Set strain for the test
        const double strain_min = -20000.0;
        const double strain_increments = 500.0;
        const int strain_steps = 41;

        Matrix result_matrix(strain_steps,7);
        Matrix strain_matrix;

        // results: p el f res yield converged eps_star damage
        if (print_results) std::cout << "\n\np\tvm*\teps_hard_ratio\teps\teps_rate\tdamage\talpha\n";

        for (size_t i = 0; i < strain_steps; ++i)
        {
            // Initialize the CL
            Vector dummy;
            material_properties.SetValue(DENSITY, density_ref);
            cl.InitializeMaterial(material_properties, geometry, dummy);

            // Set the parameters
            stress_vector.clear();
            strain_vector.clear();
            strain_vector[0] = (strain_min+ i*strain_increments) / 1e6;
            strain_matrix = ZeroMatrix(3);
            strain_matrix(0, 0) = strain_vector[0]; // exx

            //KRATOS_WATCH(strain_matrix)

            Matrix F = IdentityMatrix(3) + strain_matrix;
            const double detF = MathUtils<double>::Det(F);
            material_properties.SetValue(DENSITY, density_ref/ detF);

            cl_parameters.SetDeformationGradientF(F);
            cl_parameters.SetDeterminantF(detF);
            Matrix dummy_matrix;
            Vector dummy_vector;
            cl_parameters.SetConstitutiveMatrix(dummy_matrix);
            cl_parameters.SetShapeFunctionsValues(dummy_vector);
            cl_parameters.SetShapeFunctionsDerivatives(dummy_matrix);

            // Run test
            double value_damage, value_eps, value_eps_rate, value_vmstress,
                value_pressure, value_alpha, value_hard_ratio;
            cl.CalculateMaterialResponseKirchhoff(cl_parameters);

            // Check MP_DAMAGE
            cl.GetValue(MP_DAMAGE, value_damage);
            //KRATOS_CHECK_NEAR(ref_temperature, value, tolerance);

            // Check MP_EQUIVALENT_PLASTIC_STRAIN
            cl.GetValue(MP_HARDENING_RATIO, value_hard_ratio);

            // Check MP_EQUIVALENT_PLASTIC_STRAIN
            cl.GetValue(MP_EQUIVALENT_PLASTIC_STRAIN, value_eps);
            //KRATOS_CHECK_NEAR(ref_plastic_strain, value, tolerance * tolerance);

            // Check MP_EQUIVALENT_PLASTIC_STRAIN_RATE
            cl.GetValue(MP_EQUIVALENT_PLASTIC_STRAIN_RATE, value_eps_rate);
            //KRATOS_CHECK_NEAR(ref_plastic_strain_rate, value, tolerance * tolerance);

            // Check MP_EQUIVALENT_STRESS
            cl.GetValue(MP_EQUIVALENT_STRESS, value_vmstress);
            //KRATOS_CHECK_NEAR(ref_equivalent_stress, value, tolerance);

            // Check MP_PRESSURE
            cl.GetValue(MP_PRESSURE, value_pressure);
            //KRATOS_CHECK_NEAR(ref_equivalent_stress, value, tolerance);

            // Check MP_COMPACTION_RATIO
            cl.GetValue(MP_COMPACTION_RATIO, value_alpha);
            //KRATOS_CHECK_NEAR(ref_equivalent_stress, value, tolerance);

            KRATOS_WATCH(material_properties[DENSITY])

            result_matrix(i, 0) = value_pressure / 1e6;
            result_matrix(i, 1) = value_vmstress / material_properties[RHT_COMPRESSIVE_STRENGTH];
            result_matrix(i, 2) = value_hard_ratio;
            result_matrix(i, 3) = value_eps;
            result_matrix(i, 4) = value_eps_rate;
            result_matrix(i, 5) = value_damage;
            result_matrix(i, 6) = value_alpha;

            if (print_results) {
                for (size_t j = 0; j < result_matrix.size2(); ++j) {
                    if (j > 0) std::cout << ", ";
                    std::cout << result_matrix(i, j);
                }
                std::cout << "\n";
            }

            int test = 1;
        }



    }
} // namespace Testing
} // namespace Kratos
