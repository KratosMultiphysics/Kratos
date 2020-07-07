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

    KRATOS_TEST_CASE_IN_SUITE(ParticleConstitutiveLawRHTConcrete3DUniaxialCompression, KratosParticleMechanicsFastSuite)
    {
        const bool print_results = false;

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
        std::vector<std::vector<double>> ref_results(strain_steps);

        if (print_results) std::cout <<
            "\n\npressure\tvm/pressure\teps_hard_ratio\teps\teps_rate\tdamage\talpha\n";

        // Run the CL
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

            // Check MP_EQUIVALENT_PLASTIC_STRAIN
            cl.GetValue(MP_HARDENING_RATIO, value_hard_ratio);

            // Check MP_EQUIVALENT_PLASTIC_STRAIN
            cl.GetValue(MP_EQUIVALENT_PLASTIC_STRAIN, value_eps);

            // Check MP_EQUIVALENT_PLASTIC_STRAIN_RATE
            cl.GetValue(MP_EQUIVALENT_PLASTIC_STRAIN_RATE, value_eps_rate);

            // Check MP_EQUIVALENT_STRESS
            cl.GetValue(MP_EQUIVALENT_STRESS, value_vmstress);

            // Check MP_PRESSURE
            cl.GetValue(MP_PRESSURE, value_pressure);

            // Check MP_COMPACTION_RATIO
            cl.GetValue(MP_COMPACTION_RATIO, value_alpha);

            result_matrix(i, 0) = value_pressure / 1e6;
            result_matrix(i, 1) = value_vmstress / material_properties[RHT_COMPRESSIVE_STRENGTH];
            result_matrix(i, 2) = value_hard_ratio;
            result_matrix(i, 3) = value_eps;
            result_matrix(i, 4) = value_eps_rate;
            result_matrix(i, 5) = value_damage;
            result_matrix(i, 6) = value_alpha;

            if (print_results) {
                std::cout << "ref_results[" << i << "] = {";
                for (size_t j = 0; j < result_matrix.size2(); ++j) {
                    if (j > 0) std::cout << ", ";
                    std::cout << std::fixed << std::setprecision(9) << result_matrix(i, j);
                }
                std::cout << "};\n";
            }
        }

        // Assemble reference results
        // Result:         Pressure     , VM/pressure, eps/eps_hard, eps       , eps_rate   ,  damage    , alpha
        ref_results[0] = { 204.615514977, 4.697978205, 1.000000000, 0.010051313, 0.000000000, 1.000000000, 1.172607454 };
        ref_results[1] = { 199.946090557, 4.632285894, 1.000000000, 0.009763872, 0.000000000, 1.000000000, 1.173025007 };
        ref_results[2] = { 195.283860293, 4.566095072, 1.000000000, 0.009476780, 0.000000000, 1.000000000, 1.173442588 };
        ref_results[3] = { 190.628798380, 4.499388219, 1.000000000, 0.009190048, 0.000000000, 1.000000000, 1.173860197 };
        ref_results[4] = { 185.980879151, 4.432146799, 1.000000000, 0.008903690, 0.000000000, 1.000000000, 1.174277833 };
        ref_results[5] = { 181.340077068, 4.364351200, 1.000000000, 0.008617719, 0.000000000, 1.000000000, 1.174695497 };
        ref_results[6] = { 176.706366726, 4.297454222, 1.000000000, 0.008331120, 0.000000000, 0.977959760, 1.175113188 };
        ref_results[7] = { 172.079722849, 4.234227581, 1.000000000, 0.008041957, 0.000000000, 0.893189623, 1.175530905 };
        ref_results[8] = { 167.460120294, 4.170465575, 1.000000000, 0.007753168, 0.000000000, 0.808947286, 1.175948649 };
        ref_results[9] = { 162.847534042, 4.106147385, 1.000000000, 0.007464767, 0.000000000, 0.725253395, 1.176366419 };
        ref_results[10] = { 158.241939205, 4.041250868, 1.000000000, 0.007176771, 0.000000000, 0.642129923, 1.176784215 };
        ref_results[11] = { 153.643311020, 3.975752462, 1.000000000, 0.006889195, 0.000000000, 0.559600287, 1.177202037 };
        ref_results[12] = { 149.051624851, 3.909627044, 1.000000000, 0.006602057, 0.000000000, 0.477689484, 1.177619885 };
        ref_results[13] = { 144.466856186, 3.842847769, 1.000000000, 0.006315376, 0.000000000, 0.396424257, 1.178037758 };
        ref_results[14] = { 139.888980639, 3.775385924, 1.000000000, 0.006029172, 0.000000000, 0.315833259, 1.178455656 };
        ref_results[15] = { 135.317973946, 3.707210717, 1.000000000, 0.005743466, 0.000000000, 0.235947262, 1.178873579 };
        ref_results[16] = { 130.753811966, 3.638289057, 1.000000000, 0.005458281, 0.000000000, 0.156799388, 1.179291526 };
        ref_results[17] = { 126.196470677, 3.568585333, 1.000000000, 0.005173643, 0.000000000, 0.078425359, 1.179709499 };
        ref_results[18] = { 121.645926182, 3.498061099, 1.000000000, 0.004889578, 0.000000000, 0.000863809, 1.180127495 };
        ref_results[19] = { 117.102154701, 3.339999999, 0.976398008, 0.004666667, 0.000000000, 0.014449873, 1.180545516 };
        ref_results[20] = { 112.565132574, 3.180952380, 0.951718157, 0.004444444, 0.000000000, 0.144411475, 1.180963560 };
        ref_results[21] = { 108.034836258, 3.021904761, 0.926157641, 0.004222222, 0.000000000, 0.083624128, 1.181381628 };
        ref_results[22] = { 103.511242331, 2.862857142, 0.899642126, 0.004000000, 0.000000000, 0.023532063, 1.181799720 };
        ref_results[23] = { 98.994327482, 2.703809523, 0.872087144, 0.003777778, 0.000000000, 0.000000000, 1.182217835 };
        ref_results[24] = { 94.484068522, 2.544761904, 0.843396167, 0.003555556, 0.000000000, 0.000000000, 1.182635972 };
        ref_results[25] = { 89.980442371, 2.385714285, 0.813458205, 0.003333333, 0.000000000, 0.000000000, 1.183054133 };
        ref_results[26] = { 85.483426068, 2.226666666, 0.782144771, 0.003111111, 0.000000000, 0.000000000, 1.183472316 };
        ref_results[27] = { 80.992996762, 2.067619047, 0.749306023, 0.002888889, 0.000000000, 0.000000000, 1.183890522 };
        ref_results[28] = { 76.509131718, 1.908571429, 0.714765781, 0.002666667, 0.000000000, 0.000000000, 1.184308750 };
        ref_results[29] = { 72.031808310, 1.749523810, 0.678315024, 0.002444444, 0.000000000, 0.000000000, 1.184727001 };
        ref_results[30] = { 67.561004024, 1.590476191, 0.639703272, 0.002222222, 0.000000000, 0.000000000, 1.185145273 };
        ref_results[31] = { 63.096696458, 1.431428572, 0.598626964, 0.002000000, 0.000000000, 0.000000000, 1.185563567 };
        ref_results[32] = { 58.638863317, 1.272380953, 0.554713543, 0.001777778, 0.000000000, 0.000000000, 1.185981882 };
        ref_results[33] = { 54.187482416, 1.113333334, 0.507499161, 0.001555556, 0.000000000, 0.000000000, 1.186400219 };
        ref_results[34] = { 49.742531678, 0.954285713, 0.456396759, 0.001333333, 0.000000000, 0.000000000, 1.186818577 };
        ref_results[35] = { 45.303989134, 0.795238094, 0.400649118, 0.001111111, 0.000000000, 0.000000000, 1.187236956 };
        ref_results[36] = { 40.871832922, 0.636190476, 0.339257677, 0.000888889, 0.000000000, 0.000000000, 1.187655356 };
        ref_results[37] = { 36.446041284, 0.477142858, 0.270870653, 0.000666667, 0.000000000, 0.000000000, 1.188073776 };
        ref_results[38] = { 29.720080934, 0.706706107, 0.124273995, 0.000172960, 0.000000000, 0.000000000, 1.188400000 };
        ref_results[39] = { 14.849653534, 0.477142857, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 1.188400000 };
        ref_results[40] = { 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 0.000000000, 1.188400000 };

        // Check results
        for (size_t i = 0; i < result_matrix.size1(); i++)
        {
            for (size_t j = 0; j < result_matrix.size2(); j++)
            {
                KRATOS_CHECK_NEAR(result_matrix(i,j), ref_results[i][j], 1e-6);
            }
        }
    }
} // namespace Testing
} // namespace Kratos
