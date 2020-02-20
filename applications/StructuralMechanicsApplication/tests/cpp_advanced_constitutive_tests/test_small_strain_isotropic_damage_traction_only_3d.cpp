// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
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

// Constitutive law
#include "custom_advanced_constitutive/small_strain_isotropic_damage_traction_only_3d.h"
#include "includes/model_part.h"
#include "geometries/tetrahedra_3d_4.h"

namespace Kratos
{
namespace Testing
{
// We test the associated plasticity Constitutive laws...
typedef Node<3> NodeType;

// Check the correct calculation of the integrated stress with the CL's in small strain
KRATOS_TEST_CASE_IN_SUITE(_ConstitutiveLaw_SmallStrainIsotropicDamageTractionOnly3D, KratosStructuralMechanicsFastSuite)
{
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
    material_properties.SetValue(YOUNG_MODULUS, 3000);
    material_properties.SetValue(POISSON_RATIO, 0.3);
    material_properties.SetValue(YIELD_STRESS, 0.5);
    material_properties.SetValue(INFINITY_YIELD_STRESS, 0.7);
    Vector hardening_moduli(2);
    hardening_moduli(0) = 0.3; hardening_moduli(1) = 0.15;
    material_properties.SetValue(HARDENING_MODULI_VECTOR, hardening_moduli);
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
    SmallStrainIsotropicDamageTractionOnly3D cl = SmallStrainIsotropicDamageTractionOnly3D();

    // Set variables for the test
    const double tolerance = 1.0e-4;
    std::size_t nr_ts = 4;  // timesteps tested
    Vector ser(nr_ts), dvr(nr_ts);  // reference strain energy, damage variable
    Matrix epr(nr_ts, 6), str(nr_ts, 6);  // reference strain ("epsilon"), stress
    Matrix cmr(nr_ts, 6 * 6);  // reference constitutive matrix


    //
    // Test: check correct behavior of internal and calculated variables
    //
    KRATOS_CHECK_IS_FALSE(cl.Has(STRAIN_ENERGY));  // = False, in order to use CalculateValue())
    KRATOS_CHECK_IS_FALSE(cl.Has(DAMAGE_VARIABLE));  // = False, in order to use CalculateValue())

    //
    // Test: load - unload in traction
    //
    epr(0,0)=1.000000e-04; epr(0,1)=1.000000e-04; epr(0,2)=0.000000e+00; epr(0,3)=1.000000e-04; epr(0,4)=0.000000e+00; epr(0,5)=0.000000e+00;
    epr(1,0)=3.000000e-04; epr(1,1)=3.000000e-04; epr(1,2)=0.000000e+00; epr(1,3)=3.000000e-04; epr(1,4)=0.000000e+00; epr(1,5)=0.000000e+00;
    epr(2,0)=1.000000e-03; epr(2,1)=1.000000e-03; epr(2,2)=0.000000e+00; epr(2,3)=1.000000e-03; epr(2,4)=0.000000e+00; epr(2,5)=0.000000e+00;
    epr(3,0)=-1.000000e-03; epr(3,1)=-1.000000e-03; epr(3,2)=0.000000e+00; epr(3,3)=-1.000000e-03; epr(3,4)=0.000000e+00; epr(3,5)=0.000000e+00;

    str(0,0)=5.003084e-01; str(0,1)=5.003084e-01; str(0,2)=3.001850e-01; str(0,3)=1.000617e-01; str(0,4)=0.000000e+00; str(0,5)=0.000000e+00;
    str(1,0)=7.504626e-01; str(1,1)=7.504626e-01; str(1,2)=4.502775e-01; str(1,3)=1.500925e-01; str(1,4)=0.000000e+00; str(1,5)=0.000000e+00;
    str(2,0)=1.356232e+00; str(2,1)=1.356232e+00; str(2,2)=8.137391e-01; str(2,3)=2.712464e-01; str(2,4)=0.000000e+00; str(2,5)=0.000000e+00;
    str(3,0)=-1.356232e+00; str(3,1)=-1.356232e+00; str(3,2)=-8.137391e-01; str(3,3)=-2.712464e-01; str(3,4)=0.000000e+00; str(3,5)=0.000000e+00;

    cmr(0, 0)=2.014743e+03; cmr(0, 1)=1.350944e+01; cmr(0, 2)=6.084757e+02; cmr(0, 3)=-2.974831e+02; cmr(0, 4)=0.000000e+00; cmr(0, 5)=0.000000e+00;
    cmr(0, 6)=1.350944e+01; cmr(0, 7)=2.014743e+03; cmr(0, 8)=6.084757e+02; cmr(0, 9)=-2.974831e+02; cmr(0,10)=0.000000e+00; cmr(0,11)=0.000000e+00;
    cmr(0,12)=6.084757e+02; cmr(0,13)=6.084757e+02; cmr(0,14)=2.966689e+03; cmr(0,15)=-1.784899e+02; cmr(0,16)=0.000000e+00; cmr(0,17)=0.000000e+00;
    cmr(0,18)=-2.974831e+02; cmr(0,19)=-2.974831e+02; cmr(0,20)=-1.784899e+02; cmr(0,21)=9.411201e+02; cmr(0,22)=0.000000e+00; cmr(0,23)=0.000000e+00;
    cmr(0,24)=0.000000e+00; cmr(0,25)=0.000000e+00; cmr(0,26)=0.000000e+00; cmr(0,27)=0.000000e+00; cmr(0,28)=1.000617e+03; cmr(0,29)=0.000000e+00;
    cmr(0,30)=0.000000e+00; cmr(0,31)=0.000000e+00; cmr(0,32)=0.000000e+00; cmr(0,33)=0.000000e+00; cmr(0,34)=0.000000e+00; cmr(0,35)=1.000617e+03;
    cmr(1, 0)=1.007371e+03; cmr(1, 1)=6.754721e+00; cmr(1, 2)=3.042379e+02; cmr(1, 3)=-1.487416e+02; cmr(1, 4)=0.000000e+00; cmr(1, 5)=0.000000e+00;
    cmr(1, 6)=6.754721e+00; cmr(1, 7)=1.007371e+03; cmr(1, 8)=3.042379e+02; cmr(1, 9)=-1.487416e+02; cmr(1,10)=0.000000e+00; cmr(1,11)=0.000000e+00;
    cmr(1,12)=3.042379e+02; cmr(1,13)=3.042379e+02; cmr(1,14)=1.483344e+03; cmr(1,15)=-8.924494e+01; cmr(1,16)=0.000000e+00; cmr(1,17)=0.000000e+00;
    cmr(1,18)=-1.487416e+02; cmr(1,19)=-1.487416e+02; cmr(1,20)=-8.924494e+01; cmr(1,21)=4.705601e+02; cmr(1,22)=0.000000e+00; cmr(1,23)=0.000000e+00;
    cmr(1,24)=0.000000e+00; cmr(1,25)=0.000000e+00; cmr(1,26)=0.000000e+00; cmr(1,27)=0.000000e+00; cmr(1,28)=5.003084e+02; cmr(1,29)=0.000000e+00;
    cmr(1,30)=0.000000e+00; cmr(1,31)=0.000000e+00; cmr(1,32)=0.000000e+00; cmr(1,33)=0.000000e+00; cmr(1,34)=0.000000e+00; cmr(1,35)=5.003084e+02;
    cmr(2, 0)=7.262499e+02; cmr(2, 1)=1.837572e+02; cmr(2, 2)=2.730021e+02; cmr(2, 3)=-4.462247e+01; cmr(2, 4)=0.000000e+00; cmr(2, 5)=0.000000e+00;
    cmr(2, 6)=1.837572e+02; cmr(2, 7)=7.262499e+02; cmr(2, 8)=2.730021e+02; cmr(2, 9)=-4.462247e+01; cmr(2,10)=0.000000e+00; cmr(2,11)=0.000000e+00;
    cmr(2,12)=2.730021e+02; cmr(2,13)=2.730021e+02; cmr(2,14)=8.690418e+02; cmr(2,15)=-2.677348e+01; cmr(2,16)=0.000000e+00; cmr(2,17)=0.000000e+00;
    cmr(2,18)=-4.462247e+01; cmr(2,19)=-4.462247e+01; cmr(2,20)=-2.677348e+01; cmr(2,21)=2.623219e+02; cmr(2,22)=0.000000e+00; cmr(2,23)=0.000000e+00;
    cmr(2,24)=0.000000e+00; cmr(2,25)=0.000000e+00; cmr(2,26)=0.000000e+00; cmr(2,27)=0.000000e+00; cmr(2,28)=2.712464e+02; cmr(2,29)=0.000000e+00;
    cmr(2,30)=0.000000e+00; cmr(2,31)=0.000000e+00; cmr(2,32)=0.000000e+00; cmr(2,33)=0.000000e+00; cmr(2,34)=0.000000e+00; cmr(2,35)=2.712464e+02;
    cmr(3, 0)=9.493622e+02; cmr(3, 1)=4.068695e+02; cmr(3, 2)=4.068695e+02; cmr(3, 3)=0.000000e+00; cmr(3, 4)=0.000000e+00; cmr(3, 5)=0.000000e+00;
    cmr(3, 6)=4.068695e+02; cmr(3, 7)=9.493622e+02; cmr(3, 8)=4.068695e+02; cmr(3, 9)=0.000000e+00; cmr(3,10)=0.000000e+00; cmr(3,11)=0.000000e+00;
    cmr(3,12)=4.068695e+02; cmr(3,13)=4.068695e+02; cmr(3,14)=9.493622e+02; cmr(3,15)=0.000000e+00; cmr(3,16)=0.000000e+00; cmr(3,17)=0.000000e+00;
    cmr(3,18)=0.000000e+00; cmr(3,19)=0.000000e+00; cmr(3,20)=0.000000e+00; cmr(3,21)=2.712464e+02; cmr(3,22)=0.000000e+00; cmr(3,23)=0.000000e+00;
    cmr(3,24)=0.000000e+00; cmr(3,25)=0.000000e+00; cmr(3,26)=0.000000e+00; cmr(3,27)=0.000000e+00; cmr(3,28)=2.712464e+02; cmr(3,29)=0.000000e+00;
    cmr(3,30)=0.000000e+00; cmr(3,31)=0.000000e+00; cmr(3,32)=0.000000e+00; cmr(3,33)=0.000000e+00; cmr(3,34)=0.000000e+00; cmr(3,35)=2.712464e+02;

    ser[0]=5.503392e-05; ser[1]=2.476526e-04; ser[2]=1.491855e-03; ser[3]=1.491855e-03;
    dvr[0]=1.327988e-01; dvr[1]=5.663994e-01; dvr[2]=7.649198e-01; dvr[3]=7.649198e-01;

    // Here we must simulate the call sequence of the element
    Vector dummy;
    cl.InitializeMaterial(material_properties, geometry, dummy);
    for (std::size_t t = 0; t < nr_ts; ++t){
        for (std::size_t comp = 0; comp < 6; ++comp) {
            strain_vector[comp] = epr(t, comp);
        }
        if (cl.RequiresInitializeMaterialResponse()){
            cl.InitializeMaterialResponseCauchy(cl_parameters);
        }
        cl.CalculateMaterialResponseCauchy(cl_parameters);
        if (cl.RequiresFinalizeMaterialResponse()) {
            cl.FinalizeMaterialResponseCauchy(cl_parameters);
        }
        double value;

        // Check damage variable
        cl.CalculateValue(cl_parameters, DAMAGE_VARIABLE, value);
        // TODO(marandra): NAN values are not handled correctly by KRATOS CHECK functions
        if (std::isnan(dvr[t]/value)){
            KRATOS_CHECK_NEAR(dvr[t], value, tolerance);
        } else {
            KRATOS_CHECK_NEAR(dvr[t]/value, 1, tolerance);
        }

        // Check strain energy
        cl.CalculateValue(cl_parameters, STRAIN_ENERGY, value);
        if (std::isnan(ser[t]/value)){
            KRATOS_CHECK_NEAR(ser[t], value, tolerance);
        } else {
            KRATOS_CHECK_NEAR(ser[t]/value, 1, tolerance);
        }

        // Check stress
        for (std::size_t comp = 0; comp < 6; ++comp){
            KRATOS_CHECK_IS_FALSE(std::isnan(stress_vector[comp]));
            if (std::isnan(stress_vector[comp]/str(t, comp))){
                KRATOS_CHECK_NEAR(stress_vector[comp], str(t, comp), tolerance);
            } else {
                KRATOS_CHECK_NEAR(stress_vector[comp]/str(t, comp), 1, tolerance);
            }
        }

        // Check constitutive tensor
        for (std::size_t i = 0; i < 6; ++i){
            for (std::size_t j = 0; j < 6; ++j){
                std::size_t idx = i * 6 + j;
                KRATOS_CHECK_IS_FALSE(std::isnan(const_matrix(i, j)));
                if (std::isnan(const_matrix(i, j)/cmr(t, idx))){
                    KRATOS_CHECK_NEAR(const_matrix(i, j), cmr(t, idx), tolerance);
                } else {
                    KRATOS_CHECK_NEAR(const_matrix(i, j)/cmr(t, idx), 1, tolerance);
                    }
            }
        }
    }


    //
    // Test: load - unload in compression
    //
    epr(0,0)=-1.000000e-04; epr(0,1)=-1.000000e-04; epr(0,2)=0.000000e+00; epr(0,3)=-1.000000e-04; epr(0,4)=0.000000e+00; epr(0,5)=0.000000e+00;
    epr(1,0)=-3.000000e-04; epr(1,1)=-3.000000e-04; epr(1,2)=0.000000e+00; epr(1,3)=-3.000000e-04; epr(1,4)=0.000000e+00; epr(1,5)=0.000000e+00;
    epr(2,0)=-1.000000e-03; epr(2,1)=-1.000000e-03; epr(2,2)=0.000000e+00; epr(2,3)=-1.000000e-03; epr(2,4)=0.000000e+00; epr(2,5)=0.000000e+00;
    epr(3,0)=1.000000e-03; epr(3,1)=1.000000e-03; epr(3,2)=0.000000e+00; epr(3,3)=1.000000e-03; epr(3,4)=0.000000e+00; epr(3,5)=0.000000e+00;

    str(0,0)=-5.769231e-01; str(0,1)=-5.769231e-01; str(0,2)=-3.461538e-01; str(0,3)=-1.153846e-01; str(0,4)=0.000000e+00; str(0,5)=0.000000e+00;
    str(1,0)=-1.730769e+00; str(1,1)=-1.730769e+00; str(1,2)=-1.038462e+00; str(1,3)=-3.461538e-01; str(1,4)=0.000000e+00; str(1,5)=0.000000e+00;
    str(2,0)=-5.769231e+00; str(2,1)=-5.769231e+00; str(2,2)=-3.461538e+00; str(2,3)=-1.153846e+00; str(2,4)=0.000000e+00; str(2,5)=0.000000e+00;
    str(3,0)=1.356232e+00; str(3,1)=1.356232e+00; str(3,2)=8.137391e-01; str(3,3)=2.712464e-01; str(3,4)=0.000000e+00; str(3,5)=0.000000e+00;

    cmr(0, 0)=4.038462e+03; cmr(0, 1)=1.730769e+03; cmr(0, 2)=1.730769e+03; cmr(0, 3)=0.000000e+00; cmr(0, 4)=0.000000e+00; cmr(0, 5)=0.000000e+00;
    cmr(0, 6)=1.730769e+03; cmr(0, 7)=4.038462e+03; cmr(0, 8)=1.730769e+03; cmr(0, 9)=0.000000e+00; cmr(0,10)=0.000000e+00; cmr(0,11)=0.000000e+00;
    cmr(0,12)=1.730769e+03; cmr(0,13)=1.730769e+03; cmr(0,14)=4.038462e+03; cmr(0,15)=0.000000e+00; cmr(0,16)=0.000000e+00; cmr(0,17)=0.000000e+00;
    cmr(0,18)=0.000000e+00; cmr(0,19)=0.000000e+00; cmr(0,20)=0.000000e+00; cmr(0,21)=1.153846e+03; cmr(0,22)=0.000000e+00; cmr(0,23)=0.000000e+00;
    cmr(0,24)=0.000000e+00; cmr(0,25)=0.000000e+00; cmr(0,26)=0.000000e+00; cmr(0,27)=0.000000e+00; cmr(0,28)=1.153846e+03; cmr(0,29)=0.000000e+00;
    cmr(0,30)=0.000000e+00; cmr(0,31)=0.000000e+00; cmr(0,32)=0.000000e+00; cmr(0,33)=0.000000e+00; cmr(0,34)=0.000000e+00; cmr(0,35)=1.153846e+03;
    cmr(1, 0)=4.038462e+03; cmr(1, 1)=1.730769e+03; cmr(1, 2)=1.730769e+03; cmr(1, 3)=0.000000e+00; cmr(1, 4)=0.000000e+00; cmr(1, 5)=0.000000e+00;
    cmr(1, 6)=1.730769e+03; cmr(1, 7)=4.038462e+03; cmr(1, 8)=1.730769e+03; cmr(1, 9)=0.000000e+00; cmr(1,10)=0.000000e+00; cmr(1,11)=0.000000e+00;
    cmr(1,12)=1.730769e+03; cmr(1,13)=1.730769e+03; cmr(1,14)=4.038462e+03; cmr(1,15)=0.000000e+00; cmr(1,16)=0.000000e+00; cmr(1,17)=0.000000e+00;
    cmr(1,18)=0.000000e+00; cmr(1,19)=0.000000e+00; cmr(1,20)=0.000000e+00; cmr(1,21)=1.153846e+03; cmr(1,22)=0.000000e+00; cmr(1,23)=0.000000e+00;
    cmr(1,24)=0.000000e+00; cmr(1,25)=0.000000e+00; cmr(1,26)=0.000000e+00; cmr(1,27)=0.000000e+00; cmr(1,28)=1.153846e+03; cmr(1,29)=0.000000e+00;
    cmr(1,30)=0.000000e+00; cmr(1,31)=0.000000e+00; cmr(1,32)=0.000000e+00; cmr(1,33)=0.000000e+00; cmr(1,34)=0.000000e+00; cmr(1,35)=1.153846e+03;
    cmr(2, 0)=4.038462e+03; cmr(2, 1)=1.730769e+03; cmr(2, 2)=1.730769e+03; cmr(2, 3)=0.000000e+00; cmr(2, 4)=0.000000e+00; cmr(2, 5)=0.000000e+00;
    cmr(2, 6)=1.730769e+03; cmr(2, 7)=4.038462e+03; cmr(2, 8)=1.730769e+03; cmr(2, 9)=0.000000e+00; cmr(2,10)=0.000000e+00; cmr(2,11)=0.000000e+00;
    cmr(2,12)=1.730769e+03; cmr(2,13)=1.730769e+03; cmr(2,14)=4.038462e+03; cmr(2,15)=0.000000e+00; cmr(2,16)=0.000000e+00; cmr(2,17)=0.000000e+00;
    cmr(2,18)=0.000000e+00; cmr(2,19)=0.000000e+00; cmr(2,20)=0.000000e+00; cmr(2,21)=1.153846e+03; cmr(2,22)=0.000000e+00; cmr(2,23)=0.000000e+00;
    cmr(2,24)=0.000000e+00; cmr(2,25)=0.000000e+00; cmr(2,26)=0.000000e+00; cmr(2,27)=0.000000e+00; cmr(2,28)=1.153846e+03; cmr(2,29)=0.000000e+00;
    cmr(2,30)=0.000000e+00; cmr(2,31)=0.000000e+00; cmr(2,32)=0.000000e+00; cmr(2,33)=0.000000e+00; cmr(2,34)=0.000000e+00; cmr(2,35)=1.153846e+03;
    cmr(3, 0)=7.262499e+02; cmr(3, 1)=1.837572e+02; cmr(3, 2)=2.730021e+02; cmr(3, 3)=-4.462247e+01; cmr(3, 4)=0.000000e+00; cmr(3, 5)=0.000000e+00;
    cmr(3, 6)=1.837572e+02; cmr(3, 7)=7.262499e+02; cmr(3, 8)=2.730021e+02; cmr(3, 9)=-4.462247e+01; cmr(3,10)=0.000000e+00; cmr(3,11)=0.000000e+00;
    cmr(3,12)=2.730021e+02; cmr(3,13)=2.730021e+02; cmr(3,14)=8.690418e+02; cmr(3,15)=-2.677348e+01; cmr(3,16)=0.000000e+00; cmr(3,17)=0.000000e+00;
    cmr(3,18)=-4.462247e+01; cmr(3,19)=-4.462247e+01; cmr(3,20)=-2.677348e+01; cmr(3,21)=2.623219e+02; cmr(3,22)=0.000000e+00; cmr(3,23)=0.000000e+00;
    cmr(3,24)=0.000000e+00; cmr(3,25)=0.000000e+00; cmr(3,26)=0.000000e+00; cmr(3,27)=0.000000e+00; cmr(3,28)=2.712464e+02; cmr(3,29)=0.000000e+00;
    cmr(3,30)=0.000000e+00; cmr(3,31)=0.000000e+00; cmr(3,32)=0.000000e+00; cmr(3,33)=0.000000e+00; cmr(3,34)=0.000000e+00; cmr(3,35)=2.712464e+02;

    ser[0]=6.346154e-05; ser[1]=5.711538e-04; ser[2]=6.346154e-03; ser[3]=1.491855e-03;
    dvr[0]=0.000000e+00; dvr[1]=0.000000e+00; dvr[2]=0.000000e+00; dvr[3]=7.649198e-01;

    // Here we must simulate the call sequence of the element
    cl.InitializeMaterial(material_properties, geometry, dummy);
    for (std::size_t t = 0; t < nr_ts; ++t){
        for (std::size_t comp = 0; comp < 6; ++comp) {
            strain_vector[comp] = epr(t, comp);
        }
        cl.InitializeMaterialResponseCauchy(cl_parameters);
        cl.CalculateMaterialResponseCauchy(cl_parameters);
        cl.FinalizeMaterialResponseCauchy(cl_parameters);
        double value;

        // Check damage variable
        cl.CalculateValue(cl_parameters, DAMAGE_VARIABLE, value);
        if (std::isnan(dvr[t]/value)){
            KRATOS_CHECK_NEAR(dvr[t], value, tolerance);
        } else {
            KRATOS_CHECK_NEAR(dvr[t]/value, 1, tolerance);
        }

        // Check strain energy
        cl.CalculateValue(cl_parameters, STRAIN_ENERGY, value);
        if (std::isnan(ser[t]/value)){
            KRATOS_CHECK_NEAR(ser[t], value, tolerance);
        } else {
            KRATOS_CHECK_NEAR(ser[t]/value, 1, tolerance);
        }

        // Check stress
        for (std::size_t comp = 0; comp < 6; ++comp){
            KRATOS_CHECK_IS_FALSE(std::isnan(stress_vector[comp]));
            if (std::isnan(stress_vector[comp]/str(t, comp))){
                KRATOS_CHECK_NEAR(stress_vector[comp], str(t, comp), tolerance);
            } else {
                KRATOS_CHECK_NEAR(stress_vector[comp]/str(t, comp), 1, tolerance);
            }
        }

        // Check constitutive tensor
        for (std::size_t i = 0; i < 6; ++i){
            for (std::size_t j = 0; j < 6; ++j){
                std::size_t idx = i * 6 + j;
                KRATOS_CHECK_IS_FALSE(std::isnan(const_matrix(i, j)));
                if (std::isnan(const_matrix(i, j)/cmr(t, idx))){
                    KRATOS_CHECK_NEAR(const_matrix(i, j), cmr(t, idx), tolerance);
                } else {
                    KRATOS_CHECK_NEAR(const_matrix(i, j)/cmr(t, idx), 1, tolerance);
                }
            }
        }
    }


}

} // namespace Testing
} // namespace Kratos
