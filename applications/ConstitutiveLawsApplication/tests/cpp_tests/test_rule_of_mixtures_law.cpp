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
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "containers/model.h"
#include "utilities/read_materials_utility.h"
#include "includes/mat_variables.h"

// Application includes
#include "constitutive_laws_application_variables.h"
#include "custom_constitutive/composites/rule_of_mixtures_law.h"
#include "tests/cpp_tests/constitutive_laws_fast_suite.h"

namespace Kratos::Testing
{

// void GiDIODebugRuleMixtures(ModelPart& ThisModelPart)
// {
//     GidIO<> gid_io("TEST_RULE_MIXTURES", GiD_PostBinary, SingleFile, WriteUndeformed,  WriteElementsOnly);
//     const int nl_iter = ThisModelPart.GetProcessInfo()[NL_ITERATION_NUMBER];
//     const double label = static_cast<double>(nl_iter);

//     gid_io.InitializeMesh(label);
//     gid_io.WriteMesh(ThisModelPart.GetMesh());
//     gid_io.FinalizeMesh();
//     gid_io.InitializeResults(label, ThisModelPart.GetMesh());
//     gid_io.WriteNodalResults(DISPLACEMENT, ThisModelPart.Nodes(), label, 0);
//     gid_io.PrintOnGaussPoints(GREEN_LAGRANGE_STRAIN_VECTOR, ThisModelPart, label);
//     gid_io.PrintOnGaussPoints(PK2_STRESS_VECTOR, ThisModelPart, label);
//     gid_io.WriteNodalFlags(ACTIVE, "ACTIVE", ThisModelPart.Nodes(), label);
// }

Parameters GetTwoLayersParameters()
{
    Parameters two_layers = Parameters(R"(
    {
    "properties" : [{
        "model_part_name" : "Main",
        "properties_id"   : 1,
        "Material"        : {
            "constitutive_law" : {
                "name" : "ParallelRuleOfMixturesLaw3D",
                "combination_factors"      : [0.4, 0.6 ]
            },
            "Variables"        : {
            },
            "Tables"           : {}
        },
        "sub_properties" : [{
            "properties_id"   : 11,
            "Material"        : {
                "constitutive_law" : {
                    "name" : "LinearElastic3DLaw"
                },
                "Variables"        : {
                    "DENSITY"       : 7850.0,
                    "YOUNG_MODULUS" : 206900000000.0,
                    "POISSON_RATIO" : 0.29
                },
                "Tables"           : {}
            }
        },{
            "properties_id"   : 12,
            "Material"        : {
                "constitutive_law" : {
                    "name" : "HyperElastic3DLaw"
                },
                "Variables"        : {
                    "DENSITY"       : 2000.0,
                    "YOUNG_MODULUS" : 30000000000.0,
                    "POISSON_RATIO" : 0.49
                },
                "Tables"           : {}
            }
        }]
    }]
    })" );

    return two_layers;
}

Parameters GetThreeLayersParameters()
{
    Parameters three_layers = Parameters(R"(
    {
    "properties" : [{
        "model_part_name" : "Main",
        "properties_id"   : 1,
        "Material"        : {
            "constitutive_law" : {
                "name" : "ParallelRuleOfMixturesLaw3D",
                "combination_factors"      : [0.4, 0.3, 0.3 ]
            },
            "Variables"        : {
            },
            "Tables"           : {}
        },
        "sub_properties" : [{
            "properties_id"   : 11,
            "Material"        : {
                "constitutive_law" : {
                    "name" : "LinearElastic3DLaw"
                },
                "Variables"        : {
                    "DENSITY"       : 7850.0,
                    "YOUNG_MODULUS" : 206900000000.0,
                    "POISSON_RATIO" : 0.29
                },
                "Tables"           : {}
            }
        },{
            "properties_id"   : 12,
            "Material"        : {
                "constitutive_law" : {
                    "name" : "LinearElastic3DLaw"
                },
                "Variables"        : {
                    "DENSITY"       : 7850.0,
                    "YOUNG_MODULUS" : 206900000000000.0,
                    "POISSON_RATIO" : 0.29
                },
                "Tables"           : {}
            }
        },{
            "properties_id"   : 13,
            "Material"        : {
                "constitutive_law" : {
                    "name" : "HyperElastic3DLaw"
                },
                "Variables"        : {
                    "DENSITY"       : 2000.0,
                    "YOUNG_MODULUS" : 30000000000000.0,
                    "POISSON_RATIO" : 0.49
                },
                "Tables"           : {}
            }
        }]
    }]
    })" );

    return three_layers;
}

void Create3DGeometryHexahedraRuleOfMixtures(ModelPart& rThisModelPart, std::size_t NumberOfLayers = 2, const std::string& ElementName = "SmallDisplacementElement3D8N")
{
    rThisModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);

    auto& r_process_info = rThisModelPart.GetProcessInfo();
    r_process_info[STEP] = 1;
    r_process_info[NL_ITERATION_NUMBER] = 1;

    Parameters parameters = NumberOfLayers == 2 ? GetTwoLayersParameters() : GetThreeLayersParameters();

    // Read properties
    auto read_util = ReadMaterialsUtility(rThisModelPart.GetModel());
    read_util.ReadMaterials(parameters);

    // Create nodes and elements
    Properties::Pointer p_elem_prop = rThisModelPart.pGetProperties(1);

    // First we create the nodes
    rThisModelPart.CreateNewNode(1 , 0.0 , 1.0 , 1.0);
    rThisModelPart.CreateNewNode(2 , 0.0 , 1.0 , 0.0);
    rThisModelPart.CreateNewNode(3 , 0.0 , 0.0 , 1.0);
    rThisModelPart.CreateNewNode(4 , 1.0 , 1.0 , 1.0);
    rThisModelPart.CreateNewNode(5 , 0.0 , 0.0 , 0.0);
    rThisModelPart.CreateNewNode(6 , 1.0 , 1.0 , 0.0);
    rThisModelPart.CreateNewNode(7 , 1.0 , 0.0 , 1.0);
    rThisModelPart.CreateNewNode(8 , 1.0 , 0.0 , 0.0);

    // Now we create the elements
    rThisModelPart.CreateNewElement(ElementName, 1, {{5,8,6,2,3,7,4,1}}, p_elem_prop);

    const auto& const_process_info = rThisModelPart.GetProcessInfo();

    // Initialize elements
    for (auto& r_elem : rThisModelPart.Elements()) {
        r_elem.Initialize(const_process_info);
        r_elem.InitializeSolutionStep(const_process_info);
        r_elem.InitializeNonLinearIteration(const_process_info);
    }
}

void Create3DGeometryTetrahedraRuleOfMixtures(ModelPart& rThisModelPart, std::size_t NumberOfLayers = 2, const std::string& ElementName = "SmallDisplacementElement3D4N")
{
    rThisModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);

    auto& r_process_info = rThisModelPart.GetProcessInfo();
    r_process_info[STEP] = 1;
    r_process_info[NL_ITERATION_NUMBER] = 1;

    Parameters parameters = NumberOfLayers == 2 ? GetTwoLayersParameters() : GetThreeLayersParameters();

    // Read properties
    auto read_util = ReadMaterialsUtility(rThisModelPart.GetModel());
    read_util.ReadMaterials(parameters);

    // Create nodes and elements
    Properties::Pointer p_elem_prop = rThisModelPart.pGetProperties(1);

    // First we create the nodes
    rThisModelPart.CreateNewNode(1 , 0.0 , 1.0 , 1.0);
    rThisModelPart.CreateNewNode(2 , 0.0 , 1.0 , 0.0);
    rThisModelPart.CreateNewNode(3 , 0.0 , 0.0 , 1.0);
    rThisModelPart.CreateNewNode(4 , 1.0 , 1.0 , 1.0);
    rThisModelPart.CreateNewNode(5 , 0.0 , 0.0 , 0.0);
    rThisModelPart.CreateNewNode(6 , 1.0 , 1.0 , 0.0);

    rThisModelPart.CreateNewNode(7 , 1.0 , 0.0 , 1.0);
    rThisModelPart.CreateNewNode(8 , 1.0 , 0.0 , 0.0);
    rThisModelPart.CreateNewNode(9 , 2.0 , 1.0 , 1.0);
    rThisModelPart.CreateNewNode(10 , 2.0 , 1.0 , 0.0);
    rThisModelPart.CreateNewNode(11 , 2.0 , 0.0 , 1.0);
    rThisModelPart.CreateNewNode(12 , 2.0 , 0.0 , 0.0);

    // Now we create the elements
    rThisModelPart.CreateNewElement(ElementName, 1, {{12,10,8,9}}, p_elem_prop);
    rThisModelPart.CreateNewElement(ElementName, 2, {{4,6,9,7}}, p_elem_prop);
    rThisModelPart.CreateNewElement(ElementName, 3, {{11,7,9,8}}, p_elem_prop);
    rThisModelPart.CreateNewElement(ElementName, 4, {{5,3,8,6}}, p_elem_prop);
    rThisModelPart.CreateNewElement(ElementName, 5, {{4,6,7,3}}, p_elem_prop);
    rThisModelPart.CreateNewElement(ElementName, 6, {{2,3,5,6}}, p_elem_prop);
    rThisModelPart.CreateNewElement(ElementName, 7, {{10,9,6,8}}, p_elem_prop);
    rThisModelPart.CreateNewElement(ElementName, 8, {{7,8,3,6}}, p_elem_prop);
    rThisModelPart.CreateNewElement(ElementName, 9, {{7,8,6,9}}, p_elem_prop);
    rThisModelPart.CreateNewElement(ElementName, 10, {{4,1,6,3}}, p_elem_prop);
    rThisModelPart.CreateNewElement(ElementName, 11, {{9,12,11,8}}, p_elem_prop);
    rThisModelPart.CreateNewElement(ElementName, 12, {{3,2,1,6}}, p_elem_prop);

    const auto& const_process_info = rThisModelPart.GetProcessInfo();

    // Initialize elements
    for (auto& r_elem : rThisModelPart.Elements()) {
        r_elem.Initialize(const_process_info);
        r_elem.InitializeSolutionStep(const_process_info);
        r_elem.InitializeNonLinearIteration(const_process_info);
    }
}

/**
* Check the correct work of Rule of Mixtures (Hexahedron 2 layers)
*/
KRATOS_TEST_CASE_IN_SUITE(RuleOfMixturesConstitutiveLawHexahedronTwoLayers, KratosConstitutiveLawsWithStructuralElementsSuite)
{
    Model this_model;
    ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);
    Create3DGeometryHexahedraRuleOfMixtures(r_model_part);

    const array_1d<double, 3> zero = ZeroVector(3);
    array_1d<double, 3> delta = ZeroVector(3);
    delta[0] = 1.0e-2;
    for (auto& node : r_model_part.Nodes()) {
        if (node.X() < 1.0e-3) {
            node.Fix(DISPLACEMENT_X);
            node.Fix(DISPLACEMENT_Y);
            node.Fix(DISPLACEMENT_Z);
            node.FastGetSolutionStepValue(DISPLACEMENT) = zero;
        } else {
            node.FastGetSolutionStepValue(DISPLACEMENT) = delta;
            node.Coordinates() += delta;
        }
    }

    /// Tolerance
    const double tolerance = 1.0e-6;

    ProcessInfo& process_info = r_model_part.GetProcessInfo();
    for (auto& elem : r_model_part.Elements()) {
        std::vector<Vector> solution;
        elem.CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, solution, process_info);

        for (auto& sol : solution) {
            KRATOS_EXPECT_LE((sol[0] - 4.09597e+09)/4.09597e+09, tolerance);
        }
    }
}

/**
* Check the correct work of Rule of Mixtures (Hexahedron 3 layers)
*/
KRATOS_TEST_CASE_IN_SUITE(RuleOfMixturesConstitutiveLawHexahedronThreeLayers, KratosConstitutiveLawsWithStructuralElementsSuite)
{
    Model this_model;
    ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);
    Create3DGeometryHexahedraRuleOfMixtures(r_model_part, 3);

    const array_1d<double, 3> zero = ZeroVector(3);
    array_1d<double, 3> delta = ZeroVector(3);
    delta[0] = 1.0e-2;
    for (auto& node : r_model_part.Nodes()) {
        if (node.X() < 1.0e-3) {
            node.Fix(DISPLACEMENT_X);
            node.Fix(DISPLACEMENT_Y);
            node.Fix(DISPLACEMENT_Z);
            node.FastGetSolutionStepValue(DISPLACEMENT) = zero;
        } else {
            node.FastGetSolutionStepValue(DISPLACEMENT) = delta;
            node.Coordinates() += delta;
        }
    }

    /// Tolerance
    const double tolerance = 1.0e-6;

    ProcessInfo& process_info = r_model_part.GetProcessInfo();
    for (auto& elem : r_model_part.Elements()) {
        std::vector<Vector> solution;
        elem.CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, solution, process_info);

        for (auto& sol : solution) {
            KRATOS_EXPECT_LE((sol[0] - 2.32156e+12)/2.32156e+12, tolerance);
        }
    }
}

/**
* Check the correct work of Rule of Mixtures (Tetrahedron 2 layers)
*/
KRATOS_TEST_CASE_IN_SUITE(RuleOfMixturesConstitutiveLawTetrahedronTwoLayers, KratosConstitutiveLawsWithStructuralElementsSuite)
{
    Model this_model;
    ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);
    Create3DGeometryTetrahedraRuleOfMixtures(r_model_part);

    const array_1d<double, 3> zero = ZeroVector(3);
    array_1d<double, 3> delta = ZeroVector(3);
    delta[0] = 1.0e-2;
    for (auto& node : r_model_part.Nodes()) {
        if (node.X() < 1.0e-3) {
            node.Fix(DISPLACEMENT_X);
            node.Fix(DISPLACEMENT_Y);
            node.Fix(DISPLACEMENT_Z);
            node.FastGetSolutionStepValue(DISPLACEMENT) = zero;
        } else {
            node.FastGetSolutionStepValue(DISPLACEMENT) = delta;
            node.Coordinates() += delta;
        }
    }

    /// Tolerance
    const double tolerance = 1.0e-6;

    ProcessInfo& process_info = r_model_part.GetProcessInfo();
    for (auto& elem : r_model_part.Elements()) {
        std::vector<Vector> solution;
        elem.CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, solution, process_info);

        for (auto& sol : solution) {
            KRATOS_EXPECT_LE((sol[0] - 4.09597e+09)/4.09597e+09, tolerance);
        }
    }
}

/**
* Check the correct work of Rule of Mixtures (Tetrahedron 3 layers)
*/
KRATOS_TEST_CASE_IN_SUITE(RuleOfMixturesConstitutiveLawTetrahedronThreeLayers, KratosConstitutiveLawsWithStructuralElementsSuite)
{
    Model this_model;
    ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);
    Create3DGeometryTetrahedraRuleOfMixtures(r_model_part, 3);

    const array_1d<double, 3> zero = ZeroVector(3);
    array_1d<double, 3> delta = ZeroVector(3);
    delta[0] = 1.0e-2;
    for (auto& node : r_model_part.Nodes()) {
        if (node.X() < 1.0e-3) {
            node.Fix(DISPLACEMENT_X);
            node.Fix(DISPLACEMENT_Y);
            node.Fix(DISPLACEMENT_Z);
            node.FastGetSolutionStepValue(DISPLACEMENT) = zero;
        } else {
            node.FastGetSolutionStepValue(DISPLACEMENT) = delta;
            node.Coordinates() += delta;
        }
    }

    /// Tolerance
    const double tolerance = 1.0e-6;

    ProcessInfo& process_info = r_model_part.GetProcessInfo();
    for (auto& elem : r_model_part.Elements()) {
        std::vector<Vector> solution;
        elem.CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, solution, process_info);

        for (auto& sol : solution) {
            KRATOS_EXPECT_LE((sol[0] - 2.32156e+12)/2.32156e+12, tolerance);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(ParallelRuleOfMixturesClonePreservesIndependentLayerHistory, KratosConstitutiveLawsFastSuite)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("Main");
    const Tetrahedra3D4<Node> geometry(
        r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0),
        r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0),
        r_model_part.CreateNewNode(3, 0.0, 1.0, 0.0),
        r_model_part.CreateNewNode(4, 0.0, 0.0, 1.0));

    Properties mixture_properties(1);
    const auto& r_damage_law = KratosComponents<ConstitutiveLaw>::Get("SmallStrainIsotropicDamage3DVonMises");
    const auto make_layer_properties = [&](ModelPart::IndexType Id, double YieldStress) {
        auto p_properties = Kratos::make_shared<Properties>(Id);
        p_properties->SetValue(CONSTITUTIVE_LAW, r_damage_law.Clone());
        p_properties->SetValue(YOUNG_MODULUS, 210e9);
        p_properties->SetValue(POISSON_RATIO, 0.22);
        p_properties->SetValue(YIELD_STRESS_COMPRESSION, YieldStress);
        p_properties->SetValue(YIELD_STRESS_TENSION, YieldStress);
        p_properties->SetValue(FRACTURE_ENERGY, 1000.0);
        p_properties->SetValue(SOFTENING_TYPE, 0);
        mixture_properties.AddSubProperties(p_properties);
        return p_properties;
    };
    const auto p_first_layer = make_layer_properties(11, 3e6);
    make_layer_properties(12, 5e6);

    const auto& r_registered_mixture = KratosComponents<ConstitutiveLaw>::Get("ParallelRuleOfMixturesLaw3D");
    const auto p_source = r_registered_mixture.Create(Parameters(R"({"combination_factors": [1.0, 3.0]})"));
    p_source->InitializeMaterial(mixture_properties, geometry, ScalarVector(4, 0.25));
    KRATOS_EXPECT_EQ(p_source->Check(mixture_properties, geometry, r_model_part.GetProcessInfo()), 0);
    ConstitutiveLaw::Parameters values(geometry, mixture_properties, r_model_part.GetProcessInfo());
    Vector strain = ZeroVector(6);
    Vector stress = ZeroVector(6);
    Matrix tangent = ZeroMatrix(6, 6);
    values.SetStrainVector(strain);
    values.SetStressVector(stress);
    values.SetConstitutiveMatrix(tangent);
    values.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    values.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
    strain[0] = 8e-5;
    p_source->InitializeMaterialResponseCauchy(values);
    p_source->FinalizeMaterialResponseCauchy(values);

    const auto p_clone = p_source->Clone();
    auto* p_source_mixture = dynamic_cast<ParallelRuleOfMixturesLaw<3>*>(p_source.get());
    auto* p_clone_mixture = dynamic_cast<ParallelRuleOfMixturesLaw<3>*>(p_clone.get());
    ASSERT_NE(p_source_mixture, nullptr);
    ASSERT_NE(p_clone_mixture, nullptr);
    const auto& r_source_layers = p_source_mixture->GetConstitutiveLaws();
    const auto& r_cloned_layers = p_clone_mixture->GetConstitutiveLaws();
    ASSERT_EQ(r_source_layers.size(), 2);
    ASSERT_EQ(r_cloned_layers.size(), 2);
    const auto get_history = [](ConstitutiveLaw& rLaw) {
        Vector history;
        rLaw.GetValue(INTERNAL_VARIABLES, history);
        return history;
    };
    const Vector first_history = get_history(*r_source_layers[0]);
    const Vector second_history = get_history(*r_source_layers[1]);
    ASSERT_EQ(first_history.size(), 2);
    ASSERT_EQ(second_history.size(), 2);
    KRATOS_EXPECT_GT(first_history[0], second_history[0]);
    KRATOS_EXPECT_GT(second_history[0], 0.0);
    KRATOS_EXPECT_VECTOR_NEAR(get_history(*r_cloned_layers[0]), first_history, 1e-12);
    KRATOS_EXPECT_VECTOR_NEAR(get_history(*r_cloned_layers[1]), second_history, 1e-12);
    double damage = 0.0;
    p_clone->GetValue(DAMAGE, damage);
    KRATOS_EXPECT_NEAR(damage, 0.25 * first_history[0] + 0.75 * second_history[0], 1e-12);

    // Compare both mixtures under identical unloading before changing only one cloned layer.
    strain[0] = 4e-5;
    p_source->CalculateMaterialResponsePK2(values);
    const Vector source_stress = stress;
    p_clone->CalculateMaterialResponsePK2(values);
    KRATOS_EXPECT_VECTOR_NEAR(stress, source_stress, 1e-8);

    values.SetMaterialProperties(*p_first_layer);
    strain[0] = 1.6e-4;
    r_cloned_layers[0]->FinalizeMaterialResponsePK2(values);
    const Vector advanced_history = get_history(*r_cloned_layers[0]);
    KRATOS_EXPECT_GT(advanced_history[0], first_history[0]);
    KRATOS_EXPECT_VECTOR_NEAR(get_history(*r_cloned_layers[1]), second_history, 1e-12);
    KRATOS_EXPECT_VECTOR_NEAR(get_history(*r_source_layers[0]), first_history, 1e-12);
    KRATOS_EXPECT_VECTOR_NEAR(get_history(*r_source_layers[1]), second_history, 1e-12);
}

} // namespace Kratos::Testing
