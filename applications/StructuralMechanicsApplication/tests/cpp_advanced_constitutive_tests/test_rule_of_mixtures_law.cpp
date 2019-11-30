// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                     license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "containers/model.h"
// #include "includes/gid_io.h"
#include "utilities/read_materials_utility.h"
#include "custom_advanced_constitutive/rule_of_mixtures_law.h"
#include "includes/mat_variables.h"

namespace Kratos
{
namespace Testing
{
/// Nodetype definition
typedef Node<3> NodeType;

// void GiDIODebugRuleMixtures(ModelPart& ThisModelPart)
// {
//     GidIO<> gid_io("TEST_RULE_MIXTURES", GiD_PostBinary, SingleFile, WriteUndeformed,  WriteElementsOnly);
//     const int nl_iter = ThisModelPart.GetProcessInfo()[NL_ITERATION_NUMBER];
//     const double label = static_cast<double>(nl_iter);
//
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
                "name" : "RuleOfMixturesLaw",
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
                "name" : "RuleOfMixturesLaw",
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

void Create3DGeometryHexahedraRuleOfMixtures(ModelPart& rThisModelPart, std::size_t NumberOfLayers = 2, const std::string ElementName = "SmallDisplacementElement3D8N")
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

    // Initialize elements
    for (auto& r_elem : rThisModelPart.Elements()) {
        r_elem.Initialize();
        r_elem.InitializeSolutionStep(r_process_info);
        r_elem.InitializeNonLinearIteration(r_process_info);
    }
}

void Create3DGeometryTetrahedraRuleOfMixtures(ModelPart& rThisModelPart, std::size_t NumberOfLayers = 2, const std::string ElementName = "SmallDisplacementElement3D4N")
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

    // Initialize elements
    for (auto& r_elem : rThisModelPart.Elements()) {
        r_elem.Initialize();
        r_elem.InitializeSolutionStep(r_process_info);
        r_elem.InitializeNonLinearIteration(r_process_info);
    }
}

/**
* Check the correct work of Rule of Mixtures (Hexahedron 2 layers)
*/
KRATOS_TEST_CASE_IN_SUITE(RuleOfMixturesConstitutiveLawHexahedronTwoLayers, KratosStructuralMechanicsFastSuite)
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

//     // DEBUG
//     GiDIODebugRuleMixtures(model_part);

    /// Tolerance
    const double tolerance = 1.0e-6;

    ProcessInfo& process_info = r_model_part.GetProcessInfo();
    for (auto& elem : r_model_part.Elements()) {
        std::vector<Vector> solution;
        elem.CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, solution, process_info);

        for (auto& sol : solution) {
            KRATOS_CHECK_LESS_EQUAL((sol[0] - 4.09597e+09)/4.09597e+09, tolerance);
        }
    }
}

/**
* Check the correct work of Rule of Mixtures (Hexahedron 3 layers)
*/
KRATOS_TEST_CASE_IN_SUITE(RuleOfMixturesConstitutiveLawHexahedronThreeLayers, KratosStructuralMechanicsFastSuite)
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

//     // DEBUG
//     GiDIODebugRuleMixtures(model_part);

    /// Tolerance
    const double tolerance = 1.0e-6;

    ProcessInfo& process_info = r_model_part.GetProcessInfo();
    for (auto& elem : r_model_part.Elements()) {
        std::vector<Vector> solution;
        elem.CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, solution, process_info);

        for (auto& sol : solution) {
            KRATOS_CHECK_LESS_EQUAL((sol[0] - 2.32156e+12)/2.32156e+12, tolerance);
        }
    }
}

/**
* Check the correct work of Rule of Mixtures (Tetrahedron 2 layers)
*/
KRATOS_TEST_CASE_IN_SUITE(RuleOfMixturesConstitutiveLawTetrahedronTwoLayers, KratosStructuralMechanicsFastSuite)
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

//     // DEBUG
//     GiDIODebugRuleMixtures(model_part);

    /// Tolerance
    const double tolerance = 1.0e-6;

    ProcessInfo& process_info = r_model_part.GetProcessInfo();
    for (auto& elem : r_model_part.Elements()) {
        std::vector<Vector> solution;
        elem.CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, solution, process_info);

        for (auto& sol : solution) {
            KRATOS_CHECK_LESS_EQUAL((sol[0] - 4.09597e+09)/4.09597e+09, tolerance);
        }
    }
}

/**
* Check the correct work of Rule of Mixtures (Tetrahedron 3 layers)
*/
KRATOS_TEST_CASE_IN_SUITE(RuleOfMixturesConstitutiveLawTetrahedronThreeLayers, KratosStructuralMechanicsFastSuite)
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

//     // DEBUG
//     GiDIODebugRuleMixtures(model_part);

    /// Tolerance
    const double tolerance = 1.0e-6;

    ProcessInfo& process_info = r_model_part.GetProcessInfo();
    for (auto& elem : r_model_part.Elements()) {
        std::vector<Vector> solution;
        elem.CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, solution, process_info);

        for (auto& sol : solution) {
            KRATOS_CHECK_LESS_EQUAL((sol[0] - 2.32156e+12)/2.32156e+12, tolerance);
        }
    }
}
} // namespace Testing
} // namespace Kratos
