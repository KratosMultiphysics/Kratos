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
#include "includes/gid_io.h"
#include "utilities/read_materials_utility.h"
#include "custom_constitutive/rule_of_mixtures_law.h"
#include "includes/mat_variables.h"

namespace Kratos
{
namespace Testing
{
/// Nodetype definition
typedef Node<3> NodeType;

void GiDIODebugRuleMixtures(ModelPart& ThisModelPart)
{
    GidIO<> gid_io("TEST_RULE_MIXTURES", GiD_PostBinary, SingleFile, WriteUndeformed,  WriteElementsOnly);
    const int nl_iter = ThisModelPart.GetProcessInfo()[NL_ITERATION_NUMBER];
    const double label = static_cast<double>(nl_iter);

    gid_io.InitializeMesh(label);
    gid_io.WriteMesh(ThisModelPart.GetMesh());
    gid_io.FinalizeMesh();
    gid_io.InitializeResults(label, ThisModelPart.GetMesh());
    gid_io.WriteNodalResults(DISPLACEMENT, ThisModelPart.Nodes(), label, 0);
    gid_io.PrintOnGaussPoints(GREEN_LAGRANGE_STRAIN_VECTOR, ThisModelPart, label);
    gid_io.PrintOnGaussPoints(PK2_STRESS_VECTOR, ThisModelPart, label);
    gid_io.WriteNodalFlags(ACTIVE, "ACTIVE", ThisModelPart.Nodes(), label);
}

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
                "sub_properties_indexes"   : [ 11,  13 ],
                "combination_factors"      : [0.4, 0.6 ],
                "material_rotation_angles" : [ 0,  90 ]
            },
            "Variables"        : {
                "THICKNESS"     : null,
                "DENSITY"       : null,
                "YOUNG_MODULUS" : null,
                "POISSON_RATIO" : null
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
            "properties_id"   : 13,
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
                "sub_properties_indexes"   : [ 11,  12,  13 ],
                "combination_factors"      : [0.4, 0.3, 0.3 ],
                "material_rotation_angles" : [ 0,   45,  90 ]
            },
            "Variables"        : {
                "THICKNESS"     : null,
                "DENSITY"       : null,
                "YOUNG_MODULUS" : null,
                "POISSON_RATIO" : null
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

void Create3DGeometryHexahedra(ModelPart& rThisModelPart, std::size_t NumberOfLayers = 2, const std::string ElementName = "SmallDisplacementElement3D8N")
{
    Model this_model;
    this_model.AddModelPart(&rThisModelPart);

    rThisModelPart.SetBufferSize(2);

    rThisModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);

    auto& process_info = rThisModelPart.GetProcessInfo();
    process_info[STEP] = 1;
    process_info[NL_ITERATION_NUMBER] = 1;

    Parameters parameters = NumberOfLayers == 2 ? GetTwoLayersParameters() : GetThreeLayersParameters();

    // Read properties
    auto read_util = ReadMaterialsUtility(this_model);
    read_util.ReadMaterials(parameters);

    // Create nodes and elements
    Properties::Pointer p_elem_prop = rThisModelPart.pGetProperties(1);

    // First we create the nodes
    NodeType::Pointer p_node_1 = rThisModelPart.CreateNewNode(1 , 0.0 , 1.0 , 1.0);
    NodeType::Pointer p_node_2 = rThisModelPart.CreateNewNode(2 , 0.0 , 1.0 , 0.0);
    NodeType::Pointer p_node_3 = rThisModelPart.CreateNewNode(3 , 0.0 , 0.0 , 1.0);
    NodeType::Pointer p_node_4 = rThisModelPart.CreateNewNode(4 , 1.0 , 1.0 , 1.0);
    NodeType::Pointer p_node_5 = rThisModelPart.CreateNewNode(5 , 0.0 , 0.0 , 0.0);
    NodeType::Pointer p_node_6 = rThisModelPart.CreateNewNode(6 , 1.0 , 1.0 , 0.0);
    NodeType::Pointer p_node_7 = rThisModelPart.CreateNewNode(7 , 1.0 , 0.0 , 1.0);
    NodeType::Pointer p_node_8 = rThisModelPart.CreateNewNode(8 , 1.0 , 0.0 , 0.0);

    // Now we create the "conditions"
    std::vector<NodeType::Pointer> element_nodes (8);
    element_nodes[0] = p_node_5;
    element_nodes[1] = p_node_8;
    element_nodes[2] = p_node_6;
    element_nodes[3] = p_node_2;
    element_nodes[4] = p_node_3;
    element_nodes[5] = p_node_7;
    element_nodes[6] = p_node_4;
    element_nodes[7] = p_node_1;
    Hexahedra3D8 <NodeType> hexahedra( PointerVector<NodeType>{element_nodes} );

    Element::Pointer p_elem_0 = rThisModelPart.CreateNewElement(ElementName, 1, hexahedra, p_elem_prop);

    // Initialize elements
    for (auto& elem : rThisModelPart.Elements()) {
        elem.Initialize();
        elem.InitializeSolutionStep(process_info);
        elem.InitializeNonLinearIteration(process_info);
    }
}

void Create3DGeometryTetrahedra(ModelPart& rThisModelPart, std::size_t NumberOfLayers = 2, const std::string ElementName = "SmallDisplacementElement3D4N")
{
    Model this_model;
    this_model.AddModelPart(&rThisModelPart);

    rThisModelPart.SetBufferSize(2);

    rThisModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);

    auto& process_info = rThisModelPart.GetProcessInfo();
    process_info[STEP] = 1;
    process_info[NL_ITERATION_NUMBER] = 1;

    Parameters parameters = NumberOfLayers == 2 ? GetTwoLayersParameters() : GetThreeLayersParameters();

    // Read properties
    auto read_util = ReadMaterialsUtility(this_model);
    read_util.ReadMaterials(parameters);

    // Create nodes and elements
    Properties::Pointer p_elem_prop = rThisModelPart.pGetProperties(1);

    // First we create the nodes
    NodeType::Pointer p_node_1 = rThisModelPart.CreateNewNode(1 , 0.0 , 1.0 , 1.0);
    NodeType::Pointer p_node_2 = rThisModelPart.CreateNewNode(2 , 0.0 , 1.0 , 0.0);
    NodeType::Pointer p_node_3 = rThisModelPart.CreateNewNode(3 , 0.0 , 0.0 , 1.0);
    NodeType::Pointer p_node_4 = rThisModelPart.CreateNewNode(4 , 1.0 , 1.0 , 1.0);
    NodeType::Pointer p_node_5 = rThisModelPart.CreateNewNode(5 , 0.0 , 0.0 , 0.0);
    NodeType::Pointer p_node_6 = rThisModelPart.CreateNewNode(6 , 1.0 , 1.0 , 0.0);

    NodeType::Pointer p_node_7 = rThisModelPart.CreateNewNode(7 , 1.0 , 0.0 , 1.0);
    NodeType::Pointer p_node_8 = rThisModelPart.CreateNewNode(8 , 1.0 , 0.0 , 0.0);
    NodeType::Pointer p_node_9 = rThisModelPart.CreateNewNode(9 , 2.0 , 1.0 , 1.0);
    NodeType::Pointer p_node_10 = rThisModelPart.CreateNewNode(10 , 2.0 , 1.0 , 0.0);
    NodeType::Pointer p_node_11 = rThisModelPart.CreateNewNode(11 , 2.0 , 0.0 , 1.0);
    NodeType::Pointer p_node_12 = rThisModelPart.CreateNewNode(12 , 2.0 , 0.0 , 0.0);

    // Now we create the "conditions"
    std::vector<NodeType::Pointer> element_nodes_0 (4);
    element_nodes_0[0] = p_node_12;
    element_nodes_0[1] = p_node_10;
    element_nodes_0[2] = p_node_8;
    element_nodes_0[3] = p_node_9;
    Tetrahedra3D4 <NodeType> tetrahedra_0( PointerVector<NodeType>{element_nodes_0} );

    std::vector<NodeType::Pointer> element_nodes_1 (4);
    element_nodes_1[0] = p_node_4;
    element_nodes_1[1] = p_node_6;
    element_nodes_1[2] = p_node_9;
    element_nodes_1[3] = p_node_7;
    Tetrahedra3D4 <NodeType> tetrahedra_1( PointerVector<NodeType>{element_nodes_1} );

    std::vector<NodeType::Pointer> element_nodes_2 (4);
    element_nodes_2[0] = p_node_11;
    element_nodes_2[1] = p_node_7;
    element_nodes_2[2] = p_node_9;
    element_nodes_2[3] = p_node_8;
    Tetrahedra3D4 <NodeType> tetrahedra_2( PointerVector<NodeType>{element_nodes_2} );

    std::vector<NodeType::Pointer> element_nodes_3 (4);
    element_nodes_3[0] = p_node_5;
    element_nodes_3[1] = p_node_3;
    element_nodes_3[2] = p_node_8;
    element_nodes_3[3] = p_node_6;
    Tetrahedra3D4 <NodeType> tetrahedra_3( PointerVector<NodeType>{element_nodes_3} );

    std::vector<NodeType::Pointer> element_nodes_4 (4);
    element_nodes_4[0] = p_node_4;
    element_nodes_4[1] = p_node_6;
    element_nodes_4[2] = p_node_7;
    element_nodes_4[3] = p_node_3;
    Tetrahedra3D4 <NodeType> tetrahedra_4( PointerVector<NodeType>{element_nodes_4} );

    std::vector<NodeType::Pointer> element_nodes_5 (4);
    element_nodes_5[0] = p_node_2;
    element_nodes_5[1] = p_node_3;
    element_nodes_5[2] = p_node_5;
    element_nodes_5[3] = p_node_6;
    Tetrahedra3D4 <NodeType> tetrahedra_5( PointerVector<NodeType>{element_nodes_5} );

    std::vector<NodeType::Pointer> element_nodes_6 (4);
    element_nodes_6[0] = p_node_10;
    element_nodes_6[1] = p_node_9;
    element_nodes_6[2] = p_node_6;
    element_nodes_6[3] = p_node_8;
    Tetrahedra3D4 <NodeType> tetrahedra_6( PointerVector<NodeType>{element_nodes_6} );

    std::vector<NodeType::Pointer> element_nodes_7 (4);
    element_nodes_7[0] = p_node_7;
    element_nodes_7[1] = p_node_8;
    element_nodes_7[2] = p_node_3;
    element_nodes_7[3] = p_node_6;
    Tetrahedra3D4 <NodeType> tetrahedra_7( PointerVector<NodeType>{element_nodes_7} );

    std::vector<NodeType::Pointer> element_nodes_8 (4);
    element_nodes_8[0] = p_node_7;
    element_nodes_8[1] = p_node_8;
    element_nodes_8[2] = p_node_6;
    element_nodes_8[3] = p_node_9;
    Tetrahedra3D4 <NodeType> tetrahedra_8( PointerVector<NodeType>{element_nodes_8} );

    std::vector<NodeType::Pointer> element_nodes_9 (4);
    element_nodes_9[0] = p_node_4;
    element_nodes_9[1] = p_node_1;
    element_nodes_9[2] = p_node_6;
    element_nodes_9[3] = p_node_3;
    Tetrahedra3D4 <NodeType> tetrahedra_9( PointerVector<NodeType>{element_nodes_9} );

    std::vector<NodeType::Pointer> element_nodes_10 (4);
    element_nodes_10[0] = p_node_9;
    element_nodes_10[1] = p_node_12;
    element_nodes_10[2] = p_node_11;
    element_nodes_10[3] = p_node_8;
    Tetrahedra3D4 <NodeType> tetrahedra_10( PointerVector<NodeType>{element_nodes_10} );

    std::vector<NodeType::Pointer> element_nodes_11 (4);
    element_nodes_11[0] = p_node_3;
    element_nodes_11[1] = p_node_2;
    element_nodes_11[2] = p_node_1;
    element_nodes_11[3] = p_node_6;
    Tetrahedra3D4 <NodeType> tetrahedra_11( PointerVector<NodeType>{element_nodes_11} );

    Element::Pointer p_elem_0 = rThisModelPart.CreateNewElement(ElementName, 1, tetrahedra_0, p_elem_prop);
    Element::Pointer p_elem_1 = rThisModelPart.CreateNewElement(ElementName, 2, tetrahedra_1, p_elem_prop);
    Element::Pointer p_elem_2 = rThisModelPart.CreateNewElement(ElementName, 3, tetrahedra_2, p_elem_prop);
    Element::Pointer p_elem_3 = rThisModelPart.CreateNewElement(ElementName, 4, tetrahedra_3, p_elem_prop);
    Element::Pointer p_elem_4 = rThisModelPart.CreateNewElement(ElementName, 5, tetrahedra_4, p_elem_prop);
    Element::Pointer p_elem_5 = rThisModelPart.CreateNewElement(ElementName, 6, tetrahedra_5, p_elem_prop);
    Element::Pointer p_elem_6 = rThisModelPart.CreateNewElement(ElementName, 7, tetrahedra_6, p_elem_prop);
    Element::Pointer p_elem_7 = rThisModelPart.CreateNewElement(ElementName, 8, tetrahedra_7, p_elem_prop);
    Element::Pointer p_elem_8 = rThisModelPart.CreateNewElement(ElementName, 9, tetrahedra_8, p_elem_prop);
    Element::Pointer p_elem_9 = rThisModelPart.CreateNewElement(ElementName, 10, tetrahedra_9, p_elem_prop);
    Element::Pointer p_elem_10 = rThisModelPart.CreateNewElement(ElementName, 11, tetrahedra_10, p_elem_prop);
    Element::Pointer p_elem_11 = rThisModelPart.CreateNewElement(ElementName, 12, tetrahedra_11, p_elem_prop);

    // Initialize elements
    for (auto& elem : rThisModelPart.Elements()) {
        elem.Initialize();
        elem.InitializeSolutionStep(process_info);
        elem.InitializeNonLinearIteration(process_info);
    }
}

/**
* Check the correct work of Rule of Mixtures (Hexahedron 2 layers)
*/
KRATOS_TEST_CASE_IN_SUITE(RuleOfMixturesConstitutiveLawHexahedronTwoLayers, KratosStructuralMechanicsFastSuite)
{
    ModelPart model_part("Main");
    Create3DGeometryHexahedra(model_part);

    const array_1d<double, 3> zero = ZeroVector(3);
    array_1d<double, 3> delta = ZeroVector(3);
    delta[0] = 1.0e-2;
    for (auto& node : model_part.Nodes()) {
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

    ProcessInfo& process_info = model_part.GetProcessInfo();
    for (auto& elem : model_part.Elements()) {
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
    ModelPart model_part("Main");
    Create3DGeometryHexahedra(model_part, 3);

    const array_1d<double, 3> zero = ZeroVector(3);
    array_1d<double, 3> delta = ZeroVector(3);
    delta[0] = 1.0e-2;
    for (auto& node : model_part.Nodes()) {
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

    ProcessInfo& process_info = model_part.GetProcessInfo();
    for (auto& elem : model_part.Elements()) {
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
    ModelPart model_part("Main");
    Create3DGeometryTetrahedra(model_part);

    const array_1d<double, 3> zero = ZeroVector(3);
    array_1d<double, 3> delta = ZeroVector(3);
    delta[0] = 1.0e-2;
    for (auto& node : model_part.Nodes()) {
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

    ProcessInfo& process_info = model_part.GetProcessInfo();
    for (auto& elem : model_part.Elements()) {
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
    ModelPart model_part("Main");
    Create3DGeometryTetrahedra(model_part, 3);

    const array_1d<double, 3> zero = ZeroVector(3);
    array_1d<double, 3> delta = ZeroVector(3);
    delta[0] = 1.0e-2;
    for (auto& node : model_part.Nodes()) {
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

    ProcessInfo& process_info = model_part.GetProcessInfo();
    for (auto& elem : model_part.Elements()) {
        std::vector<Vector> solution;
        elem.CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, solution, process_info);

        for (auto& sol : solution) {
            KRATOS_CHECK_LESS_EQUAL((sol[0] - 2.32156e+12)/2.32156e+12, tolerance);
        }
    }
}
} // namespace Testing
} // namespace Kratos
