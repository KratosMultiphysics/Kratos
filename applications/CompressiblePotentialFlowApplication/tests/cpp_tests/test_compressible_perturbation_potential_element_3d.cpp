//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Inigo Lopez
//
//

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "compressible_potential_flow_application_variables.h"
#include "fluid_dynamics_application_variables.h"
#include "custom_elements/compressible_perturbation_potential_flow_element.h"
#include "custom_utilities/potential_flow_utilities.h"
#include "tests/cpp_tests/test_utilities.h"

namespace Kratos {
namespace Testing {

typedef ModelPart::IndexType IndexType;
typedef ModelPart::NodeIterator NodeIteratorType;

void GenerateCompressiblePerturbationElement3D(ModelPart& rModelPart) {
    // Variables addition
    rModelPart.AddNodalSolutionStepVariable(VELOCITY_POTENTIAL);
    rModelPart.AddNodalSolutionStepVariable(AUXILIARY_VELOCITY_POTENTIAL);

    // Set the element properties
    Properties::Pointer pElemProp = rModelPart.CreateNewProperties(0);
    rModelPart.GetProcessInfo()[FREE_STREAM_DENSITY] = 1.225;
    rModelPart.GetProcessInfo()[FREE_STREAM_MACH] = 0.6;
    rModelPart.GetProcessInfo()[HEAT_CAPACITY_RATIO] = 1.4;
    rModelPart.GetProcessInfo()[SOUND_VELOCITY] = 340.3;
    rModelPart.GetProcessInfo()[MACH_LIMIT] = 0.94;

    BoundedVector<double, 3> free_stream_velocity = ZeroVector(3);
    free_stream_velocity(0) = rModelPart.GetProcessInfo().GetValue(FREE_STREAM_MACH) *
                              rModelPart.GetProcessInfo().GetValue(SOUND_VELOCITY);
    rModelPart.GetProcessInfo()[FREE_STREAM_VELOCITY] = free_stream_velocity;

    BoundedVector<double, 3> free_stream_velocity_direction = ZeroVector(3);
    free_stream_velocity_direction(0) = 1.0;
    rModelPart.GetProcessInfo()[FREE_STREAM_VELOCITY_DIRECTION] = free_stream_velocity_direction;

    BoundedVector<double, 3> wake_normal = ZeroVector(3);
    wake_normal(2) = 1.0;
    rModelPart.GetProcessInfo()[WAKE_NORMAL] = wake_normal;

    // Geometry creation
    rModelPart.CreateNewNode(1, 0.0, -0.2, -0.2);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 0.1, 1.0, 0.0);
    rModelPart.CreateNewNode(4, -0.1, 0.0, 1.0);
    std::vector<ModelPart::IndexType> elemNodes{1, 2, 3, 4};
    rModelPart.CreateNewElement(
        "CompressiblePerturbationPotentialFlowElement3D4N", 1, elemNodes, pElemProp);
}

KRATOS_TEST_CASE_IN_SUITE(PingCompressiblePerturbationPotentialFlowElementLHS3D,
                          CompressiblePotentialApplicationFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement3D(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    std::array<double, 4> potential{1.39572, 110.69275, 121.1549827, 104.284736};

    Matrix LHS_finite_diference = ZeroMatrix(number_of_nodes, number_of_nodes);
    Matrix LHS_analytical = ZeroMatrix(number_of_nodes, number_of_nodes);

    PotentialFlowTestUtilities::ComputeElementalSensitivities<4>(
        model_part, LHS_finite_diference, LHS_analytical, potential);

    KRATOS_CHECK_MATRIX_NEAR(LHS_finite_diference, LHS_analytical, 1e-10);
}

KRATOS_TEST_CASE_IN_SUITE(PingCompressiblePerturbationPotentialFlowElementLHS3DClamping,
                          CompressiblePotentialApplicationFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement3D(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    std::array<double, 4> potential{1.39572, 117.69275, 121.1549827, 104.284736};

    Matrix LHS_finite_diference = ZeroMatrix(number_of_nodes, number_of_nodes);
    Matrix LHS_analytical = ZeroMatrix(number_of_nodes, number_of_nodes);

    PotentialFlowTestUtilities::ComputeElementalSensitivities<4>(
        model_part, LHS_finite_diference, LHS_analytical, potential);

    KRATOS_CHECK_MATRIX_NEAR(LHS_finite_diference, LHS_analytical, 1e-10);
}

KRATOS_TEST_CASE_IN_SUITE(PingWakeCompressiblePerturbationPotentialFlowElementLHS3D,
                          CompressiblePotentialApplicationFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement3D(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    std::array<double, 8> potential{1.39572, 110.69275, 121.1549827, 104.284736,
                                    2.39572, 46.69275,  100.1549827, 102.284736};

    Matrix LHS_finite_diference = ZeroMatrix(2 * number_of_nodes, 2 * number_of_nodes);
    Matrix LHS_analytical = ZeroMatrix(2 * number_of_nodes, 2 * number_of_nodes);

    PotentialFlowTestUtilities::ComputeWakeElementalSensitivities<4>(
        model_part, LHS_finite_diference, LHS_analytical, potential);

    KRATOS_CHECK_MATRIX_NEAR(LHS_finite_diference, LHS_analytical, 1e-10);
}

KRATOS_TEST_CASE_IN_SUITE(PingWakeCompressiblePerturbationPotentialFlowElementLHS3DClamping,
                          CompressiblePotentialApplicationFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement3D(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    std::array<double, 8> potential{1.39572, 117.69275, 121.1549827, 104.284736,
                                    2.39572, 146.69275, 100.1549827, 102.284736};

    Matrix LHS_finite_diference = ZeroMatrix(2 * number_of_nodes, 2 * number_of_nodes);
    Matrix LHS_analytical = ZeroMatrix(2 * number_of_nodes, 2 * number_of_nodes);

    PotentialFlowTestUtilities::ComputeWakeElementalSensitivities<4>(
        model_part, LHS_finite_diference, LHS_analytical, potential);

    KRATOS_CHECK_MATRIX_NEAR(LHS_finite_diference, LHS_analytical, 1e-10);
}

KRATOS_TEST_CASE_IN_SUITE(PingWakeStructureCompressiblePerturbationPotentialFlowElementLHS3D,
                          CompressiblePotentialApplicationFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement3D(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    p_element->Set(STRUCTURE);
    p_element->GetGeometry()[number_of_nodes-1].SetValue(TRAILING_EDGE, true);

    std::array<double, 8> potential{1.39572, 110.69275, 121.1549827, 104.284736,
                                    2.39572,  46.69275, 100.1549827, 102.284736};

    Matrix LHS_finite_diference = ZeroMatrix(2 * number_of_nodes, 2 * number_of_nodes);
    Matrix LHS_analytical = ZeroMatrix(2 * number_of_nodes, 2 * number_of_nodes);

    PotentialFlowTestUtilities::ComputeWakeElementalSensitivities<4>(
        model_part, LHS_finite_diference, LHS_analytical, potential);

    KRATOS_CHECK_MATRIX_NEAR(LHS_finite_diference, LHS_analytical, 1e-10);
}

KRATOS_TEST_CASE_IN_SUITE(PingWakeStructureCompressiblePerturbationPotentialFlowElementLHS3DClamping,
                          CompressiblePotentialApplicationFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement3D(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    p_element->Set(STRUCTURE);
    p_element->GetGeometry()[number_of_nodes-1].SetValue(TRAILING_EDGE, true);

    std::array<double, 8> potential{1.39572, 117.69275, 121.1549827, 104.284736,
                                    2.39572, 146.69275, 100.1549827, 102.284736};

    Matrix LHS_finite_diference = ZeroMatrix(2 * number_of_nodes, 2 * number_of_nodes);
    Matrix LHS_analytical = ZeroMatrix(2 * number_of_nodes, 2 * number_of_nodes);

    PotentialFlowTestUtilities::ComputeWakeElementalSensitivities<4>(
        model_part, LHS_finite_diference, LHS_analytical, potential);

    KRATOS_CHECK_MATRIX_NEAR(LHS_finite_diference, LHS_analytical, 1e-10);
}

KRATOS_TEST_CASE_IN_SUITE(CompressiblePerturbationPotentialFlowElementRHS3D,
                          CompressiblePotentialApplicationFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement3D(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);

    std::array<double, 4> potential{1.39572, 110.69275, 121.1549827, 104.284736};
    PotentialFlowTestUtilities::AssignPotentialsToNormalElement<4>(*p_element, potential);

    // Compute RHS
    Vector RHS = ZeroVector(3);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->CalculateRightHandSide(RHS, r_current_process_info);

    std::vector<double> reference{71.66991905097665, -64.11826564927853,
                                  -3.932086180475159, -3.619567221222969};

    KRATOS_CHECK_VECTOR_NEAR(RHS, reference, 1e-13);
}

KRATOS_TEST_CASE_IN_SUITE(CompressiblePerturbationPotentialFlowElementLHS3D,
                          CompressiblePotentialApplicationFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement3D(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);

    std::array<double, 4> potential{1.39572, 110.69275, 121.1549827, 104.284736};
    PotentialFlowTestUtilities::AssignPotentialsToNormalElement<4>(*p_element,potential);

    // Compute LHS
    Matrix LHS = ZeroMatrix(3, 3);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->CalculateLeftHandSide(LHS, r_current_process_info);

    std::vector<double> reference{
        0.1376306838326365,   0.0248823393106424,  -0.06452385209048064,
        -0.09798917105279824, 0.0248823393106424,  0.06159497346557513,
        -0.0664293742338969,  -0.0200479385423207, -0.06452385209048064,
        -0.0664293742338969,  0.1825781102929537,  -0.0516248839685762,
        -0.09798917105279824, -0.0200479385423207, -0.0516248839685762,
        0.1696619935636952};

    for (unsigned int i = 0; i < LHS.size1(); i++) {
        for (unsigned int j = 0; j < LHS.size2(); j++) {
            KRATOS_CHECK_NEAR(LHS(i, j), reference[i * 4 + j], 1e-16);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(CompressiblePerturbationPotentialFlowElementLHS3DClamping,
                          CompressiblePotentialApplicationFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement3D(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);

    std::array<double, 4> potential{1.39572, 117.69275, 121.1549827, 104.284736};
    PotentialFlowTestUtilities::AssignPotentialsToNormalElement<4>(*p_element,potential);

    // Compute LHS
    Matrix LHS = ZeroMatrix(3, 3);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->CalculateLeftHandSide(LHS, r_current_process_info);

    std::vector<double> reference{
        0.3488734111253458,   -0.1649639838036404,  -0.07578177407226058,
        -0.1081276532494449,  -0.1649639838036404,  0.2309495773250965,
        -0.05563491218475713, -0.01035068133669899, -0.07578177407226058,
        -0.05563491218475713, 0.1821535081663726,   -0.05073682190935494,
        -0.1081276532494449,  -0.01035068133669899, -0.05073682190935494,
        0.1692151564954988};

    for (unsigned int i = 0; i < LHS.size1(); i++) {
        for (unsigned int j = 0; j < LHS.size2(); j++) {
            KRATOS_CHECK_NEAR(LHS(i, j), reference[i * 4 + j], 1e-16);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(WakeCompressiblePerturbationPotentialFlowElementRHS3D,
                          CompressiblePotentialApplicationFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement3D(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);

    BoundedVector<double,4> distances = PotentialFlowTestUtilities::AssignDistancesToElement<4>();

    p_element->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
    p_element->GetValue(WAKE) = true;

    std::array<double, 8> potential{1.39572, 110.69275, 121.1549827, 104.284736,
                                    2.39572,  46.69275, 100.1549827, 102.284736};
    PotentialFlowTestUtilities::AssignPotentialsToWakeElement<4>(*p_element,distances,potential);

    // Compute RHS and LHS
    Vector RHS = ZeroVector(6);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->CalculateRightHandSide(RHS, r_current_process_info);

    std::vector<double> reference{11.25952380952381, -14.46333333333333,
                                  2.251904761904762, -3.619567221222969,
                                  68.65551596318301, -58.62766030853704,
                                  -4.30462713896052, -0.9519047619047626};

    KRATOS_CHECK_VECTOR_NEAR(RHS, reference, 1e-13);
}

KRATOS_TEST_CASE_IN_SUITE(WakeStructureCompressiblePerturbationPotentialFlowElementRHS3D,
                          CompressiblePotentialApplicationFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement3D(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    BoundedVector<double,4> distances = PotentialFlowTestUtilities::AssignDistancesToElement<4>();

    p_element->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
    p_element->GetValue(WAKE) = true;
    p_element->Set(STRUCTURE);
    p_element->GetGeometry()[number_of_nodes-1].SetValue(TRAILING_EDGE, true);

    std::array<double, 8> potential{1.39572, 110.69275, 121.1549827, 104.284736,
                                    2.39572,  46.69275, 100.1549827, 102.284736};
    PotentialFlowTestUtilities::AssignPotentialsToWakeElement<4>(*p_element,distances,potential);

    // Compute RHS and LHS
    Vector RHS = ZeroVector(6);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->CalculateRightHandSide(RHS, r_current_process_info);

    std::vector<double> reference{11.25952380952381, -14.46333333333333,
                                  2.251904761904762, -0.4524459026528712,
                                  68.65551596318301, -58.62766030853704,
                                  -4.30462713896052, -5.007824951224748};

    KRATOS_CHECK_VECTOR_NEAR(RHS, reference, 1e-13);
}

KRATOS_TEST_CASE_IN_SUITE(WakeCompressiblePerturbationPotentialFlowElementLHS3D,
                          CompressiblePotentialApplicationFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement3D(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);

    BoundedVector<double,4> distances = PotentialFlowTestUtilities::AssignDistancesToElement<4>();

    p_element->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
    p_element->GetValue(WAKE) = true;

    std::array<double, 8> potential{1.39572, 110.69275, 121.1549827, 104.284736,
                                    2.39572,  46.69275, 100.1549827, 102.284736};
    PotentialFlowTestUtilities::AssignPotentialsToWakeElement<4>(*p_element,distances,potential);

    // Compute LHS
    Matrix LHS = ZeroMatrix(6, 6);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->CalculateLeftHandSide(LHS, r_current_process_info);

    // Check the LHS values
    std::vector<double> reference{0.263095238095238,-0.185,0.05261904761904762,-0.1307142857142857,
    -0.263095238095238,0.185,-0.05261904761904762,0.1307142857142857,-0.185,0.2356666666666666,
    -0.03699999999999999,-0.01366666666666666,0.185,-0.2356666666666666,0.03699999999999999,
    0.01366666666666666,0.05261904761904762,-0.03699999999999999,0.01052380952380952,
    -0.02614285714285715,-0.05261904761904762,0.03699999999999999,-0.01052380952380952,
    0.02614285714285715,-0.09798917105279824,-0.0200479385423207,-0.0516248839685762,
    0.1696619935636952,0,0,0,0,0,0,0,0,0.2513125629816388,-0.05867500741543832,-0.07898156406683947,
    -0.1136559914993611,0,0,0,0,-0.05867500741543832,0.1557535754651047,-0.07370192884392225,
    -0.0233766392057441,0,0,0,0,-0.07898156406683947,-0.07370192884392225,0.2130140967141911,
    -0.06033060380342937,0.1307142857142857,0.01366666666666666,0.02614285714285715,
    -0.1705238095238095,-0.1307142857142857,-0.01366666666666666,-0.02614285714285715,
    0.1705238095238095};

    for (unsigned int i = 0; i < LHS.size1(); i++)
    {
        for (unsigned int j = 0; j < LHS.size2(); j++)
        {
            KRATOS_CHECK_NEAR(LHS(i, j), reference[8 * i + j], 1e-16);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(WakeCompressiblePerturbationPotentialFlowElementLHS3DClamping,
                          CompressiblePotentialApplicationFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement3D(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);

    BoundedVector<double,4> distances = PotentialFlowTestUtilities::AssignDistancesToElement<4>();

    p_element->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
    p_element->GetValue(WAKE) = true;

    std::array<double, 8> potential{1.39572, 117.69275, 121.1549827, 104.284736,
                                    2.39572, 146.69275, 100.1549827, 102.284736};
    PotentialFlowTestUtilities::AssignPotentialsToWakeElement<4>(*p_element,distances,potential);

    // Compute LHS
    Matrix LHS = ZeroMatrix(6, 6);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->CalculateLeftHandSide(LHS, r_current_process_info);

    // Check the LHS values
    std::vector<double> reference{0.263095238095238,-0.185,0.05261904761904762,
    -0.1307142857142857,-0.263095238095238,0.185,-0.05261904761904762,
    0.1307142857142857,-0.185,0.2356666666666666,-0.03699999999999999,
    -0.01366666666666666,0.185,-0.2356666666666666,0.03699999999999999,
    0.01366666666666666,0.05261904761904762,-0.03699999999999999,
    0.01052380952380952,-0.02614285714285715,-0.05261904761904762,
    0.03699999999999999,-0.01052380952380952,0.02614285714285715,
    -0.1081276532494449,-0.01035068133669899,-0.05073682190935494,
    0.1692151564954988,0,0,0,0,0,0,0,0,0.3488734111253458,-0.1649639838036404,
    -0.07578177407226058,-0.1081276532494449,0,0,0,0,-0.1649639838036404,
    0.2309495773250965,-0.05563491218475713,-0.01035068133669899,0,0,0,0,
    -0.07578177407226058,-0.05563491218475713,0.1821535081663726,
    -0.05073682190935494,0.1307142857142857,0.01366666666666666,
    0.02614285714285715,-0.1705238095238095,-0.1307142857142857,
    -0.01366666666666666,-0.02614285714285715,0.1705238095238095};

    for (unsigned int i = 0; i < LHS.size1(); i++)
    {
        for (unsigned int j = 0; j < LHS.size2(); j++)
        {
            KRATOS_CHECK_NEAR(LHS(i, j), reference[8 * i + j], 1e-16);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(WakeStructureCompressiblePerturbationPotentialFlowElementLHS3D,
                          CompressiblePotentialApplicationFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement3D(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    BoundedVector<double,4> distances = PotentialFlowTestUtilities::AssignDistancesToElement<4>();

    p_element->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
    p_element->GetValue(WAKE) = true;
    p_element->Set(STRUCTURE);
    p_element->GetGeometry()[number_of_nodes-1].SetValue(TRAILING_EDGE, true);

    std::array<double, 8> potential{1.39572, 110.69275, 121.1549827, 104.284736,
                                    2.39572,  46.69275, 100.1549827, 102.284736};
    PotentialFlowTestUtilities::AssignPotentialsToWakeElement<4>(*p_element,distances,potential);

    // Compute LHS
    Matrix LHS = ZeroMatrix(6, 6);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->CalculateLeftHandSide(LHS, r_current_process_info);

    // Check the LHS values
    std::vector<double> reference{0.263095238095238,-0.185,0.05261904761904762,-0.1307142857142857,
    -0.263095238095238,0.185,-0.05261904761904762,0.1307142857142857,-0.185,0.2356666666666666,
    -0.03699999999999999,-0.01366666666666666,0.185,-0.2356666666666666,0.03699999999999999,
    0.01366666666666666,0.05261904761904762,-0.03699999999999999,0.01052380952380952,
    -0.02614285714285715,-0.05261904761904762,0.03699999999999999,-0.01052380952380952,
    0.02614285714285715,-0.01224864638159978,-0.002505992317790087,-0.006453110496072026,
    0.02120774919546189,0,0,0,0,0,0,0,0,0.2513125629816388,-0.05867500741543832,
    -0.07898156406683947,-0.1136559914993611,0,0,0,0,-0.05867500741543832,0.1557535754651047,
    -0.07370192884392225,-0.0233766392057441,0,0,0,0,-0.07898156406683947,-0.07370192884392225,
    0.2130140967141911,-0.06033060380342937,0,0,0,0,-0.09944899256194101,-0.0204545593050261,
    -0.05278927832800071,0.1726928301949679};

    for (unsigned int i = 0; i < LHS.size1(); i++)
    {
        for (unsigned int j = 0; j < LHS.size2(); j++)
        {
            KRATOS_CHECK_NEAR(LHS(i, j), reference[8 * i + j], 1e-16);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(WakeStructureCompressiblePerturbationPotentialFlowElementLHS3DClamping,
                          CompressiblePotentialApplicationFastSuite)
{
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressiblePerturbationElement3D(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    BoundedVector<double,4> distances = PotentialFlowTestUtilities::AssignDistancesToElement<4>();

    p_element->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
    p_element->GetValue(WAKE) = true;
    p_element->Set(STRUCTURE);
    p_element->GetGeometry()[number_of_nodes-1].SetValue(TRAILING_EDGE, true);

    std::array<double, 8> potential{1.39572, 117.69275, 121.1549827, 104.284736,
                                    2.39572, 146.69275, 100.1549827, 102.284736};
    PotentialFlowTestUtilities::AssignPotentialsToWakeElement<4>(*p_element,distances,potential);

    // Compute LHS
    Matrix LHS = ZeroMatrix(6, 6);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->CalculateLeftHandSide(LHS, r_current_process_info);

    // Check the LHS values
    std::vector<double> reference{0.263095238095238,-0.185,0.05261904761904762,
    -0.1307142857142857,-0.263095238095238,0.185,-0.05261904761904762,
    0.1307142857142857,-0.185,0.2356666666666666,-0.03699999999999999,
    -0.01366666666666666,0.185,-0.2356666666666666,0.03699999999999999,
    0.01366666666666666,0.05261904761904762,-0.03699999999999999,
    0.01052380952380952,-0.02614285714285715,-0.05261904761904762,
    0.03699999999999999,-0.01052380952380952,0.02614285714285715,
    -0.01351595665618061,-0.001293835167087374,-0.006342102738669368,
    0.02115189456193736,0,0,0,0,0,0,0,0,0.3488734111253458,-0.1649639838036404,
    -0.07578177407226058,-0.1081276532494449,0,0,0,0,-0.1649639838036404,
    0.2309495773250965,-0.05563491218475713,-0.01035068133669899,0,0,0,0,
    -0.07578177407226058,-0.05563491218475713,0.1821535081663726,
    -0.05073682190935494,0,0,0,0,-0.09461169659326432,-0.009056846169611622,
    -0.04439471917068559,0.1480632619335615};

    for (unsigned int i = 0; i < LHS.size1(); i++)
    {
        for (unsigned int j = 0; j < LHS.size2(); j++)
        {
            KRATOS_CHECK_NEAR(LHS(i, j), reference[8 * i + j], 1e-16);
        }
    }
}

} // namespace Testing
} // namespace Kratos.
