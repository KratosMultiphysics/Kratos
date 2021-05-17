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
#include "custom_elements/compressible_potential_flow_element.h"
#include "custom_elements/embedded_compressible_potential_flow_element.h"
#include "tests/cpp_tests/test_utilities.h"

namespace Kratos {
namespace Testing {

typedef ModelPart::IndexType IndexType;
typedef ModelPart::NodeIterator NodeIteratorType;

void GenerateCompressibleElement(ModelPart& rModelPart) {
    // Variables addition
    rModelPart.AddNodalSolutionStepVariable(VELOCITY_POTENTIAL);
    rModelPart.AddNodalSolutionStepVariable(AUXILIARY_VELOCITY_POTENTIAL);

    // Set the element properties
    Properties::Pointer pElemProp = rModelPart.CreateNewProperties(0);
    BoundedVector<double, 3> v_inf = ZeroVector(3);
    v_inf(0) = 34.0;

    rModelPart.GetProcessInfo()[FREE_STREAM_VELOCITY] = v_inf;
    rModelPart.GetProcessInfo()[FREE_STREAM_DENSITY] = 1.225;
    rModelPart.GetProcessInfo()[FREE_STREAM_MACH] = 0.1;
    rModelPart.GetProcessInfo()[HEAT_CAPACITY_RATIO] = 1.4;
    rModelPart.GetProcessInfo()[SOUND_VELOCITY] = 340.0;
    rModelPart.GetProcessInfo()[MACH_LIMIT] = 0.94;

    // Geometry creation
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    std::vector<ModelPart::IndexType> elemNodes{1, 2, 3};
    rModelPart.CreateNewElement("CompressiblePotentialFlowElement2D3N", 1, elemNodes, pElemProp);
}

void GenerateCompressibleEmbeddedElement(ModelPart& rModelPart) {
    // Variables addition
    rModelPart.AddNodalSolutionStepVariable(VELOCITY_POTENTIAL);
    rModelPart.AddNodalSolutionStepVariable(AUXILIARY_VELOCITY_POTENTIAL);
    rModelPart.AddNodalSolutionStepVariable(GEOMETRY_DISTANCE);

    // Set the element properties
    rModelPart.CreateNewProperties(0);
    Properties::Pointer pElemProp = rModelPart.pGetProperties(0);
    BoundedVector<double, 3> v_inf = ZeroVector(3);
    v_inf(0) = 34.0;

    rModelPart.GetProcessInfo()[FREE_STREAM_VELOCITY] = v_inf;
    rModelPart.GetProcessInfo()[FREE_STREAM_DENSITY] = 1.0;
    rModelPart.GetProcessInfo()[FREE_STREAM_MACH] = 0.1;
    rModelPart.GetProcessInfo()[HEAT_CAPACITY_RATIO] = 1.4;
    rModelPart.GetProcessInfo()[SOUND_VELOCITY] = 340.0;
    rModelPart.GetProcessInfo()[MACH_LIMIT] = 0.94;

    // Geometry creation
    rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
    rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
    rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
    std::vector<ModelPart::IndexType> elemNodes{1, 2, 3};
    rModelPart.CreateNewElement("EmbeddedCompressiblePotentialFlowElement2D3N", 1, elemNodes, pElemProp);
}

/** Checks the CompressiblePotentialFlowElement element.
 * Checks the LHS and RHS computation.
 */
KRATOS_TEST_CASE_IN_SUITE(CompressiblePotentialFlowElementRHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressibleElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);

    // Define the nodal values
    std::array<double, 3> potential;
    potential[0] = 1.0;
    potential[1] = 2.0;
    potential[2] = 3.0;

    for (unsigned int i = 0; i < 3; i++) {
        p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = potential[i];
    }
    // Compute RHS and LHS
    Vector RHS = ZeroVector(3);
    Matrix LHS = ZeroMatrix(3, 3);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->CalculateLocalSystem(LHS, RHS, r_current_process_info);

    std::vector<double> reference({0.615561780, 0.0, -0.615561780});

    for (unsigned int i = 0; i < RHS.size(); i++) {
        KRATOS_CHECK_NEAR(RHS(i), reference[i], 1e-6);
    }
}

/** Checks the CompressiblePotentialFlowElement element.
 * Checks the LHS and RHS computation.
 */
KRATOS_TEST_CASE_IN_SUITE(CompressiblePotentialFlowElementLHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressibleElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);

    // Define the nodal values
    std::array<double, 3> potential;
    potential[0] = 1.0;
    potential[1] = 2.0;
    potential[2] = 3.0;

    for (unsigned int i = 0; i < 3; i++) {
        p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = potential[i];
    }
    // Compute RHS and LHS
    Vector RHS = ZeroVector(3);
    Matrix LHS = ZeroMatrix(3, 3);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->CalculateLocalSystem(LHS, RHS, r_current_process_info);

    std::array<double, 9> reference{0.615556466, -0.615561780, 5.314318652e-06, -0.615561780, 1.231123561, -0.615561780, 5.314318652e-06, -0.615561780, 0.615556466};

    for (unsigned int i = 0; i < LHS.size1(); i++) {
        for (unsigned int j = 0; j < LHS.size2(); j++) {
            KRATOS_CHECK_NEAR(LHS(i, j), reference[i * 3 + j], 1e-6);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(PingCompressiblePotentialFlowElementLHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressibleElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    const std::array<double, 3> potential{1.0, 20.0, 50.0};

    Matrix LHS_finite_diference = ZeroMatrix(number_of_nodes, number_of_nodes);
    Matrix LHS_analytical = ZeroMatrix(number_of_nodes, number_of_nodes);

    PotentialFlowTestUtilities::ComputeElementalSensitivities<3>(model_part, LHS_finite_diference, LHS_analytical, potential);

    KRATOS_CHECK_MATRIX_NEAR(LHS_finite_diference, LHS_analytical, 1e-10);
}

KRATOS_TEST_CASE_IN_SUITE(PingCompressiblePotentialFlowElementLHSClamping, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressibleElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    std::array<double, 3> potential{1.2495, 794.1948, 582.149583};

    Matrix LHS_finite_diference = ZeroMatrix(number_of_nodes, number_of_nodes);
    Matrix LHS_analytical = ZeroMatrix(number_of_nodes, number_of_nodes);

    PotentialFlowTestUtilities::ComputeElementalSensitivities<3>(model_part, LHS_finite_diference, LHS_analytical, potential);

    KRATOS_CHECK_MATRIX_NEAR(LHS_finite_diference, LHS_analytical, 1e-10);
}

KRATOS_TEST_CASE_IN_SUITE(EmbeddedCompressiblePotentialFlowElementCalculateLocalSystemRHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressibleEmbeddedElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);

    // Define the nodal values
    std::array<double, 3> potential{1.0, 2.0, 3.0};
    // Define the distance values
    std::array<double, 3> level_set{1.0, -1.0, -1.0};

    for (unsigned int i = 0; i < 3; i++) {
        p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = potential[i];
        p_element->GetGeometry()[i].FastGetSolutionStepValue(GEOMETRY_DISTANCE) = level_set[i];
    }

    // Compute RHS and LHS
    Vector RHS = ZeroVector(3);
    Matrix LHS = ZeroMatrix(3, 3);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->CalculateLocalSystem(LHS, RHS, r_current_process_info);

    std::vector<double> reference{0.125625, 0.0, -0.125625};

    KRATOS_CHECK_VECTOR_NEAR(RHS, reference, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(EmbeddedCompressiblePotentialFlowElementCalculateLocalSystemLHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressibleEmbeddedElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);

    // Define the nodal values
    std::array<double, 3> potential{1.0, 2.0, 3.0};
    // Define the distance values
    std::array<double, 3> level_set{1.0, -1.0, -1.0};
    for (unsigned int i = 0; i < 3; i++) {
        p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = potential[i];
        p_element->GetGeometry()[i].FastGetSolutionStepValue(GEOMETRY_DISTANCE) = level_set[i];
    }

    // Compute RHS and LHS
    Vector RHS = ZeroVector(3);
    Matrix LHS = ZeroMatrix(3, 3);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->CalculateLocalSystem(LHS, RHS, r_current_process_info);

    std::array<double, 9> reference_array{0.125624, -0.125625, 1.08455e-06, -0.125625, 0.25125, -0.125625, 1.08455e-06, -0.125625, 0.125624};
    // Copying to a 3x3 matrix to check against LHS
    Matrix reference(3, 3);
    for (unsigned int i = 0; i < reference.size1(); i++) {
        for (unsigned int j = 0; j < reference.size2(); j++) {
            reference(i, j) = reference_array[i * reference.size1() + j];
        }
    }

    KRATOS_CHECK_MATRIX_NEAR(LHS, reference, 1e-6);
}

KRATOS_TEST_CASE_IN_SUITE(PingEmbeddedCompressiblePotentialFlowElementLHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressibleEmbeddedElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    // Define the distance values
    std::array<double, 3> level_set{1.0, -1.0, -1.0};
    for (unsigned int i = 0; i < 3; i++) {
        p_element->GetGeometry()[i].FastGetSolutionStepValue(GEOMETRY_DISTANCE) = level_set[i];
    }

    const std::array<double, 3> potential{1.0, 20.0, 50.0};

    Matrix LHS_finite_diference = ZeroMatrix(number_of_nodes, number_of_nodes);
    Matrix LHS_analytical = ZeroMatrix(number_of_nodes, number_of_nodes);

    PotentialFlowTestUtilities::ComputeElementalSensitivities<3>(model_part, LHS_finite_diference, LHS_analytical, potential);

    KRATOS_CHECK_MATRIX_NEAR(LHS_finite_diference, LHS_analytical, 1e-10);
}

KRATOS_TEST_CASE_IN_SUITE(PingEmbeddedCompressiblePotentialFlowElementLHSPenalty, CompressiblePotentialApplicationFastSuite) {
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateCompressibleEmbeddedElement(model_part);
      model_part.GetProcessInfo()[PENALTY_COEFFICIENT] = 100.0;
      Element::Pointer p_element = model_part.pGetElement(1);
      const unsigned int number_of_nodes = p_element->GetGeometry().size();

      // Define the distance values
      std::array<double, 3> level_set{1.0, -1.0, -1.0};
      for (unsigned int i = 0; i < 3; i++) {
          p_element->GetGeometry()[i].FastGetSolutionStepValue(GEOMETRY_DISTANCE) = level_set[i];
      }

      const std::array<double, 3> potential{1.0, 20.0, 50.0};

      Matrix LHS_finite_diference = ZeroMatrix(number_of_nodes, number_of_nodes);
      Matrix LHS_analytical = ZeroMatrix(number_of_nodes, number_of_nodes);

      PotentialFlowTestUtilities::ComputeElementalSensitivities<3>(model_part, LHS_finite_diference, LHS_analytical, potential);

      KRATOS_CHECK_MATRIX_NEAR(LHS_finite_diference, LHS_analytical, 1e-10);
}

KRATOS_TEST_CASE_IN_SUITE(PingWakeEmbeddedCompressiblePotentialFlowElementLHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressibleEmbeddedElement(model_part);
    model_part.GetProcessInfo()[PENALTY_COEFFICIENT] = 100.0;
    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();
    Vector distances(3);
    distances(0) = 1.0;
    distances(1) = -1.0;
    distances(2) = -1.0;
    p_element->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
    p_element->GetValue(WAKE) = true;

    // Define the distance values
    std::array<double, 3> level_set{1.0, -1.0, -1.0};
    for (unsigned int i = 0; i < 3; i++) {
        p_element->GetGeometry()[i].FastGetSolutionStepValue(GEOMETRY_DISTANCE) = level_set[i];
    }

    const std::array<double, 6> potential{1.0, 40.0, 35.0, 6.0, 26.0, 14.0};

    Matrix LHS_finite_diference = ZeroMatrix(2*number_of_nodes, 2*number_of_nodes);
    Matrix LHS_analytical = ZeroMatrix(2*number_of_nodes, 2*number_of_nodes);

    PotentialFlowTestUtilities::ComputeWakeElementalSensitivities<3>(model_part, LHS_finite_diference, LHS_analytical, potential);
    KRATOS_WATCH(LHS_analytical)
    KRATOS_CHECK_MATRIX_NEAR(LHS_finite_diference, LHS_analytical, 1e-10);
}


KRATOS_TEST_CASE_IN_SUITE(PingEmbeddedCompressiblePotentialFlowElementLHSClamping, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressibleEmbeddedElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    // Define the distance values
    std::array<double, 3> level_set{1.0, -1.0, -1.0};
    for (unsigned int i = 0; i < 3; i++) {
        p_element->GetGeometry()[i].FastGetSolutionStepValue(GEOMETRY_DISTANCE) = level_set[i];
    }

    std::array<double, 3> potential{1.2495, 794.1948, 582.149583};

    Matrix LHS_finite_diference = ZeroMatrix(number_of_nodes, number_of_nodes);
    Matrix LHS_analytical = ZeroMatrix(number_of_nodes, number_of_nodes);

    PotentialFlowTestUtilities::ComputeElementalSensitivities<3>(model_part, LHS_finite_diference, LHS_analytical, potential);

    KRATOS_CHECK_MATRIX_NEAR(LHS_finite_diference, LHS_analytical, 1e-10);
}

KRATOS_TEST_CASE_IN_SUITE(CompressiblePotentialFlowElementRHSWake, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressibleElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);

    // Define the nodal values
    std::array<double, 3> potential;
    potential[0] = 1.0;
    potential[1] = 2.0;
    potential[2] = 3.0;

    Vector distances(3);
    distances(0) = 1.0;
    distances(1) = -1.0;
    distances(2) = -1.0;

    p_element->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
    p_element->GetValue(WAKE) = true;

    for (unsigned int i = 0; i < 3; i++) {
        if (distances(i) > 0.0) {
            p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = potential[i];
        }
        else {
            p_element->GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) = potential[i];
        }
    }
    for (unsigned int i = 0; i < 3; i++) {
        if (distances(i) < 0.0) {
            p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = potential[i] + 5;
        }
        else {
            p_element->GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) = potential[i] + 5;
        }
    }

    // Compute RHS and LHS
    Vector RHS = ZeroVector(6);
    Matrix LHS = ZeroMatrix(6, 6);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->CalculateLocalSystem(LHS, RHS, r_current_process_info);

    std::array<double, 6> reference{0.615561780, 0.0, 0.0, 0.0, 0.0, -0.615561780};

    for (unsigned int i = 0; i < RHS.size(); i++) {
        KRATOS_CHECK_NEAR(RHS(i), reference[i], 1e-6);
    }
}

KRATOS_TEST_CASE_IN_SUITE(CompressiblePotentialFlowElementLHSWake, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressibleElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);

    // Define the nodal values
    std::array<double, 3> potential;
    potential[0] = 1.0;
    potential[1] = 2.0;
    potential[2] = 3.0;

    Vector distances(3);
    distances(0) = 1.0;
    distances(1) = -1.0;
    distances(2) = -1.0;

    p_element->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
    p_element->GetValue(WAKE) = true;

    for (unsigned int i = 0; i < 3; i++) {
        if (distances(i) > 0.0) {
            p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = potential[i];
        }
        else {
            p_element->GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) = potential[i];
        }
    }
    for (unsigned int i = 0; i < 3; i++) {
        if (distances(i) < 0.0) {
            p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = potential[i] + 5;
        }
        else {
            p_element->GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) = potential[i] + 5;
        }
    }

    // Compute RHS and LHS
    Vector RHS = ZeroVector(6);
    Matrix LHS = ZeroMatrix(6, 6);

    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->CalculateLocalSystem(LHS, RHS, r_current_process_info);

    // Check the RHS values (the RHS is computed as the LHS x previous_solution,
    // hence, it is assumed that if the RHS is correct, the LHS is correct as well)
    std::array<double,36> reference{0.6155564666297989,-0.6155617809484512,5.314318652279458e-06,0,0,0,
                                    -0.6125,1.225,-0.6125,0.6125,-1.225,0.6125,
                                    0,-0.6125,0.6125,-0,0.6125,-0.6125,
                                    -0.6125,0.6125,-0,0.6125,-0.6125,0,
                                    0,0,0,-0.6155617809484512,1.231123561896902,-0.6155617809484512,
                                    0,0,0,5.314318652279458e-06,-0.6155617809484512,0.6155564666297989};

    for (unsigned int i = 0; i < LHS.size1(); i++) {
        for (unsigned int j = 0; j < LHS.size2(); j++) {
            KRATOS_CHECK_NEAR(LHS(i, j), reference[6 * i + j], 1e-6);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(PingWakeCompressiblePotentialFlowElementLHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressibleElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    const std::array<double, 6> potential{1.0, 40.0, 35.0, 6.0, 26.0, 14.0};

    Matrix LHS_finite_diference = ZeroMatrix(2*number_of_nodes, 2*number_of_nodes);
    Matrix LHS_analytical = ZeroMatrix(2*number_of_nodes, 2*number_of_nodes);

    PotentialFlowTestUtilities::ComputeWakeElementalSensitivities<3>(model_part, LHS_finite_diference, LHS_analytical, potential);

    KRATOS_CHECK_MATRIX_NEAR(LHS_finite_diference, LHS_analytical, 1e-10);
}

KRATOS_TEST_CASE_IN_SUITE(PingWakeCompressiblePotentialFlowElementLHSClamping, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressibleElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    const std::array<double, 6> potential{1.285837, 570.29384, 635.1583, 6.0, 796.345, 814.254};

    Matrix LHS_finite_diference = ZeroMatrix(2*number_of_nodes, 2*number_of_nodes);
    Matrix LHS_analytical = ZeroMatrix(2*number_of_nodes, 2*number_of_nodes);

    PotentialFlowTestUtilities::ComputeWakeElementalSensitivities<3>(model_part, LHS_finite_diference, LHS_analytical, potential);

    KRATOS_CHECK_MATRIX_NEAR(LHS_finite_diference, LHS_analytical, 1e-10);
}

KRATOS_TEST_CASE_IN_SUITE(PingWakeStructureCompressiblePotentialFlowElementLHS, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressibleElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    p_element->Set(STRUCTURE);
    p_element->GetGeometry()[number_of_nodes-1].SetValue(TRAILING_EDGE, true);

    const std::array<double, 6> potential{1.285837, 30.29384, 35.1583, 6.0, 46.345, 64.0};

    Matrix LHS_finite_diference = ZeroMatrix(2*number_of_nodes, 2*number_of_nodes);
    Matrix LHS_analytical = ZeroMatrix(2*number_of_nodes, 2*number_of_nodes);

    PotentialFlowTestUtilities::ComputeWakeElementalSensitivities<3>(model_part, LHS_finite_diference, LHS_analytical, potential);

    KRATOS_CHECK_MATRIX_NEAR(LHS_finite_diference, LHS_analytical, 1e-10);
}

KRATOS_TEST_CASE_IN_SUITE(PingWakeStructureCompressiblePotentialFlowElementLHSClamping, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressibleElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    p_element->Set(STRUCTURE);
    p_element->GetGeometry()[number_of_nodes-1].SetValue(TRAILING_EDGE, true);

    const std::array<double, 6> potential{1.285837, 170.29384, 135.1583, 6.0, 196.345, 114.0};

    Matrix LHS_finite_diference = ZeroMatrix(2*number_of_nodes, 2*number_of_nodes);
    Matrix LHS_analytical = ZeroMatrix(2*number_of_nodes, 2*number_of_nodes);

    PotentialFlowTestUtilities::ComputeWakeElementalSensitivities<3>(model_part, LHS_finite_diference, LHS_analytical, potential);

    KRATOS_CHECK_MATRIX_NEAR(LHS_finite_diference, LHS_analytical, 1e-10);
}

/** Checks the CompressiblePotentialFlowElement element.
* Checks the EquationIdVector.
*/
KRATOS_TEST_CASE_IN_SUITE(CompressiblePotentialFlowElementEquationIdVector, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressibleElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);

    for (unsigned int i = 0; i < 3; i++) {
        p_element->GetGeometry()[i].AddDof(VELOCITY_POTENTIAL);
    }

    Element::DofsVectorType ElementalDofList;
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->GetDofList(ElementalDofList, r_current_process_info);

    for (int i = 0; i < 3; i++) {
        ElementalDofList[i]->SetEquationId(i);
    }

    Element::EquationIdVectorType EquationIdVector;
    p_element->EquationIdVector(EquationIdVector, r_current_process_info);

    // Check the EquationIdVector values
    for (unsigned int i = 0; i < EquationIdVector.size(); i++) {
        KRATOS_CHECK(EquationIdVector[i] == i);
    }
}

/** Checks the CompressiblePotentialFlowElement element.
* Checks the EquationIdVector for the Wake.
*/
KRATOS_TEST_CASE_IN_SUITE(CompressiblePotentialFlowElementEquationIdVectorWake, CompressiblePotentialApplicationFastSuite) {
    Model this_model;
    ModelPart& model_part = this_model.CreateModelPart("Main", 3);

    GenerateCompressibleElement(model_part);
    Element::Pointer p_element = model_part.pGetElement(1);
    p_element->SetValue(WAKE, true);

    Vector distances(3);
    distances(0) = -0.5;
    distances(1) = -0.5;
    distances(2) = 0.5;
    p_element->SetValue(WAKE_ELEMENTAL_DISTANCES, distances);

    for (unsigned int i = 0; i < 3; i++) {
        p_element->GetGeometry()[i].AddDof(VELOCITY_POTENTIAL);
        p_element->GetGeometry()[i].AddDof(AUXILIARY_VELOCITY_POTENTIAL);
    }

    Element::DofsVectorType ElementalDofList;
    const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
    p_element->GetDofList(ElementalDofList, r_current_process_info);

    for (int i = 0; i < 6; i++) {
        ElementalDofList[i]->SetEquationId(i);
    }

    Element::EquationIdVectorType EquationIdVector;
    p_element->EquationIdVector(EquationIdVector, r_current_process_info);

    // Check the EquationIdVector values
    for (unsigned int i = 0; i < EquationIdVector.size(); i++) {
        KRATOS_CHECK(EquationIdVector[i] == i);
    }
}

} // namespace Testing
} // namespace Kratos.
