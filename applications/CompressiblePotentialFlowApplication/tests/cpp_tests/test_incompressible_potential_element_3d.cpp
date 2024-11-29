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
#include "tests/cpp_tests/compressible_potential_flow_fast_suite.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "compressible_potential_flow_application_variables.h"
#include "custom_elements/incompressible_potential_flow_element.h"
#include "custom_utilities/potential_flow_utilities.h"
#include "tests/cpp_tests/test_utilities.h"

namespace Kratos {
  namespace Testing {

    typedef ModelPart::IndexType IndexType;
    typedef ModelPart::NodeIterator NodeIteratorType;

    void GenerateIncompressibleElement3D(ModelPart& rModelPart)
    {
      // Variables addition
      rModelPart.AddNodalSolutionStepVariable(VELOCITY_POTENTIAL);
      rModelPart.AddNodalSolutionStepVariable(AUXILIARY_VELOCITY_POTENTIAL);

      // Set the element properties
      rModelPart.CreateNewProperties(0);
      Properties::Pointer pElemProp = rModelPart.pGetProperties(0);
      BoundedVector<double, 3> free_stream_velocity = ZeroVector(3);
      free_stream_velocity(0) = 10.0;

      BoundedVector<double, 3> free_stream_velocity_direction = ZeroVector(3);
      free_stream_velocity_direction(0) = 1.0;

      BoundedVector<double, 3> wake_normal = ZeroVector(3);
      wake_normal(2) = 1.0;

      rModelPart.GetProcessInfo()[FREE_STREAM_VELOCITY] = free_stream_velocity;
      rModelPart.GetProcessInfo()[FREE_STREAM_DENSITY] = 1.0;
      rModelPart.GetProcessInfo()[FREE_STREAM_VELOCITY_DIRECTION] = free_stream_velocity_direction;
      rModelPart.GetProcessInfo()[WAKE_NORMAL] = wake_normal;

      // Geometry creation
      rModelPart.CreateNewNode(1, 0.0, -0.2, -0.2);
      rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
      rModelPart.CreateNewNode(3, 0.1, 1.0, 0.0);
      rModelPart.CreateNewNode(4, -0.1, 0.0, 1.0);
      std::vector<ModelPart::IndexType> elemNodes{1, 2, 3, 4};
      rModelPart.CreateNewElement("IncompressiblePotentialFlowElement3D4N", 1, elemNodes, pElemProp);
    }

    // Checks the LHS of the 3D IncompressiblePotentialFlowElement element.
    KRATOS_TEST_CASE_IN_SUITE(PingIncompressiblePotentialFlowElementLHS3D,
                              CompressiblePotentialApplicationFastSuite)
    {
        Model this_model;
        ModelPart& model_part = this_model.CreateModelPart("Main", 3);

        GenerateIncompressibleElement3D(model_part);
        Element::Pointer p_element = model_part.pGetElement(1);
        const unsigned int number_of_nodes = p_element->GetGeometry().size();

        std::array<double, 4> potential{1.39572, 143.39275, 151.1549827, 134.284736};
        Matrix LHS_finite_diference = ZeroMatrix(number_of_nodes, number_of_nodes);
        Matrix LHS_analytical = ZeroMatrix(number_of_nodes, number_of_nodes);

        PotentialFlowTestUtilities::ComputeElementalSensitivities<4>(
            model_part, LHS_finite_diference, LHS_analytical, potential);

        KRATOS_EXPECT_MATRIX_NEAR(LHS_finite_diference, LHS_analytical, 1e-10);
    }

    // Checks the LHS of the 3D Wake IncompressiblePotentialFlowElement element.
    KRATOS_TEST_CASE_IN_SUITE(PingWakeIncompressiblePotentialFlowElementLHS3D,
                              CompressiblePotentialApplicationFastSuite)
    {
        Model this_model;
        ModelPart& model_part = this_model.CreateModelPart("Main", 3);

        GenerateIncompressibleElement3D(model_part);
        Element::Pointer p_element = model_part.pGetElement(1);
        const unsigned int number_of_nodes = p_element->GetGeometry().size();

        std::array<double, 8> potential{1.39572,     110.69275, 121.1549827,
                                        104.284736,  2.39572,   46.69275,
                                        100.1549827, 102.284736};
        Matrix LHS_finite_diference = ZeroMatrix(2*number_of_nodes, 2*number_of_nodes);
        Matrix LHS_analytical = ZeroMatrix(2*number_of_nodes, 2*number_of_nodes);

        PotentialFlowTestUtilities::ComputeWakeElementalSensitivities<4>(
            model_part, LHS_finite_diference, LHS_analytical, potential);

        KRATOS_EXPECT_MATRIX_NEAR(LHS_finite_diference, LHS_analytical, 1e-10);
    }

    // Checks the RHS of the 3D Wake IncompressiblePotentialFlowElement element.
    KRATOS_TEST_CASE_IN_SUITE(WakeIncompressiblePotentialFlowElementRHS3D,
                              CompressiblePotentialApplicationFastSuite)
    {
        Model this_model;
        ModelPart& model_part = this_model.CreateModelPart("Main", 3);

        GenerateIncompressibleElement3D(model_part);
        Element::Pointer p_element = model_part.pGetElement(1);

        BoundedVector<double,4> distances = PotentialFlowTestUtilities::AssignDistancesToElement<4>();

        p_element->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
        p_element->GetValue(WAKE) = true;

        std::array<double, 8> potential{1.39572,     110.69275, 121.1549827,
                                        104.284736,  2.39572,   46.69275,
                                        100.1549827, 102.284736};
        PotentialFlowTestUtilities::AssignPotentialsToWakeElement<4>(*p_element,distances,potential);

        // Compute RHS
        Vector RHS = ZeroVector(4);

        const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
        p_element->CalculateRightHandSide(RHS, r_current_process_info);

        std::vector<double> reference{11.25952380952381,-14.46333333333333,2.251904761904762,-10.51435102035238,26.29551835085714,-3.872345907866666,-10.5883452321619,-0.9519047619047605};

        KRATOS_EXPECT_VECTOR_NEAR(RHS, reference, 1e-13);
    }

    // Checks the LHS of the 3D Wake IncompressiblePotentialFlowElement element.
    KRATOS_TEST_CASE_IN_SUITE(WakeIncompressiblePotentialFlowElementLHS3D,
                              CompressiblePotentialApplicationFastSuite)
    {
        Model this_model;
        ModelPart& model_part = this_model.CreateModelPart("Main", 3);

        GenerateIncompressibleElement3D(model_part);
        Element::Pointer p_element = model_part.pGetElement(1);

        BoundedVector<double,4> distances = PotentialFlowTestUtilities::AssignDistancesToElement<4>();

        p_element->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
        p_element->GetValue(WAKE) = true;

        std::array<double, 8> potential{1.39572,     110.69275, 121.1549827,
                                        104.284736,  2.39572,   46.69275,
                                        100.1549827, 102.284736};
        PotentialFlowTestUtilities::AssignPotentialsToWakeElement<4>(*p_element,distances,potential);

        // Compute LHS
        Matrix LHS = ZeroMatrix(4, 4);

        const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
        p_element->CalculateLeftHandSide(LHS, r_current_process_info);

        std::vector<double> reference{0.263095238095238,-0.185,0.05261904761904762,-0.1307142857142857,-0.263095238095238,
        0.185,-0.05261904761904762,0.1307142857142857,-0.185,0.2356666666666666,-0.03699999999999999,-0.01366666666666666,
        0.185,-0.2356666666666666,0.03699999999999999,0.01366666666666666,0.05261904761904762,-0.03699999999999999,
        0.01052380952380952,-0.02614285714285715,-0.05261904761904762,0.03699999999999999,-0.01052380952380952,
        0.02614285714285715,-0.1114285714285714,-0.01066666666666666,-0.05228571428571429,0.1743809523809524,0,0,0,0,0,
        0,0,0,0.3595238095238095,-0.17,-0.07809523809523812,-0.1114285714285714,0,0,0,0,-0.17,0.238,-0.05733333333333333,
        -0.01066666666666666,0,0,0,0,-0.07809523809523812,-0.05733333333333333,0.1877142857142857,-0.05228571428571429,
        0.1307142857142857,0.01366666666666666,0.02614285714285715,-0.1705238095238095,-0.1307142857142857,
        -0.01366666666666666,-0.02614285714285715,0.1705238095238095};

        for (unsigned int i = 0; i < LHS.size1(); i++) {
            for (unsigned int j = 0; j < LHS.size2(); j++) {
                KRATOS_EXPECT_NEAR(LHS(i, j), reference[i * 8 + j], 1e-13);
            }
        }
    }

    // Checks the RHS of the 3D IncompressiblePotentialFlowElement element.
    KRATOS_TEST_CASE_IN_SUITE(IncompressiblePotentialFlowElementRHS3D,
                              CompressiblePotentialApplicationFastSuite)
    {
        Model this_model;
        ModelPart& model_part = this_model.CreateModelPart("Main", 3);

        GenerateIncompressibleElement3D(model_part);
        Element::Pointer p_element = model_part.pGetElement(1);

        std::array<double, 4> potential{1.39572, 143.39275, 151.1549827, 134.284736};
        PotentialFlowTestUtilities::AssignPotentialsToNormalElement<4>(*p_element,potential);

        // Compute RHS
        Vector RHS = ZeroVector(4);

        const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
        p_element->CalculateRightHandSide(RHS, r_current_process_info);

        std::vector<double> reference{50.64261358895238,-23.79161257453333,-13.02259285120952,-13.82840816320952};

        KRATOS_EXPECT_VECTOR_NEAR(RHS, reference, 1e-13);
    }

    // Checks the LHS of the 3D IncompressiblePotentialFlowElement element.
    KRATOS_TEST_CASE_IN_SUITE(IncompressiblePotentialFlowElementLHS3D,
                              CompressiblePotentialApplicationFastSuite)
    {
        Model this_model;
        ModelPart& model_part = this_model.CreateModelPart("Main", 3);

        GenerateIncompressibleElement3D(model_part);
        Element::Pointer p_element = model_part.pGetElement(1);

        std::array<double, 4> potential{1.39572, 143.39275, 151.1549827, 134.284736};
        PotentialFlowTestUtilities::AssignPotentialsToNormalElement<4>(*p_element,potential);

        // Compute LHS
        Matrix LHS = ZeroMatrix(4, 4);

        const ProcessInfo& r_current_process_info = model_part.GetProcessInfo();
        p_element->CalculateLeftHandSide(LHS, r_current_process_info);

        std::vector<double> reference{0.3595238095238095,-0.17,-0.07809523809523812,-0.1114285714285714,-0.17,0.238,-0.05733333333333333,-0.01066666666666666,-0.07809523809523812,-0.05733333333333333,0.1877142857142857,-0.05228571428571429,-0.1114285714285714,-0.01066666666666666,-0.05228571428571429,0.1743809523809524};

        for (unsigned int i = 0; i < LHS.size1(); i++) {
            for (unsigned int j = 0; j < LHS.size2(); j++) {
                KRATOS_EXPECT_NEAR(LHS(i, j), reference[i * 4 + j], 1e-13);
            }
        }
    }

  } // namespace Testing
}  // namespace Kratos.
