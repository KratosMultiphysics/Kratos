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
#include "custom_elements/compressible_potential_flow_element.h"
#include "custom_elements/embedded_compressible_potential_flow_element.h"


namespace Kratos {
  namespace Testing {

    typedef ModelPart::IndexType IndexType;
    typedef ModelPart::NodeIterator NodeIteratorType;

    void GenerateCompressibleElement(ModelPart& rModelPart)
    {
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

      // Geometry creation
      rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
      rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
      rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
      std::vector<ModelPart::IndexType> elemNodes{ 1, 2, 3 };
      rModelPart.CreateNewElement("CompressiblePotentialFlowElement2D3N", 1, elemNodes, pElemProp);
    }

    void GenerateCompressibleEmbeddedElement(ModelPart& rModelPart)
    {
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

      // Geometry creation
      rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
      rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0);
      rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0);
      std::vector<ModelPart::IndexType> elemNodes{ 1, 2, 3 };
      rModelPart.CreateNewElement("EmbeddedCompressiblePotentialFlowElement2D3N", 1, elemNodes, pElemProp);
    }

    /** Checks the IncompressiblePotentialFlowElement element.
     * Checks the LHS and RHS computation.
     */
    KRATOS_TEST_CASE_IN_SUITE(CompressiblePotentialFlowElementRHS, CompressiblePotentialApplicationFastSuite)
    {
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateCompressibleElement(model_part);
      Element::Pointer pElement = model_part.pGetElement(1);

      // Define the nodal values
      std::array<double,3> potential;
      potential[0] = 1.0;
      potential[1] = 2.0;
      potential[2] = 3.0;

      for (unsigned int i = 0; i < 3; i++){
        pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = potential[i];
      }
      // Compute RHS and LHS
      Vector RHS = ZeroVector(3);
      Matrix LHS = ZeroMatrix(3, 3);

      pElement->CalculateLocalSystem(LHS, RHS, model_part.GetProcessInfo());

      std::vector<double> reference({0.615561780, 0.0, -0.615561780});

      for (unsigned int i = 0; i < RHS.size(); i++) {
        KRATOS_CHECK_NEAR(RHS(i), reference[i], 1e-6);
      }
    }

    /** Checks the IncompressiblePotentialFlowElement element.
     * Checks the LHS and RHS computation.
     */
    KRATOS_TEST_CASE_IN_SUITE(CompressiblePotentialFlowElementLHS, CompressiblePotentialApplicationFastSuite)
    {
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateCompressibleElement(model_part);
      Element::Pointer pElement = model_part.pGetElement(1);

      // Define the nodal values
      std::array<double,3> potential;
      potential[0] = 1.0;
      potential[1] = 2.0;
      potential[2] = 3.0;

      for (unsigned int i = 0; i < 3; i++){
        pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = potential[i];
      }
      // Compute RHS and LHS
      Vector RHS = ZeroVector(3);
      Matrix LHS = ZeroMatrix(3, 3);

      pElement->CalculateLocalSystem(LHS, RHS, model_part.GetProcessInfo());

      std::array<double,9> reference({0.615556466,-0.615561780,5.314318652e-06,
                                      -0.615561780,1.231123561,-0.615561780,
                                      5.314318652e-06,-0.615561780, 0.615556466});

      for (unsigned int i = 0; i < LHS.size1(); i++) {
        for (unsigned int j = 0; j < LHS.size2(); j++) {
          KRATOS_CHECK_NEAR(LHS(i,j), reference[i*3+j], 1e-6);
        }
      }
    }

    KRATOS_TEST_CASE_IN_SUITE(EmbeddedCompressiblePotentialFlowElementCalculateLocalSystemRHS, CompressiblePotentialApplicationFastSuite)
    {
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateCompressibleEmbeddedElement(model_part);
      Element::Pointer pElement = model_part.pGetElement(1);

      // Define the nodal values
      std::array<double,3> potential({1.0, 2.0, 3.0});
      // Define the distance values
      std::array<double,3> level_set({1.0, -1.0, -1.0});

      for (unsigned int i = 0; i < 3; i++){
        pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = potential[i];
        pElement->GetGeometry()[i].FastGetSolutionStepValue(GEOMETRY_DISTANCE) = level_set[i];
      }

      // Compute RHS and LHS
      Vector RHS = ZeroVector(3);
      Matrix LHS = ZeroMatrix(3, 3);

      pElement->CalculateLocalSystem(LHS, RHS, model_part.GetProcessInfo());

      std::vector<double> reference({0.125625, 0.0, -0.125625});

      KRATOS_CHECK_VECTOR_NEAR(RHS, reference, 1e-6);
    }

    KRATOS_TEST_CASE_IN_SUITE(EmbeddedCompressiblePotentialFlowElementCalculateLocalSystemLHS, CompressiblePotentialApplicationFastSuite)
    {
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateCompressibleEmbeddedElement(model_part);
      Element::Pointer pElement = model_part.pGetElement(1);

      // Define the nodal values
      std::array<double,3> potential({1.0, 2.0, 3.0});
      // Define the distance values
      std::array<double,3> level_set({1.0, -1.0, -1.0});
      for (unsigned int i = 0; i < 3; i++){
        pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = potential[i];
        pElement->GetGeometry()[i].FastGetSolutionStepValue(GEOMETRY_DISTANCE) = level_set[i];
      }

      // Compute RHS and LHS
      Vector RHS = ZeroVector(3);
      Matrix LHS = ZeroMatrix(3, 3);

      pElement->CalculateLocalSystem(LHS, RHS, model_part.GetProcessInfo());

      std::array<double, 9> reference_array({0.251249, -0.25125, 1.08455e-06, -0.25125, 0.502499, -0.25125, 1.08455e-06, -0.25125, 0.251249});
      // Copying to a 3x3 matrix to check against LHS
      Matrix reference(3,3);
      for (unsigned int i = 0; i < reference.size1(); i++) {
        for (unsigned int j = 0; j < reference.size2(); j++) {
          reference(i,j) =  reference_array[i*reference.size1()+j];
        }
      }

      KRATOS_CHECK_MATRIX_NEAR(LHS, reference, 1e-6);
    }

    KRATOS_TEST_CASE_IN_SUITE(CompressiblePotentialFlowElementRHSWake, CompressiblePotentialApplicationFastSuite)
    {
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateCompressibleElement(model_part);
      Element::Pointer pElement = model_part.pGetElement(1);

      // Define the nodal values
      std::array<double,3> potential;
      potential[0] = 1.0;
      potential[1] = 2.0;
      potential[2] = 3.0;

      Vector distances(3);
      distances(0) = 1.0;
      distances(1) = -1.0;
      distances(2) = -1.0;

      pElement->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
      pElement->GetValue(WAKE) = true;

      for (unsigned int i = 0; i < 3; i++){
        if (distances(i) > 0.0){
          pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = potential[i];
        }
        else{
          pElement->GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) = potential[i];
        }
      }
      for (unsigned int i = 0; i < 3; i++){
        if (distances(i) < 0.0){
          pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = potential[i]+5;
        }
        else{
          pElement->GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) = potential[i]+5;
        }
      }

      // Compute RHS and LHS
      Vector RHS = ZeroVector(6);
      Matrix LHS = ZeroMatrix(6, 6);

      pElement->CalculateLocalSystem(LHS, RHS, model_part.GetProcessInfo());

      std::array<double,6> reference({0.615561780, 0.0, 0.0, 0.0, 0.0, -0.615561780});

      for (unsigned int i = 0; i < RHS.size(); i++) {
        KRATOS_CHECK_NEAR(RHS(i), reference[i], 1e-6);
      }
    }

    KRATOS_TEST_CASE_IN_SUITE(CompressiblePotentialFlowElementLHSWake, CompressiblePotentialApplicationFastSuite)
    {
      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateCompressibleElement(model_part);
      Element::Pointer pElement = model_part.pGetElement(1);

      // Define the nodal values
      std::array<double,3> potential;
      potential[0] = 1.0;
      potential[1] = 2.0;
      potential[2] = 3.0;

      Vector distances(3);
      distances(0) = 1.0;
      distances(1) = -1.0;
      distances(2) = -1.0;

      pElement->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
      pElement->GetValue(WAKE) = true;

      for (unsigned int i = 0; i < 3; i++){
        if (distances(i) > 0.0){
          pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = potential[i];
        }
        else{
          pElement->GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) = potential[i];
        }
      }
      for (unsigned int i = 0; i < 3; i++){
        if (distances(i) < 0.0){
          pElement->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = potential[i]+5;
        }
        else{
          pElement->GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) = potential[i]+5;
        }
      }

      // Compute RHS and LHS
      Vector RHS = ZeroVector(6);
      Matrix LHS = ZeroMatrix(6, 6);

      pElement->CalculateLocalSystem(LHS, RHS, model_part.GetProcessInfo());

      // Check the RHS values (the RHS is computed as the LHS x previous_solution,
      // hence, it is assumed that if the RHS is correct, the LHS is correct as well)
      std::array<double,36> reference({0.615556466,-0.615561780,5.314318652e-06,0.0,0.0,0.0,
                                  -0.615561780,1.231123561,-0.615561780,0.615561780,-1.231123561,0.615561780,
                                  5.314318652e-06,-0.615561780, 0.615556466,-5.314318652e-06,0.615561780, -0.615556466,
                                  -0.615556466, 0.615561780,-5.314318652e-06,0.615556466, -0.615561780,5.314318652e-06,
                                  0.0,0.0,0.0,-0.615561780,1.231123561,-0.615561780,
                                  0.0,0.0,0.0,5.314318652e-06,-0.615561780,0.615556466});

      for (unsigned int i = 0; i < LHS.size1(); i++) {
        for (unsigned int j = 0; j < LHS.size2(); j++) {
          KRATOS_CHECK_NEAR(LHS(i,j), reference[6*i+j], 1e-6);
        }
      }
    }

    /** Checks the IncompressiblePotentialFlowElement element.
 * Checks the EquationIdVector.
 */
    KRATOS_TEST_CASE_IN_SUITE(CompressiblePotentialFlowElementEquationIdVector, CompressiblePotentialApplicationFastSuite)
    {

      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateCompressibleElement(model_part);
      Element::Pointer pElement = model_part.pGetElement(1);

      for (unsigned int i = 0; i < 3; i++)
        pElement->GetGeometry()[i].AddDof(VELOCITY_POTENTIAL);

      Element::DofsVectorType ElementalDofList;
      pElement->GetDofList(ElementalDofList, model_part.GetProcessInfo());

      for (int i = 0; i < 3; i++){
        ElementalDofList[i]->SetEquationId(i);
      }

      Element::EquationIdVectorType EquationIdVector;
      pElement->EquationIdVector(EquationIdVector, model_part.GetProcessInfo());

      // Check the EquationIdVector values
      for (unsigned int i = 0; i < EquationIdVector.size(); i++) {
        KRATOS_CHECK(EquationIdVector[i] == i);
      }
    }

    /** Checks the IncompressiblePotentialFlowElement element.
 * Checks the EquationIdVector for the Wake.
 */
    KRATOS_TEST_CASE_IN_SUITE(CompressiblePotentialFlowElementEquationIdVectorWake, CompressiblePotentialApplicationFastSuite)
    {

      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateCompressibleElement(model_part);
      Element::Pointer pElement = model_part.pGetElement(1);
      pElement->SetValue(WAKE, true);

      Vector distances(3);
      distances(0) = -0.5;
      distances(1) = -0.5;
      distances(2) = 0.5;
      pElement->SetValue(WAKE_ELEMENTAL_DISTANCES, distances);

      for (unsigned int i = 0; i < 3; i++) {
        pElement->GetGeometry()[i].AddDof(VELOCITY_POTENTIAL);
        pElement->GetGeometry()[i].AddDof(AUXILIARY_VELOCITY_POTENTIAL);
      }

      Element::DofsVectorType ElementalDofList;
      pElement->GetDofList(ElementalDofList, model_part.GetProcessInfo());

      for (int i = 0; i < 6; i++){
        ElementalDofList[i]->SetEquationId(i);
      }

      Element::EquationIdVectorType EquationIdVector;
      pElement->EquationIdVector(EquationIdVector, model_part.GetProcessInfo());

      //Check the EquationIdVector values
      for (unsigned int i = 0; i < EquationIdVector.size(); i++) {
        KRATOS_CHECK(EquationIdVector[i] == i);
      }
    }
  } // namespace Testing
}  // namespace Kratos.
