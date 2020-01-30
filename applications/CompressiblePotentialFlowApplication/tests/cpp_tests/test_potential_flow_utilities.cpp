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


  } // namespace Testing
}  // namespace Kratos.
