//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Marc Núñez
//
//

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "custom_processes/apply_embedded_flags_process.h"

namespace Kratos {
  namespace Testing {

    typedef ModelPart::IndexType IndexType;
    typedef ModelPart::NodeIterator NodeIteratorType;

    void GenerateElements(ModelPart& rModelPart)
    {
      // Variables addition
      rModelPart.AddNodalSolutionStepVariable(GEOMETRY_DISTANCE);

      // Set the element properties
      rModelPart.CreateNewProperties(0);
      Properties::Pointer pElemProp = rModelPart.pGetProperties(0);

      // Geometry creation
      rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0);
      rModelPart.CreateNewNode(2, 1.0, 1.0, 0.0);
      rModelPart.CreateNewNode(3, 2.0, 0.0, 0.0);
      rModelPart.CreateNewNode(4, 3.0, 1.0, 0.0);
      rModelPart.CreateNewNode(5, 4.0, 0.0, 0.0);
      rModelPart.CreateNewNode(6, 5.0, 1.0, 0.0);

      std::vector<ModelPart::IndexType> elemNodes1{ 1, 3, 2 };
      std::vector<ModelPart::IndexType> elemNodes2{ 2, 3, 4 };
      std::vector<ModelPart::IndexType> elemNodes3{ 3, 5, 4 };
      std::vector<ModelPart::IndexType> elemNodes4{ 5, 6, 4 };

      rModelPart.CreateNewElement("Element2D3N", 1, elemNodes1, pElemProp);
      rModelPart.CreateNewElement("Element2D3N", 2, elemNodes2, pElemProp);
      rModelPart.CreateNewElement("Element2D3N", 3, elemNodes3, pElemProp);
      rModelPart.CreateNewElement("Element2D3N", 4, elemNodes4, pElemProp);

    }

    KRATOS_TEST_CASE_IN_SUITE(ApplyEmbeddedFlags, CompressiblePotentialApplicationFastSuite)
    {
      /* This test computes the distance to a vertical line located at x=2.5.
      Then, it checks which elements are cut by this vertical line and which
      elements are deactivated*/

      Model this_model;
      ModelPart& model_part = this_model.CreateModelPart("Main", 3);

      GenerateElements(model_part);

      double vertical_x = 2.5;

      for(auto it=model_part.NodesBegin(); it!=model_part.NodesEnd(); ++it)
      {
        double x_value = it->X();
        it->FastGetSolutionStepValue(GEOMETRY_DISTANCE)= vertical_x-x_value;
      }
      ApplyEmbeddedFlagsProcess apply_embedded_flags_processes(model_part);
      apply_embedded_flags_processes.Execute();

      KRATOS_CHECK(model_part.GetElement(1).Is(ACTIVE) && model_part.GetElement(1).IsNot(BOUNDARY));
      KRATOS_CHECK(model_part.GetElement(2).Is(ACTIVE) && model_part.GetElement(2).Is(BOUNDARY));
      KRATOS_CHECK(model_part.GetElement(3).Is(ACTIVE) && model_part.GetElement(3).Is(BOUNDARY));
      KRATOS_CHECK(model_part.GetElement(4).IsNot(ACTIVE) && model_part.GetElement(4).IsNot(BOUNDARY));
    }
  } // namespace Testing
}  // namespace Kratos.
