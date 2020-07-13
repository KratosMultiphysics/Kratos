//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Hoang-Giang Bui
//

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "includes/model_part.h"

namespace Kratos {
  namespace Testing {


    KRATOS_TEST_CASE_IN_SUITE(CommunicatorSetAndAddColors, KratosCoreFastSuite)
    {
        Model current_model;

        ModelPart& r_model_part = current_model.CreateModelPart("Main");

        r_model_part.GetCommunicator().SetNumberOfColors(2);

        KRATOS_CHECK_EQUAL(r_model_part.GetCommunicator().GetNumberOfColors(), 2);

        r_model_part.GetCommunicator().AddColors(2);

        KRATOS_CHECK_EQUAL(r_model_part.GetCommunicator().GetNumberOfColors(), 4);
    }

  }  // namespace Testing.
}  // namespace Kratos.
