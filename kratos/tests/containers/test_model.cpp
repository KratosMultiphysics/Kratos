//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//
//


// Project includes
#include "testing/testing.h"
#include "containers/model.h"

namespace Kratos {
	namespace Testing {

		KRATOS_TEST_CASE_IN_SUITE(ModelGetModelPart, KratosCoreFastSuite)
		{
            auto model_part = Kratos::make_shared<ModelPart>("Main");

			model_part->CreateSubModelPart("Inlet1");

            Model model;

            model.AddModelPart(model_part.get());

            KRATOS_CHECK_EQUAL(model.GetModelPart("Main").Name(), model_part->Name());

            ModelPart& smp = model_part->GetSubModelPart("Inlet1");
            KRATOS_CHECK_EQUAL(model.GetModelPart("Main.Inlet1").Name(), smp.Name());

            KRATOS_CHECK_EXCEPTION_IS_THROWN(model.GetModelPart("Main.Random"),
                "Error: The ModelPart named : \"Random\" was not found as SubModelPart of : \"Main\". The total input string was \"Main.Random\"");

            // TODO this should throw in the future
            // KRATOS_CHECK_EXCEPTION_IS_THROWN(model.GetModelPart("Inlet1"),
            //     "Error: The ModelPart named : \"Inlet1\" was not found as root-ModelPart. The total input string was \"Inlet1\"");

            KRATOS_CHECK_EXCEPTION_IS_THROWN(model.GetModelPart("Maiiiiin"),
                "Error: The ModelPart named : \"Maiiiiin\" was not found as root-ModelPart. The total input string was \"Maiiiiin\"");
		}

		KRATOS_TEST_CASE_IN_SUITE(ModelHasModelPart, KratosCoreFastSuite)
		{
            auto model_part = Kratos::make_shared<ModelPart>("Main");

			model_part->CreateSubModelPart("Inlet1");

            Model model;

            model.AddModelPart(model_part.get());

            KRATOS_CHECK(model.HasModelPart("Main"));
            KRATOS_CHECK(model.HasModelPart("Main.Inlet1"));

            KRATOS_CHECK_IS_FALSE(model.HasModelPart("Inlet1"));

            ModelPart& smp = model_part->GetSubModelPart("Inlet1");
            smp.CreateSubModelPart("SubInlet");

            KRATOS_CHECK(model.HasModelPart("Main.Inlet1.SubInlet"));
            KRATOS_CHECK_IS_FALSE(model.HasModelPart("Main.SubInlet"));

            KRATOS_CHECK_IS_FALSE(model.HasModelPart("Random"));
		}

	}   // namespace Testing
}  // namespace Kratos.


