//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//
//

// Project includes
#include "testing/testing.h"
#include "includes/model_part.h"


namespace Kratos {
	namespace Testing {

		KRATOS_TEST_CASE_IN_SUITE(ModelPartSubModelPartsIterator, KratosCoreFastSuite)
		{
			ModelPart model_part("Main");

			model_part.CreateSubModelPart("Inlet");
			model_part.CreateSubModelPart("Inlet2");
			model_part.CreateSubModelPart("Outlet");
			model_part.CreateSubModelPart("Outlet2");

			std::size_t id = 1;
			for(auto i_SubModelPart = model_part.SubModelPartsBegin() ; i_SubModelPart != model_part.SubModelPartsEnd() ; i_SubModelPart++){
				i_SubModelPart->CreateNewNode(id++, 0.00,0.00,0.00);
			}
			KRATOS_CHECK_EQUAL(model_part.NumberOfNodes(), 4);
			KRATOS_CHECK_EQUAL(model_part.GetSubModelPart("Inlet").NumberOfNodes(), 1);
			KRATOS_CHECK_EQUAL(model_part.GetSubModelPart("Outlet").NumberOfNodes(), 1);
			KRATOS_CHECK_EQUAL(model_part.GetSubModelPart("Outlet").GetNode(2).Id(), 2);
			// self.assertEqual(model_part.NumberOfNodes(), 3)
			// self.assertEqual(model_part.GetSubModelPart("Inlet").NumberOfNodes(), 1)
			// self.assertEqual(model_part.GetSubModelPart("Outlet2").NumberOfNodes(), 1)
			// self.assertEqual(model_part.GetSubModelPart("Outlet").NumberOfNodes(), 1)
			// self.assertEqual(model_part.GetSubModelPart("Outlet").GetNode(2).Id, 2)
		}


	}
}  // namespace Kratos.
