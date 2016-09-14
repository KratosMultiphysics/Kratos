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

// System includes


// External includes


// Project includes
#include "testing/testing.h"
#include "includes/kratos_parameters.h"

namespace Kratos {
	namespace Testing {

		KRATOS_TEST_CASE_IN_SUITE(ParametersIterator, KratosCoreFastSuite)
		{
			Parameters parameters(R"input(
			{
			   "int_value" : 10,   "double_value": 2.0,   "bool_value" : true,   "string_value" : "hello",
			   "level1":
			   {
			     "list_value":[ 3, "hi", false],
			     "tmp" : 5.0
			   }
			}
			)input");

			auto i_parameter = parameters.begin();
			KRATOS_CHECK_EQUAL(i_parameter.name(), std::string("int_value"));
			KRATOS_CHECK(i_parameter->IsInt());
			KRATOS_CHECK_EQUAL(i_parameter->GetInt(), 10);
			KRATOS_CHECK((++i_parameter)->IsDouble());
			KRATOS_CHECK_EQUAL((i_parameter++)->GetDouble(), 2.0);
			KRATOS_CHECK((*i_parameter).IsBool());
			KRATOS_CHECK_EQUAL((*i_parameter).GetBool(), true);

			unsigned int size = 0;

			for(auto i = parameters.begin() ; i != parameters.end() ; i++)
				size++;

			KRATOS_CHECK_EQUAL(size, 5);


		}
	}
}  // namespace Kratos.
