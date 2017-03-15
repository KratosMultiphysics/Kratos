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
#include "includes/properties.h"

namespace Kratos {
	namespace Testing {

		KRATOS_TEST_CASE_IN_SUITE(PropertiesHasTable, KratosCoreFastSuite)
		{
			Properties properties(0);
			properties.HasTable(TEMPERATURE, VISCOSITY);
		}
	}
}  // namespace Kratos.
