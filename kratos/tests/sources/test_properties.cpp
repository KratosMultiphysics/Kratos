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
#include "includes/checks.h"

namespace Kratos {
	namespace Testing {

		KRATOS_TEST_CASE_IN_SUITE(PropertiesHasTable, KratosCoreFastSuite)
		{
			Properties properties(0);
			KRATOS_CHECK_IS_FALSE(properties.HasTable(TEMPERATURE, VISCOSITY));

			Table<double> table;
			properties.SetTable(TEMPERATURE, VISCOSITY, table);

			KRATOS_CHECK(properties.HasTable(TEMPERATURE, VISCOSITY));
			KRATOS_CHECK_IS_FALSE(properties.HasTable(TEMPERATURE, DISPLACEMENT_X));
			KRATOS_CHECK_IS_FALSE(properties.HasTable(VISCOSITY, TEMPERATURE));
		}
	}
}  // namespace Kratos.
