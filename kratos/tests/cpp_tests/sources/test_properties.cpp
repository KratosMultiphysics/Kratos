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
#include "includes/variables.h"

namespace Kratos {
	namespace Testing {

		KRATOS_TEST_CASE_IN_SUITE(PropertiesAssignment, KratosCoreFastSuite)
		{
			Properties property_1(1);
			Properties property_2(2);

			Table<double> table;
			property_1.SetTable(TEMPERATURE, VISCOSITY, table);

			property_1.SetValue(TEMPERATURE, 1.0);
			property_1.SetValue(DENSITY, 2.0);

			KRATOS_CHECK(property_2.IsEmpty());

			property_2 = property_1;

			KRATOS_CHECK_EQUAL(property_2.GetValue(TEMPERATURE), 1.0);
			KRATOS_CHECK_EQUAL(property_2.GetValue(DENSITY), 2.0);
			KRATOS_CHECK(property_2.HasTables());
			KRATOS_CHECK(property_2.HasTable(TEMPERATURE, VISCOSITY));
			KRATOS_CHECK_IS_FALSE(property_2.HasTable(TEMPERATURE, DISPLACEMENT_X));
			KRATOS_CHECK_IS_FALSE(property_2.HasTable(VISCOSITY, TEMPERATURE));
		}

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

		KRATOS_TEST_CASE_IN_SUITE(PropertiesHasVariables, KratosCoreFastSuite)
		{
			Properties property(0);
			KRATOS_CHECK_IS_FALSE(property.HasVariables());

			property.SetValue(TEMPERATURE, 1.0);

			KRATOS_CHECK(property.HasVariables());
		}

		KRATOS_TEST_CASE_IN_SUITE(PropertiesHasTables, KratosCoreFastSuite)
		{
			Properties property(0);
			KRATOS_CHECK_IS_FALSE(property.HasTables());

			Table<double> table;
			property.SetTable(TEMPERATURE, VISCOSITY, table);

			KRATOS_CHECK(property.HasTables());
		}

		KRATOS_TEST_CASE_IN_SUITE(PropertiesIsEmpty, KratosCoreFastSuite)
		{
			Properties property(0);
			KRATOS_CHECK(property.IsEmpty());

			Table<double> table;
			property.SetTable(TEMPERATURE, VISCOSITY, table);

			KRATOS_CHECK_IS_FALSE(property.IsEmpty());

			Properties property1(1);
			KRATOS_CHECK(property1.IsEmpty());

			property1.SetValue(TEMPERATURE, 1.0);

			KRATOS_CHECK_IS_FALSE(property1.IsEmpty());

			Properties property2(2);
			KRATOS_CHECK(property2.IsEmpty());

			property2.SetTable(TEMPERATURE, VISCOSITY, table);

			property2.SetValue(TEMPERATURE, 1.0);

			KRATOS_CHECK_IS_FALSE(property2.IsEmpty());
		}

	}
}  // namespace Kratos.
