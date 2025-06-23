//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// Project includes
#include "testing/testing.h"
#include "includes/table.h"
#include "includes/expect.h"
#include "includes/variables.h"

namespace Kratos::Testing {

KRATOS_TEST_CASE_IN_SUITE(BaseTable, KratosCoreFastSuite)
{
    Table<double> table;
    for (std::size_t i = 0; i < 6; ++i)
        table.PushBack(static_cast<double>(i), 2.0 * static_cast<double>(i));
            
    double nearest = (table.GetNearestRow(2.1))[0];
    KRATOS_EXPECT_DOUBLE_EQ(nearest, 4.0);
    KRATOS_EXPECT_DOUBLE_EQ(table.GetValue(2.1), 4.2);
    KRATOS_EXPECT_DOUBLE_EQ(table(2.1), 4.2);
    KRATOS_EXPECT_DOUBLE_EQ(table.GetDerivative(2.1), 2.0);

    auto& r_data = table.Data();
    KRATOS_EXPECT_EQ(r_data.size(), 6);
            
    // Clear database
    table.Clear();
    KRATOS_EXPECT_EQ(r_data.size(), 0);

    // Inverse filling with insert
    for (std::size_t i = 6; i > 0; --i)
        table.insert(static_cast<double>(i), 2.0 * static_cast<double>(i));

    KRATOS_EXPECT_EQ(r_data.size(), 6);

    nearest = (table.GetNearestRow(2.1))[0];
    KRATOS_EXPECT_DOUBLE_EQ(nearest, 4.0);
    KRATOS_EXPECT_DOUBLE_EQ(table.GetValue(2.1), 4.2);
    KRATOS_EXPECT_DOUBLE_EQ(table(2.1), 4.2);
    KRATOS_EXPECT_DOUBLE_EQ(table.GetDerivative(2.1), 2.0);
}

KRATOS_TEST_CASE_IN_SUITE(NamesOfXAndYInTable, KratosCoreFastSuite)
{
    Table<double> table; // uses the class template

    // New tables shouldn't have any names set
    KRATOS_EXPECT_TRUE(table.NameOfX().empty());
    KRATOS_EXPECT_TRUE(table.NameOfY().empty());

    table.SetNameOfX("Foo");
    KRATOS_EXPECT_EQ(table.NameOfX(), "Foo");

    table.SetNameOfY("Bar");
    KRATOS_EXPECT_EQ(table.NameOfY(), "Bar");
}

KRATOS_TEST_CASE_IN_SUITE(NamesOfXAndYInTableSpecialization, KratosCoreFastSuite)
{
    Table<double, double> table; // uses the template specialization

    // New tables shouldn't have any names set
    KRATOS_EXPECT_TRUE(table.NameOfX().empty());
    KRATOS_EXPECT_TRUE(table.NameOfY().empty());

    table.SetNameOfX("Foo");
    KRATOS_EXPECT_EQ(table.NameOfX(), "Foo");

    table.SetNameOfY("Bar");
    KRATOS_EXPECT_EQ(table.NameOfY(), "Bar");
}

KRATOS_TEST_CASE_IN_SUITE(TableDerivativeUsesExtrapolationWhenOutsideOfDomain, KratosCoreFastSuite)
{
    Table<double, double> table;
    table.PushBack(0.0, 0.0);
    table.PushBack(1.0, 2.0);
    table.PushBack(3.0, 3.0);

    constexpr auto abs_tolerance = 1.0e-08;
    KRATOS_EXPECT_NEAR(table.GetDerivative(-2.0), 2.0, abs_tolerance); // before first point in table
    KRATOS_EXPECT_NEAR(table.GetDerivative(0.0), 2.0, abs_tolerance); // at first point in table
    KRATOS_EXPECT_NEAR(table.GetDerivative(3.0), 0.5, abs_tolerance); // at last point in table
    KRATOS_EXPECT_NEAR(table.GetDerivative(5.0), 0.5, abs_tolerance); // beyond last point in table
}

}  // namespace Kratos::Testing.
