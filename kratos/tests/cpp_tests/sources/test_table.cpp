//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

// Project includes
#include "testing/testing.h"
#include "includes/table.h"
#include "includes/checks.h"
#include "includes/variables.h"

namespace Kratos {
    namespace Testing {

        KRATOS_TEST_CASE_IN_SUITE(BaseTable, KratosCoreFastSuite)
        {
            Table<double> table;
            for (std::size_t i = 0; i < 6; ++i)
                table.PushBack(static_cast<double>(i), 2.0 * static_cast<double>(i));

            double nearest = (table.GetNearestRow(2.1))[0];
            KRATOS_CHECK_DOUBLE_EQUAL(nearest, 4.0);
            KRATOS_CHECK_DOUBLE_EQUAL(table.GetValue(2.1), 4.2);
            KRATOS_CHECK_DOUBLE_EQUAL(table(2.1), 4.2);
            KRATOS_CHECK_DOUBLE_EQUAL(table(2.0), 4.0);
            KRATOS_CHECK_DOUBLE_EQUAL(table(1.0), 2.0);
            KRATOS_CHECK_DOUBLE_EQUAL(table(7.5), 15.0);
            KRATOS_CHECK_DOUBLE_EQUAL(table(-1.5), -3.0);
            KRATOS_CHECK_DOUBLE_EQUAL(table.GetDerivative(2.1), 2.0);
            KRATOS_CHECK_DOUBLE_EQUAL(table.GetDerivative(1.0), 2.0);
            auto& r_data = table.Data();
            KRATOS_CHECK_EQUAL(r_data.size(), 6);
            // Clear database
            table.Clear();
            KRATOS_CHECK_EQUAL(r_data.size(), 0);

            for (std::size_t i = 0; i < 6; ++i){
                table.PushBack(static_cast<double>(i), 2.0 * static_cast<double>(i));
                table.PushBack(static_cast<double>(i), 2.0 * static_cast<double>(i));
            }
            KRATOS_CHECK_DOUBLE_EQUAL(table(2.1), 4.2);
            KRATOS_CHECK_DOUBLE_EQUAL(table(2.0), 4.0);
            KRATOS_CHECK_DOUBLE_EQUAL(table(5.0), 10.0);
            KRATOS_CHECK_DOUBLE_EQUAL(table(7.5), 10.0);
            KRATOS_CHECK_DOUBLE_EQUAL(table(-1.5), 0.0);
            // Clear database
            table.Clear();
            KRATOS_CHECK_EQUAL(r_data.size(), 0);

            for (std::size_t i = 0; i < 6; ++i){
                table.PushBack(static_cast<double>(i), 2.0 * static_cast<double>(i));
                table.PushBack(static_cast<double>(i), 2.0 * static_cast<double>(i)+1.0);
            }
            KRATOS_CHECK_DOUBLE_EQUAL(table(2.1), 5.1);
            KRATOS_CHECK_DOUBLE_EQUAL(table(2.0), 5.0);
            KRATOS_CHECK_DOUBLE_EQUAL(table(5.0), 11.0);
            KRATOS_CHECK_DOUBLE_EQUAL(table(7.5), 11.0);
            KRATOS_CHECK_DOUBLE_EQUAL(table(-1.5), 0.0);
            // Clear database
            table.Clear();
            KRATOS_CHECK_EQUAL(r_data.size(), 0);

            // Inverse filling with insert
            for (std::size_t i = 6; i > 0; --i)
                table.insert(static_cast<double>(i), 2.0 * static_cast<double>(i));
            KRATOS_CHECK_EQUAL(r_data.size(), 6);
            nearest = (table.GetNearestRow(2.1))[0];
            KRATOS_CHECK_DOUBLE_EQUAL(nearest, 4.0);
            KRATOS_CHECK_DOUBLE_EQUAL(table.GetValue(2.1), 4.2);
            KRATOS_CHECK_DOUBLE_EQUAL(table(2.1), 4.2);
            KRATOS_CHECK_DOUBLE_EQUAL(table.GetDerivative(2.1), 2.0);
        }

    }
}  // namespace Kratos.
