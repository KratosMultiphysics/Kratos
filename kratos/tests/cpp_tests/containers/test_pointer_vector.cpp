//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Daniel Diez
//
//

// System includes

// External includes

// Project includes
#include "includes/element.h"
#include "containers/pointer_vector.h"
#include "testing/testing.h"

namespace Kratos::Testing {

    KRATOS_TEST_CASE_IN_SUITE(PointerVectorPushBack, KratosCoreFastSuite)
    {
        PointerVector<const Element> test_container;
        auto p_element_1 = Kratos::make_intrusive<Element>(1);
        auto p_element_2 = Kratos::make_intrusive<Element>(2);
        auto p_element_3 = Kratos::make_intrusive<Element>(3);
        test_container.push_back(p_element_1);
        test_container.push_back(p_element_2);
        test_container.push_back(p_element_3);

        KRATOS_EXPECT_EQ(test_container.begin()->Id(), 1);
        KRATOS_EXPECT_EQ((test_container.begin()+1)->Id(), 2);
        KRATOS_EXPECT_EQ((test_container.end()-1)->Id(), 3);
    }
    KRATOS_TEST_CASE_IN_SUITE(PointerVectorEmplaceBack, KratosCoreFastSuite)
    {
        PointerVector<const Element> test_container;
        test_container.emplace_back(Kratos::make_intrusive<Element>(1));
        test_container.emplace_back(Kratos::make_intrusive<Element>(2));
        test_container.emplace_back(Kratos::make_intrusive<Element>(3));

        KRATOS_EXPECT_EQ(test_container.begin()->Id(), 1);
        KRATOS_EXPECT_EQ((test_container.begin()+1)->Id(), 2);
        KRATOS_EXPECT_EQ((test_container.end()-1)->Id(), 3);
    }

} // namespace Kratos.
