//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// System includes

// External includes

// Project includes
#include "includes/element.h"
#include "containers/model.h"
#include "containers/pointer_vector_set.h"
#include "testing/testing.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(PointerVectorSetCBeginAndCEnd, KratosCoreFastSuite)
{
    PointerVectorSet<const Element> test_container;
    auto p_element_1 = Kratos::make_intrusive<Element>(1);
    auto p_element_2 = Kratos::make_intrusive<Element>(2);
    test_container.push_back(p_element_1);
    test_container.push_back(p_element_2);

    KRATOS_CHECK_EQUAL(test_container.cbegin()->Id(), 1);
    KRATOS_CHECK_EQUAL((test_container.cend()-1)->Id(), 2);
}

} // namespace Testing.
} // namespace Kratos.
