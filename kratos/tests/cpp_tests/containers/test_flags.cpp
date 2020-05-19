//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//
//

// Project includes
#include "containers/flags.h"
#include "includes/kratos_flags.h"
#include "testing/testing.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(KratosFlagsSetResetClear, KratosCoreFastSuite) {

    Kratos::Flags flags;

    // Definition and setting
    KRATOS_CHECK_EQUAL(flags.IsDefined(INLET), false);
    KRATOS_CHECK_EQUAL(flags.IsDefined(OUTLET), false);
    KRATOS_CHECK_EQUAL(flags.IsDefined(INTERFACE), false);
    KRATOS_CHECK_EQUAL(flags.Is(INLET), false);
    KRATOS_CHECK_EQUAL(flags.Is(OUTLET), false);
    KRATOS_CHECK_EQUAL(flags.Is(INTERFACE), false);

    flags.Set(INLET); // default set (sets to value bit of argument, which is always true for named flags)

    KRATOS_CHECK_EQUAL(flags.IsDefined(INLET), true);
    KRATOS_CHECK_EQUAL(flags.Is(INLET), true);
    // at least once, check that the rest of the flag remains untouched
    KRATOS_CHECK_EQUAL(flags.IsDefined(OUTLET), false);
    KRATOS_CHECK_EQUAL(flags.IsDefined(INTERFACE), false);
    KRATOS_CHECK_EQUAL(flags.Is(OUTLET), false);
    KRATOS_CHECK_EQUAL(flags.Is(INTERFACE), false);

    flags.Set(OUTLET, false); // explicit set

    KRATOS_CHECK_EQUAL(flags.IsDefined(OUTLET), true);
    KRATOS_CHECK_EQUAL(flags.Is(OUTLET), false);

    flags.Set(INTERFACE, true); // explicit set

    KRATOS_CHECK_EQUAL(flags.IsDefined(INTERFACE), true);
    KRATOS_CHECK_EQUAL(flags.Is(INTERFACE), true);

    flags.Reset(INLET);

    // Check that the Reset worked only on INLET
    KRATOS_CHECK_EQUAL(flags.IsDefined(INLET), false);
    KRATOS_CHECK_EQUAL(flags.IsDefined(OUTLET), true);
    KRATOS_CHECK_EQUAL(flags.IsDefined(INTERFACE), true);
    KRATOS_CHECK_EQUAL(flags.Is(INLET), false);
    KRATOS_CHECK_EQUAL(flags.Is(OUTLET), false);
    KRATOS_CHECK_EQUAL(flags.Is(INTERFACE), true);

    flags.Clear();

    // Clear unsets everything
    KRATOS_CHECK_EQUAL(flags.IsDefined(INLET), false);
    KRATOS_CHECK_EQUAL(flags.IsDefined(OUTLET), false);
    KRATOS_CHECK_EQUAL(flags.IsDefined(INTERFACE), false);
    KRATOS_CHECK_EQUAL(flags.Is(INLET), false);
    KRATOS_CHECK_EQUAL(flags.Is(OUTLET), false);
    KRATOS_CHECK_EQUAL(flags.Is(INTERFACE), false);
}

KRATOS_TEST_CASE_IN_SUITE(KratosFlagsSetMultiple, KratosCoreFastSuite) {

    Kratos::Flags flags;

    KRATOS_CHECK_EQUAL(flags.IsDefined(INLET), false);
    KRATOS_CHECK_EQUAL(flags.IsDefined(OUTLET), false);
    KRATOS_CHECK_EQUAL(flags.IsDefined(INTERFACE), false);
    KRATOS_CHECK_EQUAL(flags.IsDefined(VISITED), false);
    KRATOS_CHECK_EQUAL(flags.Is(INLET), false);
    KRATOS_CHECK_EQUAL(flags.Is(OUTLET), false);
    KRATOS_CHECK_EQUAL(flags.Is(INTERFACE), false);
    KRATOS_CHECK_EQUAL(flags.Is(VISITED), false);

    // set by named flag uses both the 'defined' and 'value' bits of the argument
    flags.Set(INLET | OUTLET); // multiple set with 'or'

    KRATOS_CHECK_EQUAL(flags.IsDefined(INLET), true);
    KRATOS_CHECK_EQUAL(flags.IsDefined(OUTLET), true);
    KRATOS_CHECK_EQUAL(flags.Is(INLET), true);
    KRATOS_CHECK_EQUAL(flags.Is(OUTLET), true);

    flags.Set(INTERFACE & VISITED); // 'and' of named flags also defines both of them (but the value bit is false, so they get set to false)

    KRATOS_CHECK_EQUAL(flags.IsDefined(INTERFACE), true);
    KRATOS_CHECK_EQUAL(flags.IsDefined(VISITED), true);
    KRATOS_CHECK_EQUAL(flags.Is(INTERFACE), false);
    KRATOS_CHECK_EQUAL(flags.Is(VISITED), false);

    flags.Clear();

    // set by flag and bool uses 'defined' bits of the flag and the sets the 'value' bits according to the bool argument
    flags.Set(INLET | OUTLET, false); // sets both named flags to false

    KRATOS_CHECK_EQUAL(flags.IsDefined(INLET), true);
    KRATOS_CHECK_EQUAL(flags.IsDefined(OUTLET), true);
    KRATOS_CHECK_EQUAL(flags.Is(INLET), false);
    KRATOS_CHECK_EQUAL(flags.Is(OUTLET), false);

    flags.Set(INTERFACE & VISITED, true); // sets both named flags to true

    KRATOS_CHECK_EQUAL(flags.IsDefined(INTERFACE), true);
    KRATOS_CHECK_EQUAL(flags.IsDefined(VISITED), true);
    KRATOS_CHECK_EQUAL(flags.Is(INTERFACE), true);
    KRATOS_CHECK_EQUAL(flags.Is(VISITED), true);
}

KRATOS_TEST_CASE_IN_SUITE(KratosFlagsEquality, KratosCoreFastSuite) {

    KRATOS_CHECK_EQUAL( INLET == INLET.AsFalse(), false );
    KRATOS_CHECK_EQUAL( INLET != INLET.AsFalse(), true );
    // trivial check (it is acutally 0 == 0) just to see that everything works when the first argument of == is not a lhs
    KRATOS_CHECK_EQUAL( (INLET.AsFalse() & INLET) == INLET.AsFalse(), true );
}

KRATOS_TEST_CASE_IN_SUITE(KratosFlagsOperators, KratosCoreFastSuite) {

    const Kratos::Flags flags1(INLET);
    const Kratos::Flags flags2(INLET.AsFalse() | OUTLET);

    Kratos::Flags flags_or = flags1 | flags2;

    KRATOS_CHECK_EQUAL(flags_or.IsDefined(INLET), true);
    KRATOS_CHECK_EQUAL(flags_or.IsDefined(OUTLET), true);
    KRATOS_CHECK_EQUAL(flags_or.Is(INLET), true);
    KRATOS_CHECK_EQUAL(flags_or.Is(OUTLET), true);

    flags_or.Clear();
    flags_or = (flags1.AsFalse()) | flags2;

    KRATOS_CHECK_EQUAL(flags_or.IsDefined(INLET), true);
    KRATOS_CHECK_EQUAL(flags_or.IsDefined(OUTLET), true);
    KRATOS_CHECK_EQUAL(flags_or.Is(INLET), false);
    KRATOS_CHECK_EQUAL(flags_or.Is(OUTLET), true);


    Kratos::Flags flags_and = flags1 & flags2;

    KRATOS_CHECK_EQUAL(flags_and.IsDefined(INLET), true);
    KRATOS_CHECK_EQUAL(flags_and.IsDefined(OUTLET), true);
    KRATOS_CHECK_EQUAL(flags_and.Is(INLET), false);
    KRATOS_CHECK_EQUAL(flags_and.Is(OUTLET), false);

    flags_and.Clear();
    flags_and = (~flags2) & flags1;

    KRATOS_CHECK_EQUAL(flags_and.IsDefined(INLET), true);
    KRATOS_CHECK_EQUAL(flags_and.IsDefined(OUTLET), true);
    KRATOS_CHECK_EQUAL(flags_and.Is(INLET), true);
    KRATOS_CHECK_EQUAL(flags_and.Is(OUTLET), false);

    Kratos::Flags flags3(OUTLET.AsFalse() | VISITED);

    flags3 |= flags2;

    KRATOS_CHECK_EQUAL(flags3.IsDefined(INLET), true);
    KRATOS_CHECK_EQUAL(flags3.IsDefined(OUTLET), true);
    KRATOS_CHECK_EQUAL(flags3.IsDefined(VISITED), true);
    KRATOS_CHECK_EQUAL(flags3.Is(INLET), false);
    KRATOS_CHECK_EQUAL(flags3.Is(OUTLET), true);
    KRATOS_CHECK_EQUAL(flags3.Is(VISITED), true);

    Kratos::Flags flags4(OUTLET | VISITED);

    flags4 &= flags2;

    KRATOS_CHECK_EQUAL(flags4.IsDefined(INLET), true);
    KRATOS_CHECK_EQUAL(flags4.IsDefined(OUTLET), true);
    KRATOS_CHECK_EQUAL(flags4.IsDefined(VISITED), true);
    KRATOS_CHECK_EQUAL(flags4.Is(INLET), false);
    KRATOS_CHECK_EQUAL(flags4.Is(OUTLET), true);
    KRATOS_CHECK_EQUAL(flags4.Is(VISITED), false);
}

}
}