//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//


// System includes


// External includes


// Project includes
#include "testing/testing.h"
#include "includes/geometrical_object.h"
#include "includes/kratos_flags.h"


namespace Kratos::Testing {

KRATOS_TEST_CASE(GeometricalObject_IsActive)
{
    GeometricalObject geom_obj;

    KRATOS_CHECK(geom_obj.IsActive()); // active by default

    geom_obj.Set(ACTIVE, true);
    KRATOS_CHECK(geom_obj.IsActive());

    geom_obj.Set(ACTIVE, false);
    KRATOS_CHECK_IS_FALSE(geom_obj.IsActive());

    geom_obj.Reset(ACTIVE);
    KRATOS_CHECK(geom_obj.IsActive()); // active by default
}

} // namespace Kratos::Testing
