//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license:
// kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//
//

// System includes
#include <unordered_set>

// External includes

// Project includes
#include "containers/data_value_container.h"
#include "includes/variables.h"
#include "testing/testing.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(DataValueContainerHas, KratosCoreFastSuite) {
    DataValueContainer container;
    double area = 0.0;
    container.SetValue(NODAL_AREA, area);
    const bool check_true = container.Has(NODAL_AREA);
    const bool check_false = container.Has(DISPLACEMENT);

    KRATOS_CHECK_EQUAL(check_true, true);
    KRATOS_CHECK_EQUAL(check_false, false);
}

KRATOS_TEST_CASE_IN_SUITE(DataValueContainer, KratosCoreFastSuite) {
    DataValueContainer container;
    Vector original_distances(4);
    original_distances[0] = 0.00;
    original_distances[1] = 0.10;
    original_distances[2] = 0.20;
    original_distances[3] = 0.30;
    container.SetValue(ELEMENTAL_DISTANCES, original_distances);
    auto& distances = container.GetValue(ELEMENTAL_DISTANCES);

    for (std::size_t i = 0; i < distances.size(); i++)
        KRATOS_CHECK_EQUAL(distances[i], original_distances[i]);
}

KRATOS_TEST_CASE_IN_SUITE(DataValueContainerMerge, KratosCoreFastSuite) {
    DataValueContainer container_origin;
    DataValueContainer container_target;

    const double density = 1000.0;
    const double viscosity_1 = 1e-3;
    const double viscosity_2 = 2e-3;

    container_origin.SetValue(DENSITY, density);
    container_origin.SetValue(VISCOSITY, viscosity_1);
    container_target.SetValue(VISCOSITY, viscosity_2);
    Flags options;
    //First case: do not overwrite old values
    options.Set(DataValueContainer::OVERWRITE_OLD_VALUES, false);
    container_target.Merge(container_origin, options);

    KRATOS_CHECK_EQUAL(container_target.GetValue(DENSITY), density);
    KRATOS_CHECK_EQUAL(container_target.GetValue(VISCOSITY), viscosity_2);

    //Second case: do overwrite old values
    options.Set(DataValueContainer::OVERWRITE_OLD_VALUES, true);
    container_target.Merge(container_origin, options);
    KRATOS_CHECK_EQUAL(container_target.GetValue(DENSITY), density);
    KRATOS_CHECK_EQUAL(container_target.GetValue(VISCOSITY), viscosity_1);

}
}
} // namespace Kratos.
