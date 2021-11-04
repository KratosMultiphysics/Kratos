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

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/geometry_data.h"

namespace Kratos {
namespace Testing {

/// Test number of integration methods
KRATOS_TEST_CASE_IN_SUITE(GeometryDataNumberOfIntegrationMethods, KratosCoreGeometriesFastSuite)
{
    const std::size_t n_integration_methods = GeometryData::IntegrationMethod::NumberOfIntegrationMethods;
    // auto static_cast<GeometryData::IntegrationMethod>(10)

    KRATOS_CHECK_EQUAL(n_integration_methods, 10);
}


} // namespace Testing.
} // namespace Kratos.
