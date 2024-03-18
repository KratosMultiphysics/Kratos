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

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "includes/model_part.h"

// Utilities
#include "utilities/model_part_utils.h"
#include "utilities/cpp_tests_utilities.h"

namespace Kratos::Testing {

KRATOS_TEST_CASE_IN_SUITE(ModelPartUtilsFromConnectivityGenerateElementsSimple, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");
    ModelPart& r_copy_model_part = current_model.CreateModelPart("Copy");
    CppTestsUtilities::Create2DGeometry(r_model_part);
    const std::size_t number_of_elements = r_model_part.NumberOfElements();
    std::vector<std::vector<std::size_t>> connectivities(number_of_elements);
    auto it_elem_begin = r_model_part.ElementsBegin();
    for (std::size_t i = 0; i < number_of_elements; ++i) {
        auto it_elem = it_elem_begin + i;
        const auto& r_geom = it_elem->GetGeometry();
        connectivities[i].resize(3);
        connectivities[i][0] = r_geom[0].Id();
        connectivities[i][1] = r_geom[1].Id();
        connectivities[i][2] = r_geom[2].Id();
    }

    ModelPartUtils::GenerateEntitiesFromConnectivities<Element>(
        "Element2D3N",
        connectivities,
        r_model_part.Nodes(),
        r_copy_model_part.Elements(),
        nullptr
    );

    // Sort
    r_copy_model_part.Elements().Sort();

    // Initial checks
    KRATOS_EXPECT_EQ(0, r_copy_model_part.NumberOfNodes());
    KRATOS_EXPECT_EQ(number_of_elements, r_copy_model_part.NumberOfElements());

    // Connectivity checks
    auto it_elem_copied_begin = r_copy_model_part.ElementsBegin();
    for (std::size_t i = 0; i < number_of_elements; ++i) {
        auto it_elem = it_elem_copied_begin + i;
        KRATOS_EXPECT_EQ(i + 1, it_elem->Id());
        const auto& r_geom = it_elem->GetGeometry();
        KRATOS_EXPECT_EQ(connectivities[i][0], r_geom[0].Id());
        KRATOS_EXPECT_EQ(connectivities[i][1], r_geom[1].Id());
        KRATOS_EXPECT_EQ(connectivities[i][2], r_geom[2].Id());
    }
}

KRATOS_TEST_CASE_IN_SUITE(ModelPartUtilsFromConnectivityGenerateElements, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");
    ModelPart& r_copy_model_part = current_model.CreateModelPart("Copy");
    CppTestsUtilities::Create2DGeometry(r_model_part);
    const std::size_t number_of_elements = r_model_part.NumberOfElements();
    std::vector<std::size_t> properties_ids(number_of_elements);
    std::vector<std::size_t> elements_ids(number_of_elements);
    std::vector<std::vector<std::size_t>> connectivities(number_of_elements);
    auto it_elem_begin = r_model_part.ElementsBegin();
    for (std::size_t i = 0; i < number_of_elements; ++i) {
        auto it_elem = it_elem_begin + i;
        elements_ids[i] = it_elem->Id();
        properties_ids[i] = it_elem->GetProperties().Id();
        const auto& r_geom = it_elem->GetGeometry();
        connectivities[i].resize(3);
        connectivities[i][0] = r_geom[0].Id();
        connectivities[i][1] = r_geom[1].Id();
        connectivities[i][2] = r_geom[2].Id();
    }

    ModelPartUtils::GenerateEntitiesFromConnectivities<Element>(
        "Element2D3N",
        elements_ids,
        properties_ids,
        connectivities,
        r_model_part.Nodes(),
        r_model_part.rProperties(),
        r_copy_model_part.Elements()
    );

    // Sort
    r_copy_model_part.Elements().Sort();

    // Initial checks
    KRATOS_EXPECT_EQ(0, r_copy_model_part.NumberOfNodes());
    KRATOS_EXPECT_EQ(number_of_elements, r_copy_model_part.NumberOfElements());

    // Connectivity checks
    auto it_elem_copied_begin = r_copy_model_part.ElementsBegin();
    for (std::size_t i = 0; i < number_of_elements; ++i) {
        auto it_elem = it_elem_copied_begin + i;
        KRATOS_EXPECT_EQ(elements_ids[i], it_elem->Id());
        KRATOS_EXPECT_EQ(properties_ids[i], it_elem->GetProperties().Id());
        const auto& r_geom = it_elem->GetGeometry();
        KRATOS_EXPECT_EQ(connectivities[i][0], r_geom[0].Id());
        KRATOS_EXPECT_EQ(connectivities[i][1], r_geom[1].Id());
        KRATOS_EXPECT_EQ(connectivities[i][2], r_geom[2].Id());
    }
}

} // namespace Kratos::Testing.
