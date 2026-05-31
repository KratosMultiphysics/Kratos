//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Crescenzio
//
//

// System includes
#include <set>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/model_part.h"

// Application includes
#include "testing/testing.h"
#include "custom_utilities/brute_force_material_point_locator.h"
#include "mpm_application_variables.h"

namespace Kratos::Testing {

    /**
     * Test class BruteForceMaterialPointLocator, public method FindElement
     */
    KRATOS_TEST_CASE_IN_SUITE(MPMBruteForceMaterialPointElementLocator, KratosMPMFastSuite)
    {
        Model model;

        ModelPart& model_part = model.CreateModelPart("GridModelPart");

        Properties::Pointer p_prop = model_part.CreateNewProperties(0);

        model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
        model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
        model_part.CreateNewNode(4, 1.0, 1.0, 0.0);

        model_part.CreateNewElement("Element2D4N", 1, {1, 2, 4, 3}, p_prop);
        auto p_geometry = model_part.GetElement(1).pGetGeometry();

        ModelPart& mpm_model_part = model.CreateModelPart("MPMModelPart");
        const auto& process_info = mpm_model_part.GetProcessInfo();

        for (int i=1; i<6; ++i) {
            for (int j=1; j<6; ++j) {
                double xcoord = 0.2*i - 0.1;
                double ycoord = 0.2*j - 0.1;
                array_1d<double, 3> mp_coord{ xcoord, ycoord, 0.0 };
                auto p_quad = CreateQuadraturePointsUtility<Node>::CreateFromCoordinates(p_geometry, mp_coord, 1.0);
                auto p_elem = mpm_model_part.CreateNewElement("MPMUpdatedLagrangian", 5*(i-1) + j, p_quad, p_prop);
                p_elem->SetValuesOnIntegrationPoints(MP_COORD, { mp_coord }, process_info);
            }
        }

        auto locator = BruteForceMaterialPointLocator(mpm_model_part);

        KRATOS_CHECK_EQUAL(locator.FindElement({0.1, 0.1, 0.0}, 1e-6), 1);
        KRATOS_CHECK_EQUAL(locator.FindElement({0.1, 0.1000001, 0.0}, 1e-6), 1);
        KRATOS_CHECK_EQUAL(locator.FindElement({0.1, 0.15, 0.0}, 1e-6), -1);
        KRATOS_CHECK_EQUAL(locator.FindElement({0.3, 0.901, 0.0}, 1e-2), 10);
        KRATOS_CHECK_EQUAL(locator.FindElement({0.3, 0.901, 0.0}, 1e-4), -1);
        KRATOS_CHECK_EQUAL(locator.FindElement({0.5, 0.501, 0.0}, 0.2), 13);
        KRATOS_CHECK_EQUAL(locator.FindElement({0.9, 0.102, 0.0}, 0.5), 21);
        KRATOS_CHECK_EQUAL(locator.FindElement({0.9, 0.102, 0.0}, 1e-3), -1);

    } // test MPMBruteForceMaterialPointElementLocator

    /**
     * Test class BruteForceMaterialPointLocator, public method FindCondition
     */
    KRATOS_TEST_CASE_IN_SUITE(MPMBruteForceMaterialPointConditionLocator, KratosMPMFastSuite)
    {
        Model model;

        ModelPart& model_part = model.CreateModelPart("GridModelPart");

        Properties::Pointer p_prop = model_part.CreateNewProperties(0);

        model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
        model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
        model_part.CreateNewNode(4, 1.0, 1.0, 0.0);

        model_part.CreateNewElement("Element2D4N", 1, {1, 2, 4, 3}, p_prop);
        auto p_geometry = model_part.GetElement(1).pGetGeometry();

        ModelPart& mpm_model_part = model.CreateModelPart("MPMModelPart");
        const auto& process_info = mpm_model_part.GetProcessInfo();

        for (int i=1; i<6; ++i) {
            for (int j=1; j<6; ++j) {
                double xcoord = 0.2*i - 0.1;
                double ycoord = 0.2*j - 0.1;
                array_1d<double, 3> mp_coord{ xcoord, ycoord, 0.0 };
                auto p_quad = CreateQuadraturePointsUtility<Node>::CreateFromCoordinates(p_geometry, mp_coord, 1.0);
                auto p_cond = mpm_model_part.CreateNewCondition("MPMParticlePenaltyDirichletCondition", 5*(i-1) + j, p_quad, p_prop);
                p_cond->SetValuesOnIntegrationPoints(MP_COORD, { mp_coord }, process_info);
            }
        }

        auto locator = BruteForceMaterialPointLocator(mpm_model_part);

        KRATOS_CHECK_EQUAL(locator.FindCondition({0.1, 0.1, 0.0}, 1e-6), 1);
        KRATOS_CHECK_EQUAL(locator.FindCondition({0.1, 0.1000001, 0.0}, 1e-6), 1);
        KRATOS_CHECK_EQUAL(locator.FindCondition({0.1, 0.15, 0.0}, 1e-6), -1);
        KRATOS_CHECK_EQUAL(locator.FindCondition({0.3, 0.901, 0.0}, 1e-2), 10);
        KRATOS_CHECK_EQUAL(locator.FindCondition({0.3, 0.901, 0.0}, 1e-4), -1);
        KRATOS_CHECK_EQUAL(locator.FindCondition({0.5, 0.501, 0.0}, 0.2), 13);
        KRATOS_CHECK_EQUAL(locator.FindCondition({0.9, 0.102, 0.0}, 0.5), 21);
        KRATOS_CHECK_EQUAL(locator.FindCondition({0.9, 0.102, 0.0}, 1e-3), -1);

    } // test MPMBruteForceMaterialPointConditionLocator

}  // namespace Kratos::Testing
