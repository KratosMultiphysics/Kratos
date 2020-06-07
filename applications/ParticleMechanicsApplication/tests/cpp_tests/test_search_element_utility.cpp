//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//


// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "particle_mechanics_application_variables.h"
#include "containers/model.h"

#include "custom_utilities/mpm_search_element_utility.h"
#include "utilities/quadrature_points_utility.h"

namespace Kratos
{
namespace Testing
{

    // Tolerance
    static constexpr double tolerance = 1.0e-6;

    void PrepareModelPart(
        ModelPart& rModelPart,
        ModelPart& rBackgroundModelPart)
    {
        // Coordinates of Material Point
        array_1d<double, 3> mp_coordinate;
        mp_coordinate[0] = 0.0;
        mp_coordinate[1] = 0.2;
        mp_coordinate[2] = 0.0;

        // Properties
        Properties::Pointer p_elem_prop = rModelPart.CreateNewProperties(0);

        // Elements
        auto p_quad = CreateQuadraturePointsUtility<Node<3>>::CreateFromCoordinates(
            rBackgroundModelPart.GetElement(1).pGetGeometry(), mp_coordinate, 1.5);

        const Element& new_element = KratosComponents<Element>::Get("UpdatedLagrangian2D4N");
        Element::Pointer p_element = new_element.Create(
            2, p_quad, p_elem_prop);

        rModelPart.AddElement(p_element);
    }

    void PrepareBackgroundModelPart(ModelPart& rBackgroundModelPart)
    {
        // Nodes
        auto p_node_1 = rBackgroundModelPart.CreateNewNode(1, -0.5, -0.5, 0.0);
        auto p_node_2 = rBackgroundModelPart.CreateNewNode(2,  0.5, -0.5, 0.0);
        auto p_node_3 = rBackgroundModelPart.CreateNewNode(3,  0.5,  0.5, 0.0);
        auto p_node_4 = rBackgroundModelPart.CreateNewNode(4, -0.5,  0.5, 0.0);

        auto p_node_9 = rBackgroundModelPart.CreateNewNode(9, 1.5, -0.5, 0.0);
        auto p_node_10 = rBackgroundModelPart.CreateNewNode(10, 1.5, 0.5, 0.0);

        rBackgroundModelPart.CreateNewElement(
            "Element2D4N", 1, { 1, 2, 3, 4 }, nullptr);
        rBackgroundModelPart.CreateNewElement(
            "Element2D4N", 2, { 2, 9, 10, 3 }, nullptr);
    }

    /// Check if search function works properly
    KRATOS_TEST_CASE_IN_SUITE(MPMSearchElement, KratosParticleMechanicsFastSuite)
    {
        Model current_model;
        ModelPart& r_mpm_model_part = current_model.CreateModelPart("MPMModelPart");

        ModelPart& r_background_model_part = current_model.CreateModelPart("MPMBackgroundModelPart");
        PrepareBackgroundModelPart(r_background_model_part);
        PrepareModelPart(r_mpm_model_part, r_background_model_part);

        // First Coordinates of Material Point
        array_1d<double, 3> mp_coordinate;
        mp_coordinate[0] = 0.0;
        mp_coordinate[1] = 0.2;
        mp_coordinate[2] = 0.0;

        r_mpm_model_part.GetElement(2).SetValuesOnIntegrationPoints(
            MP_COORD, { mp_coordinate }, r_mpm_model_part.GetProcessInfo());

        MPMSearchElementUtility::SearchElement<2>(
            r_background_model_part, r_mpm_model_part, 1000, 1e-6);

        // Check nodes
        KRATOS_CHECK_EQUAL(r_mpm_model_part.GetElement(2).GetGeometry()[0].Id(), 1);
        KRATOS_CHECK_EQUAL(r_mpm_model_part.GetElement(2).GetGeometry()[1].Id(), 2);
        KRATOS_CHECK_EQUAL(r_mpm_model_part.GetElement(2).GetGeometry()[2].Id(), 3);
        KRATOS_CHECK_EQUAL(r_mpm_model_part.GetElement(2).GetGeometry()[3].Id(), 4);

        // New Coordinates of Material Point
        mp_coordinate[0] = 1.2;
        mp_coordinate[1] = 0.0;
        mp_coordinate[2] = 0.0;

        r_mpm_model_part.GetElement(2).SetValuesOnIntegrationPoints(
            MP_COORD, { mp_coordinate }, r_mpm_model_part.GetProcessInfo());

        MPMSearchElementUtility::SearchElement<2>(
            r_background_model_part, r_mpm_model_part, 1000, 1e-6);

        std::vector<array_1d<double, 3>> coords;
        r_mpm_model_part.GetElement(2).CalculateOnIntegrationPoints(
            MP_COORD, coords, r_mpm_model_part.GetProcessInfo());
        KRATOS_CHECK_VECTOR_NEAR(coords[0], mp_coordinate, 1e-6)
        // Check nodes
        KRATOS_CHECK_EQUAL(r_mpm_model_part.GetElement(2).GetGeometry()[0].Id(), 2);
        KRATOS_CHECK_EQUAL(r_mpm_model_part.GetElement(2).GetGeometry()[1].Id(), 9);
        KRATOS_CHECK_EQUAL(r_mpm_model_part.GetElement(2).GetGeometry()[2].Id(), 10);
        KRATOS_CHECK_EQUAL(r_mpm_model_part.GetElement(2).GetGeometry()[3].Id(), 3);
    }

} // namespace Testing
} // namespace Kratos
