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
#include "mpm_application_variables.h"
#include "containers/model.h"

#include "custom_utilities/material_point_search_utility.h"
#include "utilities/quadrature_points_utility.h"

namespace Kratos
{
namespace Testing
{
    // Tolerance
    static constexpr double tolerance_pqmpm_weight = 1.0e-4;

    void PrepareModelPart(
        ModelPart& rModelPart,
        ModelPart& rBackgroundModelPart,
        const array_1d<double, 3>& rMPCoords,
        const double IntWeight = 1.0)
    {
        // Properties
        Properties::Pointer p_elem_prop = rModelPart.CreateNewProperties(0);

        // Elements
        auto p_quad = CreateQuadraturePointsUtility<Node>::CreateFromCoordinates(
            rBackgroundModelPart.GetElement(1).pGetGeometry(), rMPCoords, IntWeight);

        const Element& new_element = KratosComponents<Element>::Get("MPMUpdatedLagrangian2D4N");
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

    void PrepareGeneralBackgroundModelPart(ModelPart& rBackgroundModelPart, const int caseIndex )
    {
        // General grid scheme:
        //  13--14--15--16
        //  |  /|  /|  /|
        //  | / | / | / |
        //  |/  |/  |/  |
        //  9---10- 11--12 (if skew then x += 0.1 in this row)
        //  |  /|  /|  /|
        //  | / | / | / |
        //  |/  |/  |/  |
        //  5---6---7---8 (if skew then x += 0.1 in this row)
        //  |  /|  /|  /|
        //  | / | / | / |
        //  |/  |/  |/  |
        //  1---2---3---4

        bool is_skew = false;
        std::string element;
        IndexType number_of_nodes = 16;
        switch (caseIndex)
        {
        case 0:
            element = "Element2D4N";
            break;
        case 1:
            is_skew = true;
            element = "Element2D4N";
            break;
        case 10:
            element = "Element2D3N";
            break;
        case 11:
            is_skew = true;
            element = "Element2D3N";
            break;
        case 20:
            element = "Element3D8N";
            number_of_nodes *= 2;
            break;
        case 21:
            is_skew = true;
            element = "Element3D8N";
            number_of_nodes *= 2;
            break;
        default:
            break;
        }


        // Nodes
        std::vector<Node::Pointer> point_vector(number_of_nodes);
        const double dx = 1.0;
        const double dy = 1.0;
        IndexType point_index = 1;
        for (size_t row = 0; row < 4; ++row) {
            for (size_t col = 0; col < 4; ++col) {
                double myskew = (col > 0 && row == 1 && is_skew) ? 0.1 : 0.0;
                point_vector[point_index-1] = rBackgroundModelPart.CreateNewNode(point_index,
                    double(col)*dx + myskew, double(row)*dy, 0.0);
                point_index += 1;
            }
        }

        if (element == "Element2D4N")
        {
            rBackgroundModelPart.CreateNewElement( element, 1, { 1, 2, 6, 5 }, nullptr);
            rBackgroundModelPart.CreateNewElement( element, 2, { 2, 3, 7, 6 }, nullptr);
            rBackgroundModelPart.CreateNewElement( element, 3, { 3, 4, 8, 7 }, nullptr);

            rBackgroundModelPart.CreateNewElement( element, 4, { 5, 6, 10, 9 }, nullptr);
            rBackgroundModelPart.CreateNewElement( element, 5, { 6, 7, 11, 10 }, nullptr);
            rBackgroundModelPart.CreateNewElement( element, 6, { 7, 8, 12, 11 }, nullptr);

            rBackgroundModelPart.CreateNewElement( element, 7, { 9, 10, 14, 13 }, nullptr);
            rBackgroundModelPart.CreateNewElement( element, 8, { 10, 11, 15, 14 }, nullptr);
            rBackgroundModelPart.CreateNewElement( element, 9, { 11, 12, 16, 15 }, nullptr);
        }
        else if (element == "Element2D3N")
        {
            rBackgroundModelPart.CreateNewElement(element, 1, { 1, 2, 6 }, nullptr);
            rBackgroundModelPart.CreateNewElement(element, 2, { 2, 3, 7 }, nullptr);
            rBackgroundModelPart.CreateNewElement(element, 3, { 3, 4, 8 }, nullptr);

            rBackgroundModelPart.CreateNewElement(element, 4, { 1, 6, 5 }, nullptr);
            rBackgroundModelPart.CreateNewElement(element, 5, { 2, 7, 6 }, nullptr);
            rBackgroundModelPart.CreateNewElement(element, 6, { 3, 8, 7 }, nullptr);

            rBackgroundModelPart.CreateNewElement(element, 7, { 5, 6, 10 }, nullptr);
            rBackgroundModelPart.CreateNewElement(element, 8, { 6, 7, 11 }, nullptr);
            rBackgroundModelPart.CreateNewElement(element, 9, { 7, 8, 12 }, nullptr);

            rBackgroundModelPart.CreateNewElement(element, 10, { 5, 10, 9 }, nullptr);
            rBackgroundModelPart.CreateNewElement(element, 11, { 6, 11, 10 }, nullptr);
            rBackgroundModelPart.CreateNewElement(element, 12, { 7, 12, 11 }, nullptr);

            rBackgroundModelPart.CreateNewElement(element, 13, { 9, 10, 14 }, nullptr);
            rBackgroundModelPart.CreateNewElement(element, 14, { 10, 11, 15 }, nullptr);
            rBackgroundModelPart.CreateNewElement(element, 15, { 11, 12, 16 }, nullptr);

            rBackgroundModelPart.CreateNewElement(element, 16, { 9, 14, 13 }, nullptr);
            rBackgroundModelPart.CreateNewElement(element, 17, { 10, 15, 14 }, nullptr);
            rBackgroundModelPart.CreateNewElement(element, 18, { 11, 16, 15 }, nullptr);
        }
        else if (element == "Element3D8N")
        {
            IndexType point_index = number_of_nodes/2 + 1;
            for (size_t row = 0; row < 4; ++row) {
                for (size_t col = 0; col < 4; ++col) {
                    double myskew = (col > 0 && row == 1 && is_skew) ? 0.1 : 0.0;
                    point_vector[point_index - 1] = rBackgroundModelPart.CreateNewNode(point_index,
                        double(col) * dx + myskew, double(row) * dy, 1.0);
                    point_index += 1;
                }
            }
            rBackgroundModelPart.CreateNewElement(element, 1, { 1, 2, 6, 5, 17,18,22,21}, nullptr);
            rBackgroundModelPart.CreateNewElement(element, 2, { 2, 3, 7, 6 , 18,19,23,22}, nullptr);

            rBackgroundModelPart.CreateNewElement(element, 4, { 5, 6, 10, 9 , 21,22,26,25}, nullptr);
            rBackgroundModelPart.CreateNewElement(element, 5, { 6, 7, 11, 10 ,22,23,27,26}, nullptr);
        }
    }

    ///// Check if search function works properly
    KRATOS_TEST_CASE_IN_SUITE(MPMSearchElementQuad2D, KratosMPMFastSuite)
    {
        // First Coordinates of Material Point
        array_1d<double, 3> mp_coordinate;
        mp_coordinate[0] = 0.0;
        mp_coordinate[1] = 0.2;
        mp_coordinate[2] = 0.0;
        const double int_weight = 1.5;

        Model current_model;
        ModelPart& r_mpm_model_part = current_model.CreateModelPart("MPMModelPart");

        ModelPart& r_background_model_part = current_model.CreateModelPart("MPMBackgroundModelPart");
        PrepareBackgroundModelPart(r_background_model_part);
        PrepareModelPart(r_mpm_model_part, r_background_model_part, mp_coordinate, int_weight);
        const ProcessInfo& r_current_process_info = r_mpm_model_part.GetProcessInfo();

        r_mpm_model_part.GetElement(2).SetValuesOnIntegrationPoints(
            MP_COORD, { mp_coordinate }, r_current_process_info);

        MPMSearchElementUtility::SearchElement<2>(
            r_background_model_part, r_mpm_model_part, 1000, 1e-6);

        // Check nodes
        KRATOS_EXPECT_EQ(r_mpm_model_part.GetElement(2).GetGeometry()[0].Id(), 1);
        KRATOS_EXPECT_EQ(r_mpm_model_part.GetElement(2).GetGeometry()[1].Id(), 2);
        KRATOS_EXPECT_EQ(r_mpm_model_part.GetElement(2).GetGeometry()[2].Id(), 3);
        KRATOS_EXPECT_EQ(r_mpm_model_part.GetElement(2).GetGeometry()[3].Id(), 4);

        // New Coordinates of Material Point
        mp_coordinate[0] = 1.2;
        mp_coordinate[1] = 0.0;
        mp_coordinate[2] = 0.0;

        r_mpm_model_part.GetElement(2).SetValuesOnIntegrationPoints(
            MP_COORD, { mp_coordinate }, r_current_process_info);

        MPMSearchElementUtility::SearchElement<2>(
            r_background_model_part, r_mpm_model_part, 1000, 1e-6);

        std::vector<array_1d<double, 3>> coords;
        r_mpm_model_part.GetElement(2).CalculateOnIntegrationPoints(
            MP_COORD, coords, r_current_process_info);
        KRATOS_EXPECT_VECTOR_NEAR(coords[0], mp_coordinate, 1e-6)
        // Check nodes
        KRATOS_EXPECT_EQ(r_mpm_model_part.GetElement(2).GetGeometry()[0].Id(), 2);
        KRATOS_EXPECT_EQ(r_mpm_model_part.GetElement(2).GetGeometry()[1].Id(), 9);
        KRATOS_EXPECT_EQ(r_mpm_model_part.GetElement(2).GetGeometry()[2].Id(), 10);
        KRATOS_EXPECT_EQ(r_mpm_model_part.GetElement(2).GetGeometry()[3].Id(), 3);
    }

    ///// Check if search function works properly
    KRATOS_TEST_CASE_IN_SUITE(MPMSearchElementTri2D, KratosMPMFastSuite)
    {
        // First Coordinates of Material Point: (0.0, 0.2, 0.0)
        array_1d<double, 3> mp_coordinate;
        mp_coordinate[0] = 0.0;
        mp_coordinate[1] = 0.2;
        mp_coordinate[2] = 0.0;
        const double int_weight = 1.5;

        const IndexType grid_case = 10;

        Model current_model;
        ModelPart& r_mpm_model_part = current_model.CreateModelPart("MPMModelPart");

        ModelPart& r_background_model_part = current_model.CreateModelPart("MPMBackgroundModelPart");
        PrepareGeneralBackgroundModelPart(r_background_model_part, grid_case);
        PrepareModelPart(r_mpm_model_part, r_background_model_part, mp_coordinate, int_weight);

        const ProcessInfo& r_current_process_info = r_mpm_model_part.GetProcessInfo();
        const double search_tolerance = 1e-7;
        const double check_tolerance = 1e-10;

        r_mpm_model_part.GetElement(2).SetValuesOnIntegrationPoints(
            MP_COORD, { mp_coordinate }, r_current_process_info);

        MPMSearchElementUtility::SearchElement<2>(
            r_background_model_part, r_mpm_model_part, 1000, search_tolerance);

        std::vector<array_1d<double, 3>> coords;
        r_mpm_model_part.GetElement(2).CalculateOnIntegrationPoints(
            MP_COORD, coords, r_current_process_info);
        KRATOS_EXPECT_VECTOR_NEAR(coords[0], mp_coordinate, 1e-6)

        // Check nodes
        KRATOS_EXPECT_EQ(r_mpm_model_part.GetElement(2).GetGeometry()[0].Id(), 1);
        KRATOS_EXPECT_EQ(r_mpm_model_part.GetElement(2).GetGeometry()[1].Id(), 6);
        KRATOS_EXPECT_EQ(r_mpm_model_part.GetElement(2).GetGeometry()[2].Id(), 5);

        Matrix shape_functions_values = r_mpm_model_part.GetElement(2).GetGeometry().ShapeFunctionsValues();
        KRATOS_EXPECT_NEAR(shape_functions_values(0,0), 0.8, check_tolerance);
        KRATOS_EXPECT_NEAR(shape_functions_values(0,1), 0.0, check_tolerance);
        KRATOS_EXPECT_NEAR(shape_functions_values(0,2), 0.2, check_tolerance);

        // New Coordinates of Material Point
        mp_coordinate[0] = 1.2;
        mp_coordinate[1] = 0.0;
        mp_coordinate[2] = 0.0;

        r_mpm_model_part.GetElement(2).SetValuesOnIntegrationPoints(
            MP_COORD, { mp_coordinate }, r_current_process_info);

        MPMSearchElementUtility::SearchElement<2>(
            r_background_model_part, r_mpm_model_part, 1000, search_tolerance);

        r_mpm_model_part.GetElement(2).CalculateOnIntegrationPoints(
            MP_COORD, coords, r_current_process_info);
        KRATOS_EXPECT_VECTOR_NEAR(coords[0], mp_coordinate, 1e-6)

        // Check nodes
        KRATOS_EXPECT_EQ(r_mpm_model_part.GetElement(2).GetGeometry()[0].Id(), 2);
        KRATOS_EXPECT_EQ(r_mpm_model_part.GetElement(2).GetGeometry()[1].Id(), 3);
        KRATOS_EXPECT_EQ(r_mpm_model_part.GetElement(2).GetGeometry()[2].Id(), 7);

        shape_functions_values = r_mpm_model_part.GetElement(2).GetGeometry().ShapeFunctionsValues();
        KRATOS_EXPECT_NEAR(shape_functions_values(0,0), 0.8, check_tolerance);
        KRATOS_EXPECT_NEAR(shape_functions_values(0,1), 0.2, check_tolerance);
        KRATOS_EXPECT_NEAR(shape_functions_values(0,2), 0.0, check_tolerance);

        // New Coordinates of Material Point
        mp_coordinate[0] = 2.6;
        mp_coordinate[1] = 2.1;
        mp_coordinate[2] = 0.0;

        r_mpm_model_part.GetElement(2).SetValuesOnIntegrationPoints(
            MP_COORD, { mp_coordinate }, r_current_process_info);

        MPMSearchElementUtility::SearchElement<2>(
            r_background_model_part, r_mpm_model_part, 1000, search_tolerance);

        r_mpm_model_part.GetElement(2).CalculateOnIntegrationPoints(
            MP_COORD, coords, r_current_process_info);
        KRATOS_EXPECT_VECTOR_NEAR(coords[0], mp_coordinate, 1e-6)

        // Check nodes
        KRATOS_EXPECT_EQ(r_mpm_model_part.GetElement(2).GetGeometry()[0].Id(), 11);
        KRATOS_EXPECT_EQ(r_mpm_model_part.GetElement(2).GetGeometry()[1].Id(), 12);
        KRATOS_EXPECT_EQ(r_mpm_model_part.GetElement(2).GetGeometry()[2].Id(), 16);

        shape_functions_values = r_mpm_model_part.GetElement(2).GetGeometry().ShapeFunctionsValues();
        KRATOS_EXPECT_NEAR(shape_functions_values(0,0), 0.4, check_tolerance);
        KRATOS_EXPECT_NEAR(shape_functions_values(0,1), 0.5, check_tolerance);
        KRATOS_EXPECT_NEAR(shape_functions_values(0,2), 0.1, check_tolerance);

        // New Coordinates of Material Point
        mp_coordinate[0] = 0.2349;
        mp_coordinate[1] = 2.7238;
        mp_coordinate[2] = 0.0;

        r_mpm_model_part.GetElement(2).SetValuesOnIntegrationPoints(
            MP_COORD, { mp_coordinate }, r_current_process_info);

        MPMSearchElementUtility::SearchElement<2>(
            r_background_model_part, r_mpm_model_part, 1000, search_tolerance);

        r_mpm_model_part.GetElement(2).CalculateOnIntegrationPoints(
            MP_COORD, coords, r_current_process_info);
        KRATOS_EXPECT_VECTOR_NEAR(coords[0], mp_coordinate, 1e-6)

        // Check nodes
        KRATOS_EXPECT_EQ(r_mpm_model_part.GetElement(2).GetGeometry()[0].Id(), 9);
        KRATOS_EXPECT_EQ(r_mpm_model_part.GetElement(2).GetGeometry()[1].Id(), 14);
        KRATOS_EXPECT_EQ(r_mpm_model_part.GetElement(2).GetGeometry()[2].Id(), 13);

        shape_functions_values = r_mpm_model_part.GetElement(2).GetGeometry().ShapeFunctionsValues();
        KRATOS_EXPECT_NEAR(shape_functions_values(0,0), 0.2762, check_tolerance);
        KRATOS_EXPECT_NEAR(shape_functions_values(0,1), 0.2349, check_tolerance);
        KRATOS_EXPECT_NEAR(shape_functions_values(0,2), 0.4889, check_tolerance);
    }

    /// PQMPM test 1 - MP domain is entirely within cell and should only make 1 mp
    KRATOS_TEST_CASE_IN_SUITE(MPMSearchElementPQMPM2DQuadWithinCell, KratosMPMFastSuite)
    {
        // First Coordinates of Material Point
        array_1d<double, 3> mp_coordinate;
        mp_coordinate[0] = 0.5;
        mp_coordinate[1] = 0.5;
        mp_coordinate[2] = 0.0;
        const double int_weight = 1.0;
        std::vector<double> mp_volume (1);
        mp_volume[0] = 1.0;

        // Case of background grid( 0=quad, 10=tri, 20=hex. +1 = skew)
        const IndexType grid_case = 0;

        Model current_model;
        ModelPart& r_mpm_model_part = current_model.CreateModelPart("MPMModelPart");

        ModelPart& r_background_model_part = current_model.CreateModelPart("MPMBackgroundModelPart");
        PrepareGeneralBackgroundModelPart(r_background_model_part, grid_case);
        PrepareModelPart(r_mpm_model_part, r_background_model_part, mp_coordinate, int_weight);

        r_background_model_part.GetProcessInfo().SetValue(IS_PQMPM, true);
        r_background_model_part.GetProcessInfo().SetValue(IS_MAKE_NORMAL_MP_IF_PQMPM_FAILS, false);
        r_background_model_part.GetProcessInfo().SetValue(PQMPM_SUBPOINT_MIN_VOLUME_FRACTION, 1e-24);

        const ProcessInfo& r_current_process_info = r_mpm_model_part.GetProcessInfo();

        r_mpm_model_part.GetElement(2).SetValuesOnIntegrationPoints(
            MP_COORD, { mp_coordinate }, r_current_process_info);
        r_mpm_model_part.GetElement(2).SetValuesOnIntegrationPoints(
            MP_VOLUME, mp_volume, r_current_process_info);

        MPMSearchElementUtility::SearchElement<2>(
            r_background_model_part, r_mpm_model_part, 1000, 1e-6);

        // Check results - MP should lie entirely within one cell
        Geometry<Node>& rGeom = r_mpm_model_part.GetElement(2).GetGeometry();
        KRATOS_EXPECT_EQ(rGeom.IntegrationPointsNumber(), 1);
        KRATOS_EXPECT_NEAR(rGeom.IntegrationPoints()[0].Weight(),1.0, std::numeric_limits<double>::epsilon());
    }

    /// PQMPM test 2 - MP domain is entirely within cell and should only make 1 mp
    KRATOS_TEST_CASE_IN_SUITE(MPMSearchElementPQMPM2DTriWithinCell, KratosMPMFastSuite)
    {
        // First Coordinates of Material Point
        array_1d<double, 3> mp_coordinate;
        mp_coordinate[0] = 0.8;
        mp_coordinate[1] = 0.1;
        mp_coordinate[2] = 0.0;
        const double int_weight = 1.0;
        std::vector<double> mp_volume (1);
        mp_volume[0] = 0.01;

        // Case of background grid( 0=quad, 10=tri, 20=hex. +1 = skew)
        const IndexType grid_case = 10;

        Model current_model;
        ModelPart& r_mpm_model_part = current_model.CreateModelPart("MPMModelPart");

        ModelPart& r_background_model_part = current_model.CreateModelPart("MPMBackgroundModelPart");
        PrepareGeneralBackgroundModelPart(r_background_model_part, grid_case);
        PrepareModelPart(r_mpm_model_part, r_background_model_part, mp_coordinate, int_weight);

        r_background_model_part.GetProcessInfo().SetValue(IS_PQMPM, true);
        r_background_model_part.GetProcessInfo().SetValue(IS_MAKE_NORMAL_MP_IF_PQMPM_FAILS, false);
        r_background_model_part.GetProcessInfo().SetValue(PQMPM_SUBPOINT_MIN_VOLUME_FRACTION, 1e-24);

        const ProcessInfo& r_current_process_info = r_mpm_model_part.GetProcessInfo();

        r_mpm_model_part.GetElement(2).SetValuesOnIntegrationPoints(
            MP_COORD, { mp_coordinate }, r_current_process_info);
        r_mpm_model_part.GetElement(2).SetValuesOnIntegrationPoints(
            MP_VOLUME, mp_volume, r_current_process_info);

        MPMSearchElementUtility::SearchElement<2>(
            r_background_model_part, r_mpm_model_part, 1000, 1e-6);

        // Check results - MP should lie entirely within one cell
        Geometry<Node>& rGeom = r_mpm_model_part.GetElement(2).GetGeometry();
        KRATOS_EXPECT_EQ(rGeom.IntegrationPointsNumber(), 1);
        KRATOS_EXPECT_NEAR(rGeom.IntegrationPoints()[0].Weight(),1.0, std::numeric_limits<double>::epsilon());
    }

    /// PQMPM test 3 - MP domain is not entirely within one cell and should make 9 mp. 2D skew quad
    KRATOS_TEST_CASE_IN_SUITE(MPMSearchElementPQMPM2DQuad, KratosMPMFastSuite)
    {
        // First Coordinates of Material Point
        array_1d<double, 3> mp_coordinate;
        mp_coordinate[0] = 1.5;
        mp_coordinate[1] = 1.5;
        mp_coordinate[2] = 0.0;
        const double int_weight = 1.0;
        std::vector<double> mp_volume(1);
        mp_volume[0] = 2.0;

        // Case of background grid( 0=quad, 10=tri, 20=hex. +1 = skew)
        const IndexType grid_case = 1;

        Model current_model;
        ModelPart& r_mpm_model_part = current_model.CreateModelPart("MPMModelPart");

        ModelPart& r_background_model_part = current_model.CreateModelPart("MPMBackgroundModelPart");
        PrepareGeneralBackgroundModelPart(r_background_model_part, grid_case);
        PrepareModelPart(r_mpm_model_part, r_background_model_part, mp_coordinate, int_weight);

        r_background_model_part.GetProcessInfo().SetValue(IS_PQMPM, true);
        r_background_model_part.GetProcessInfo().SetValue(IS_MAKE_NORMAL_MP_IF_PQMPM_FAILS, false);

        const ProcessInfo& r_current_process_info = r_mpm_model_part.GetProcessInfo();

        r_mpm_model_part.GetElement(2).SetValuesOnIntegrationPoints(
            MP_COORD, { mp_coordinate }, r_current_process_info);
        r_mpm_model_part.GetElement(2).SetValuesOnIntegrationPoints(
            MP_VOLUME, mp_volume, r_current_process_info);

        MPMSearchElementUtility::SearchElement<2>(
            r_background_model_part, r_mpm_model_part, 1000, 1e-6);

        // Check results
        Geometry<Node>& rGeom = r_mpm_model_part.GetElement(2).GetGeometry();
        //KRATOS_WATCH(rGeom.IntegrationPointsNumber())
        KRATOS_EXPECT_EQ(rGeom.IntegrationPointsNumber(), 9);
        std::vector<double> result_weight = { 0.5, 0.0307296, 0.103553, 0.0121636, 0.0785534,
            0.103553, 0.128553, 0.0214466, 0.0214466 };
        for (size_t i = 0; i < rGeom.IntegrationPointsNumber(); i++)
        {
            //KRATOS_WATCH(rGeom.IntegrationPoints()[i].Weight())
            KRATOS_EXPECT_NEAR(rGeom.IntegrationPoints()[i].Weight(), result_weight[i], tolerance_pqmpm_weight);
        }
    }

    /// PQMPM test 4 - MP domain is not entirely within one cell and should make 5 mp. 2D skew tri
    KRATOS_TEST_CASE_IN_SUITE(MPMSearchElementPQMPM2DTri, KratosMPMFastSuite)
    {
        // First Coordinates of Material Point
        array_1d<double, 3> mp_coordinate;
        mp_coordinate[0] = 1.5;
        mp_coordinate[1] = 1.9;
        mp_coordinate[2] = 0.0;
        const double int_weight = 1.0;
        std::vector<double> mp_volume(1);
        mp_volume[0] = 1.0;

        // Case of background grid( 0=quad, 10=tri, 20=hex. +1 = skew)
        const IndexType grid_case = 11;

        Model current_model;
        ModelPart& r_mpm_model_part = current_model.CreateModelPart("MPMModelPart");

        ModelPart& r_background_model_part = current_model.CreateModelPart("MPMBackgroundModelPart");
        PrepareGeneralBackgroundModelPart(r_background_model_part, grid_case);
        PrepareModelPart(r_mpm_model_part, r_background_model_part, mp_coordinate, int_weight);

        r_background_model_part.GetProcessInfo().SetValue(IS_PQMPM, true);
        r_background_model_part.GetProcessInfo().SetValue(IS_MAKE_NORMAL_MP_IF_PQMPM_FAILS, false);

        const ProcessInfo& r_current_process_info = r_mpm_model_part.GetProcessInfo();

        r_mpm_model_part.GetElement(2).SetValuesOnIntegrationPoints(
            MP_COORD, { mp_coordinate }, r_current_process_info);
        r_mpm_model_part.GetElement(2).SetValuesOnIntegrationPoints(
            MP_VOLUME, mp_volume, r_current_process_info);

        MPMSearchElementUtility::SearchElement<2>(
            r_background_model_part, r_mpm_model_part, 1000, 1e-6);

        // Check results
        Geometry<Node>& rGeom = r_mpm_model_part.GetElement(2).GetGeometry();
        //KRATOS_WATCH(rGeom.IntegrationPointsNumber())
        KRATOS_EXPECT_EQ(rGeom.IntegrationPointsNumber(), 5);
        std::vector<double> result_weight = { 0.42, 0.018, 0.162, 0.32, 0.08};
        for (size_t i = 0; i < rGeom.IntegrationPointsNumber(); i++)
        {
            //KRATOS_WATCH(rGeom.IntegrationPoints()[i].Weight())
            KRATOS_EXPECT_NEAR(rGeom.IntegrationPoints()[i].Weight(), result_weight[i], tolerance_pqmpm_weight);
        }
    }

    /// PQMPM test 5 - MP domain is not entirely within one cell and should make 4 mp. 3D hex
    KRATOS_TEST_CASE_IN_SUITE(MPMSearchElementPQMPM3DHex, KratosMPMFastSuite)
    {
        // First Coordinates of Material Point
        array_1d<double, 3> mp_coordinate;
        mp_coordinate[0] = 0.8;
        mp_coordinate[1] = 0.8;
        mp_coordinate[2] = 0.5;
        const double int_weight = 1.0;
        std::vector<double> mp_volume(1);
        mp_volume[0] = 0.8;

        // Case of background grid( 0=quad, 10=tri, 20=hex. +1 = skew)
        const IndexType grid_case = 20;

        Model current_model;
        ModelPart& r_mpm_model_part = current_model.CreateModelPart("MPMModelPart");

        ModelPart& r_background_model_part = current_model.CreateModelPart("MPMBackgroundModelPart");
        PrepareGeneralBackgroundModelPart(r_background_model_part, grid_case);
        PrepareModelPart(r_mpm_model_part, r_background_model_part, mp_coordinate, int_weight);

        r_background_model_part.GetProcessInfo().SetValue(IS_PQMPM, true);
        r_background_model_part.GetProcessInfo().SetValue(IS_MAKE_NORMAL_MP_IF_PQMPM_FAILS, false);

        const ProcessInfo& r_current_process_info = r_mpm_model_part.GetProcessInfo();

        r_mpm_model_part.GetElement(2).SetValuesOnIntegrationPoints(
            MP_COORD, { mp_coordinate }, r_current_process_info);
        r_mpm_model_part.GetElement(2).SetValuesOnIntegrationPoints(
            MP_VOLUME, mp_volume, r_current_process_info);

        MPMSearchElementUtility::SearchElement<3>(
            r_background_model_part, r_mpm_model_part, 1000, 1e-6);

        // Check results
        Geometry<Node>& rGeom = r_mpm_model_part.GetElement(2).GetGeometry();
        //KRATOS_WATCH(rGeom.IntegrationPointsNumber())
        KRATOS_EXPECT_EQ(rGeom.IntegrationPointsNumber(), 4);
        std::vector<double> result_weight = { 0.511859 , 0.203584 , 0.203584 , 0.0809724 };
        for (size_t i = 0; i < rGeom.IntegrationPointsNumber(); i++)
        {
            //KRATOS_WATCH(rGeom.IntegrationPoints()[i].Weight())
            KRATOS_EXPECT_NEAR(rGeom.IntegrationPoints()[i].Weight(), result_weight[i], tolerance_pqmpm_weight);
        }
    }

    // TEST DISABLED - KRATOS_EXPECT_EXCEPTION_IS_THROWN doesn't work in parallel!
    /// PQMPM test 6 - Check pqmpm fails with unstructured 3D mesh
    //KRATOS_TEST_CASE_IN_SUITE(MPMSearchElementPQMPM3DHexError, KratosMPMFastSuite)
    //{
    //    // First Coordinates of Material Point
    //    array_1d<double, 3> mp_coordinate;
    //    mp_coordinate[0] = 0.8;
    //    mp_coordinate[1] = 0.8;
    //    mp_coordinate[2] = 0.5;
    //    const double int_weight = 1.0;
    //    std::vector<double> mp_volume(1);
    //    mp_volume[0] = 0.8;
    //
    //    // Case of background grid( 0=quad, 10=tri, 20=hex. +1 = skew)
    //    const IndexType grid_case = 21;
    //
    //    Model current_model;
    //    ModelPart& r_mpm_model_part = current_model.CreateModelPart("MPMModelPart");
    //
    //    ModelPart& r_background_model_part = current_model.CreateModelPart("MPMBackgroundModelPart");
    //    PrepareGeneralBackgroundModelPart(r_background_model_part, grid_case);
    //    PrepareModelPart(r_mpm_model_part, r_background_model_part, mp_coordinate, int_weight);
    //
    //    r_background_model_part.GetProcessInfo().SetValue(IS_PQMPM, true);
    //    r_background_model_part.GetProcessInfo().SetValue(IS_MAKE_NORMAL_MP_IF_PQMPM_FAILS, false);

    // const ProcessInfo& r_current_process_info = r_mpm_model_part.GetProcessInfo();
    //
    //    r_mpm_model_part.GetElement(2).SetValuesOnIntegrationPoints(
    //        MP_COORD, { mp_coordinate }, r_current_process_info);
    //    r_mpm_model_part.GetElement(2).SetValuesOnIntegrationPoints(
    //        MP_VOLUME, mp_volume, r_current_process_info);
    //
    //    KRATOS_EXPECT_EXCEPTION_IS_THROWN(MPMSearchElementUtility::SearchElement<2>(
    //        r_background_model_part, r_mpm_model_part, 1000, 1e-6) ,
    //        "Error: ERROR")
    //}

    // TEST DISABLED - KRATOS_EXPECT_EXCEPTION_IS_THROWN doesn't work in parallel!
    ///// PQMPM test 7 - Check pqmpm fails if point tries to split outside mesh
    //KRATOS_TEST_CASE_IN_SUITE(MPMSearchElementPQMPM2DOutsideSplit, KratosMPMFastSuite)
    //{
    //    // First Coordinates of Material Point
    //    array_1d<double, 3> mp_coordinate;
    //    mp_coordinate[0] = 0.1;
    //    mp_coordinate[1] = 1.1;
    //    mp_coordinate[2] = 0.0;
    //    const double int_weight = 1.0;
    //    std::vector<double> mp_volume(1);
    //    mp_volume[0] = 1.0;
    //
    //    // Case of background grid( 0=quad, 10=tri, 20=hex. +1 = skew)
    //    const IndexType grid_case = 0;
    //
    //    Model current_model;
    //    ModelPart& r_mpm_model_part = current_model.CreateModelPart("MPMModelPart");
    //
    //    ModelPart& r_background_model_part = current_model.CreateModelPart("MPMBackgroundModelPart");
    //    PrepareGeneralBackgroundModelPart(r_background_model_part, grid_case);
    //    PrepareModelPart(r_mpm_model_part, r_background_model_part, mp_coordinate, int_weight);
    //
    //    r_background_model_part.GetProcessInfo().SetValue(IS_PQMPM, true);
    //    r_background_model_part.GetProcessInfo().SetValue(IS_MAKE_NORMAL_MP_IF_PQMPM_FAILS, false);

    // const ProcessInfo& r_current_process_info = r_mpm_model_part.GetProcessInfo();
    //
    //    r_mpm_model_part.GetElement(2).SetValuesOnIntegrationPoints(
    //        MP_COORD, { mp_coordinate }, r_current_process_info);
    //    r_mpm_model_part.GetElement(2).SetValuesOnIntegrationPoints(
    //        MP_VOLUME, mp_volume, r_current_process_info);
    //
    //    KRATOS_EXPECT_EXCEPTION_IS_THROWN(MPMSearchElementUtility::SearchElement<2>(
    //        r_background_model_part, r_mpm_model_part, 1000, 1e-6) ,
    //        "Error: ERROR")
    //}

    /// PQMPM test 8 - Check pqmpm falls back to normal mpm if point tries to split outside mesh
    KRATOS_TEST_CASE_IN_SUITE(MPMSearchElementPQMPM2DFallbackToMPM, KratosMPMFastSuite)
    {
        // First Coordinates of Material Point
        array_1d<double, 3> mp_coordinate;
        mp_coordinate[0] = 0.1;
        mp_coordinate[1] = 0.1;
        mp_coordinate[2] = 0.0;
        const double int_weight = 1.0;
        std::vector<double> mp_volume(1);
        mp_volume[0] = 0.8;

        // Case of background grid( 0=quad, 10=tri, 20=hex. +1 = skew)
        const IndexType grid_case = 0;

        Model current_model;
        ModelPart& r_mpm_model_part = current_model.CreateModelPart("MPMModelPart");

        ModelPart& r_background_model_part = current_model.CreateModelPart("MPMBackgroundModelPart");
        PrepareGeneralBackgroundModelPart(r_background_model_part, grid_case);
        PrepareModelPart(r_mpm_model_part, r_background_model_part, mp_coordinate, int_weight);

        r_background_model_part.GetProcessInfo().SetValue(IS_PQMPM, true);
        r_background_model_part.GetProcessInfo().SetValue(IS_MAKE_NORMAL_MP_IF_PQMPM_FAILS, true);

        const ProcessInfo& r_current_process_info = r_mpm_model_part.GetProcessInfo();

        r_mpm_model_part.GetElement(2).SetValuesOnIntegrationPoints(
            MP_COORD, { mp_coordinate }, r_current_process_info);
        r_mpm_model_part.GetElement(2).SetValuesOnIntegrationPoints(
            MP_VOLUME, mp_volume, r_current_process_info);

        Geometry<Node>& rGeom = r_mpm_model_part.GetElement(2).GetGeometry();
        KRATOS_EXPECT_EQ(rGeom.IntegrationPointsNumber(), 1);
        KRATOS_EXPECT_NEAR(rGeom.IntegrationPoints()[0].Weight(), 1.0, std::numeric_limits<double>::epsilon());
    }

    /// PQMPM test 9 - Check pqmpm does not split across a BC
    KRATOS_TEST_CASE_IN_SUITE(MPMSearchElementPQMPM2DQuadBC, KratosMPMFastSuite)
    {
        // First Coordinates of Material Point
        array_1d<double, 3> mp_coordinate;
        mp_coordinate[0] = 0.9;
        mp_coordinate[1] = 0.9;
        mp_coordinate[2] = 0.0;
        const double int_weight = 1.0;
        std::vector<double> mp_volume(1);
        mp_volume[0] = 1.0;

        // Case of background grid( 0=quad, 10=tri, 20=hex. +1 = skew)
        const IndexType grid_case = 0;

        Model current_model;
        ModelPart& r_mpm_model_part = current_model.CreateModelPart("MPMModelPart");

        ModelPart& r_background_model_part = current_model.CreateModelPart("MPMBackgroundModelPart");
        r_background_model_part.AddNodalSolutionStepVariable(DISPLACEMENT_X);
        r_background_model_part.AddNodalSolutionStepVariable(DISPLACEMENT_Y);
        PrepareGeneralBackgroundModelPart(r_background_model_part, grid_case);
        PrepareModelPart(r_mpm_model_part, r_background_model_part, mp_coordinate, int_weight);
        r_background_model_part.pGetNode(6)->Fix(DISPLACEMENT_X);

        r_background_model_part.GetProcessInfo().SetValue(IS_PQMPM, true);
        r_background_model_part.GetProcessInfo().SetValue(IS_MAKE_NORMAL_MP_IF_PQMPM_FAILS, false);

        const ProcessInfo& r_current_process_info = r_mpm_model_part.GetProcessInfo();

        r_mpm_model_part.GetElement(2).SetValuesOnIntegrationPoints(
            MP_COORD, { mp_coordinate }, r_current_process_info);
        r_mpm_model_part.GetElement(2).SetValuesOnIntegrationPoints(
            MP_VOLUME, mp_volume, r_current_process_info);

        MPMSearchElementUtility::SearchElement<2>(
            r_background_model_part, r_mpm_model_part, 1000, 1e-6);

        // Check results - MP should lie entirely within one cell
        Geometry<Node>& rGeom = r_mpm_model_part.GetElement(2).GetGeometry();
        KRATOS_EXPECT_EQ(rGeom.IntegrationPointsNumber(), 1);
        KRATOS_EXPECT_NEAR(rGeom.IntegrationPoints()[0].Weight(), 1.0, std::numeric_limits<double>::epsilon());
    }

    /// PQMPM test 10 - Check pqmpm reverts back to normal MPM if min_pqmpm fraction is specified and not fulfilled
    KRATOS_TEST_CASE_IN_SUITE(MPMSearchElementPQMPM2DFraction, KratosMPMFastSuite)
    {
        // First Coordinates of Material Point
        array_1d<double, 3> mp_coordinate;
        mp_coordinate[0] = 0.6;
        mp_coordinate[1] = 0.6;
        mp_coordinate[2] = 0.0;
        const double int_weight = 1.0;
        std::vector<double> mp_volume(1);
        mp_volume[0] = 1.0;

        // Case of background grid( 0=quad, 10=tri, 20=hex. +1 = skew)
        const IndexType grid_case = 0;

        Model current_model;
        ModelPart& r_mpm_model_part = current_model.CreateModelPart("MPMModelPart");

        ModelPart& r_background_model_part = current_model.CreateModelPart("MPMBackgroundModelPart");
        PrepareGeneralBackgroundModelPart(r_background_model_part, grid_case);
        PrepareModelPart(r_mpm_model_part, r_background_model_part, mp_coordinate, int_weight);

        r_background_model_part.GetProcessInfo().SetValue(IS_PQMPM, true);
        r_background_model_part.GetProcessInfo().SetValue(IS_MAKE_NORMAL_MP_IF_PQMPM_FAILS, true);
        r_background_model_part.GetProcessInfo().SetValue(PQMPM_SUBPOINT_MIN_VOLUME_FRACTION, 0.5);

        const ProcessInfo& r_current_process_info = r_mpm_model_part.GetProcessInfo();

        r_mpm_model_part.GetElement(2).SetValuesOnIntegrationPoints(
            MP_COORD, { mp_coordinate }, r_current_process_info);
        r_mpm_model_part.GetElement(2).SetValuesOnIntegrationPoints(
            MP_VOLUME, mp_volume, r_current_process_info);

        MPMSearchElementUtility::SearchElement<2>(
            r_background_model_part, r_mpm_model_part, 1000, 1e-6);

        // Check results - MP should lie entirely within one cell
        Geometry<Node>& rGeom = r_mpm_model_part.GetElement(2).GetGeometry();
        KRATOS_EXPECT_EQ(rGeom.IntegrationPointsNumber(), 1);
        KRATOS_EXPECT_NEAR(rGeom.IntegrationPoints()[0].Weight(), 1.0, std::numeric_limits<double>::epsilon());
    }

} // namespace Testing
} // namespace Kratos
