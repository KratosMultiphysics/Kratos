//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Andi Katili
//


// System includes
#include <vector>

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "mpm_application_variables.h"
#include "custom_utilities/material_point_search_utility.h"
#include "custom_utilities/mapping_utilities/mpm_flip_particle_mapping_utility.hpp"


namespace Kratos::Testing
{
    void Prepare2D1EBackgroundModelPart(ModelPart& rBackgroundModelPart)
    {
        // Grid scheme:
        //  13---14---15---16
        //  |    |    |    |
        //  |    |    |    |
        //  |    |    |    |
        //  9----10---11----12 
        //  |    |    |    |
        //  |    |    |    |
        //  |    |    |    |
        //  5----6----7----8 
        //  |    |    |    |
        //  |    |    |    |
        //  |    |    |    |
        //  1----2----3----4
        // Nodes
        auto p_node_1 = rBackgroundModelPart.CreateNewNode( 1,  0.0 ,  0.0 , 0.0);
        auto p_node_2 = rBackgroundModelPart.CreateNewNode( 2,  1.0 ,  0.0 , 0.0);
        auto p_node_3 = rBackgroundModelPart.CreateNewNode( 3,  1.0 ,  1.0 , 0.0);
        auto p_node_4 = rBackgroundModelPart.CreateNewNode( 4,  0.0 ,  1.0 , 0.0);

        rBackgroundModelPart.CreateNewElement(
            "Element2D4N", 1, { 1, 2, 3, 4 }, nullptr);
    }

    template <SizeType TDimension>
    void CreateMP(ModelPart& rModelPart, ModelPart& rBackgroundModelPart, Properties::Pointer pProperties, const array_1d<double, 3>& rMPCoordinate, const double& rMPVolume)
    {
        // Create new material point element
        unsigned int new_element_id = rModelPart.NumberOfElements() + 1;

        BinBasedFastPointLocator<TDimension> SearchStructure(rBackgroundModelPart);
        SearchStructure.UpdateSearchDatabase();
        typename BinBasedFastPointLocator<TDimension>::ResultContainerType results(100);
        typename BinBasedFastPointLocator<TDimension>::ResultIteratorType result_begin = results.begin();
        Element::Pointer pelem;
        Vector N;

        bool is_found = SearchStructure.FindPointOnMesh(rMPCoordinate, N, pelem, result_begin);
        if (!is_found) KRATOS_WARNING("MaterialPointGeneratorUtility") << "::search failed." << std::endl;

        auto p_new_mp_geometry = CreateQuadraturePointsUtility<Node>::CreateFromCoordinates(
                            pelem->pGetGeometry(), rMPCoordinate, rMPVolume);

        const Element& element = KratosComponents<Element>::Get("MPMUpdatedLagrangian");
        Element::Pointer p_new_element = element.Create(new_element_id, p_new_mp_geometry, pProperties);

        p_new_element->SetValuesOnIntegrationPoints(MP_COORD, {rMPCoordinate}, rModelPart.GetProcessInfo());
        rModelPart.AddElement(p_new_element);
    }
    template <SizeType TDimension>
    void PrepareMP(
        ModelPart& rModelPart,
        ModelPart& rBackgroundModelPart, std::vector<array_1d<double, 3>>& rMPCoordinates, 
        const double& rMPVolume)
    {
        // Properties
        Properties::Pointer p_elem_prop = rModelPart.pGetProperties(0);

        // Elements
        for (auto& mp_coordinate : rMPCoordinates)
        {
            // array_1d<double, 3> mp_coordinate1{0.211324865,0.211324865, 0.0};
            CreateMP<TDimension>(rModelPart, rBackgroundModelPart, p_elem_prop, mp_coordinate, rMPVolume);
        }
    }

    void Prepare2D1EModelPart(ModelPart& rModelPart, ModelPart& rBackgroundModelPart)
    {
        // Properties
        Properties::Pointer p_elem_prop = rModelPart.CreateNewProperties(0);

        // Elements
        array_1d<double, 3> mp_coordinate1{0.211324865,0.211324865, 0.0};
        array_1d<double, 3> mp_coordinate2{0.788675135,0.211324865, 0.0};
        array_1d<double, 3> mp_coordinate3{0.788675135,0.788675135, 0.0};
        array_1d<double, 3> mp_coordinate4{0.211324865,0.788675135, 0.0};
        std::vector<array_1d<double, 3>> mp_coordinates = {mp_coordinate1, mp_coordinate2, mp_coordinate3, mp_coordinate4};


        PrepareMP<2>(rModelPart, rBackgroundModelPart, mp_coordinates, 0.25);

        auto pElement1 = rModelPart.pGetElement(1);
        auto pElement2 = rModelPart.pGetElement(2);
        auto pElement3 = rModelPart.pGetElement(3);
        auto pElement4 = rModelPart.pGetElement(4);

        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // mass
        double mp_mass1{ 2.5 };
        pElement1->SetValuesOnIntegrationPoints(MP_MASS, { mp_mass1 }, r_current_process_info);
        double mp_mass2{ 2.5 };
        pElement2->SetValuesOnIntegrationPoints(MP_MASS, { mp_mass2 }, r_current_process_info);
        double mp_mass3{ 2.5 };
        pElement3->SetValuesOnIntegrationPoints(MP_MASS, { mp_mass3 }, r_current_process_info);
        double mp_mass4{ 2.5 };
        pElement4->SetValuesOnIntegrationPoints(MP_MASS, { mp_mass4 }, r_current_process_info);

        // volume
        double mp_volume1{ 0.25 };
        pElement1->SetValuesOnIntegrationPoints(MP_VOLUME, { mp_volume1 }, r_current_process_info);
        double mp_volume2{ 0.25 };
        pElement2->SetValuesOnIntegrationPoints(MP_VOLUME, { mp_volume2 }, r_current_process_info);
        double mp_volume3{ 0.25 };
        pElement3->SetValuesOnIntegrationPoints(MP_VOLUME, { mp_volume3 }, r_current_process_info);
        double mp_volume4{ 0.25 };
        pElement4->SetValuesOnIntegrationPoints(MP_VOLUME, { mp_volume4 }, r_current_process_info);
        
        // mp velocity
        array_1d<double, 3> velocity1{1.211324865,0.0, 0.0};
        array_1d<double, 3> velocity2{1.788675135,0.0, 0.0};
        array_1d<double, 3> velocity3{1.788675135,0.0, 0.0};
        array_1d<double, 3> velocity4{1.211324865,0.0, 0.0};
        pElement1->SetValuesOnIntegrationPoints(MP_VELOCITY, { velocity1 }, r_current_process_info);
        pElement2->SetValuesOnIntegrationPoints(MP_VELOCITY, { velocity2 }, r_current_process_info);
        pElement3->SetValuesOnIntegrationPoints(MP_VELOCITY, { velocity3 }, r_current_process_info);
        pElement4->SetValuesOnIntegrationPoints(MP_VELOCITY, { velocity4 }, r_current_process_info);

        // mp acceleration
        array_1d<double, 3> acceleration1{2.211324865, 0.0, 0.0};
        array_1d<double, 3> acceleration2{2.788675135, 0.0, 0.0};
        array_1d<double, 3> acceleration3{2.788675135, 0.0, 0.0};
        array_1d<double, 3> acceleration4{2.211324865, 0.0, 0.0};
        pElement1->SetValuesOnIntegrationPoints(MP_ACCELERATION, { acceleration1 }, r_current_process_info);
        pElement2->SetValuesOnIntegrationPoints(MP_ACCELERATION, { acceleration2 }, r_current_process_info);
        pElement3->SetValuesOnIntegrationPoints(MP_ACCELERATION, { acceleration3 }, r_current_process_info);
        pElement4->SetValuesOnIntegrationPoints(MP_ACCELERATION, { acceleration4 }, r_current_process_info);
    }

    void SetVelocityAndAccelerationToCoordinate(ModelPart& rModelPart)
    {
        ProcessInfo& rProcessInfo = rModelPart.GetProcessInfo();
        std::vector<array_1d<double, 3>> coordinate;
        for (auto& element : rModelPart.Elements())
        {
            element.CalculateOnIntegrationPoints(MP_COORD       , coordinate, rProcessInfo);
            for (auto& component : coordinate[0])
            {
                component = component + 1.0;
            }
            element.SetValuesOnIntegrationPoints(MP_VELOCITY    , coordinate, rProcessInfo);
            for (auto& component : coordinate[0])
            {
                component = -(component + 1.0);
            }
            element.SetValuesOnIntegrationPoints(MP_ACCELERATION, coordinate, rProcessInfo);
        }
    }
    void SetUniformMPMass(ModelPart& rModelPart, double mass)
    {
        ProcessInfo& rProcessInfo = rModelPart.GetProcessInfo();
        for (auto& element : rModelPart.Elements())
        {
            element.SetValuesOnIntegrationPoints(MP_MASS, {mass}, rProcessInfo);
        }
    }
    void SetUniformMPVolume(ModelPart& rModelPart, double volume)
    {
        ProcessInfo& rProcessInfo = rModelPart.GetProcessInfo();
        for (auto& element : rModelPart.Elements())
        {
            element.SetValuesOnIntegrationPoints(MP_MASS, {volume}, rProcessInfo);
        }
    }
    void Prepare2D9EBackgroundModelPart(ModelPart& rBackgroundModelPart)
    {
        // Grid scheme:
        //  13---14---15---16
        //  |    |    |    |
        //  | 7  | 8  | 9  |
        //  |    |    |    |
        //  9----10---11---12 
        //  |    |    |    |
        //  | 4  | 5  | 6  |
        //  |    |    |    |
        //  5----6----7----8 
        //  |    |    |    |
        //  | 1  | 2  | 3  |
        //  |    |    |    |
        //  1----2----3----4

        // Nodes
        auto p_node_1  = rBackgroundModelPart.CreateNewNode( 1,  0.0 ,  0.0 , 0.0);
        auto p_node_2  = rBackgroundModelPart.CreateNewNode( 2,  1.0 ,  0.0 , 0.0);
        auto p_node_3  = rBackgroundModelPart.CreateNewNode( 3,  2.0 ,  0.0 , 0.0);
        auto p_node_4  = rBackgroundModelPart.CreateNewNode( 4,  3.0 ,  0.0 , 0.0);

        auto p_node_5  = rBackgroundModelPart.CreateNewNode( 5,  0.0 ,  1.0 , 0.0);
        auto p_node_6  = rBackgroundModelPart.CreateNewNode( 6,  1.0 ,  1.0 , 0.0);
        auto p_node_7  = rBackgroundModelPart.CreateNewNode( 7,  2.0 ,  1.0 , 0.0);
        auto p_node_8  = rBackgroundModelPart.CreateNewNode( 8,  3.0 ,  1.0 , 0.0);

        auto p_node_9  = rBackgroundModelPart.CreateNewNode( 9 ,  0.0 ,  2.0 , 0.0);
        auto p_node_10 = rBackgroundModelPart.CreateNewNode( 10,  1.0 ,  2.0 , 0.0);
        auto p_node_11 = rBackgroundModelPart.CreateNewNode( 11,  2.0 ,  2.0 , 0.0);
        auto p_node_12 = rBackgroundModelPart.CreateNewNode( 12,  3.0 ,  2.0 , 0.0);

        auto p_node_13 = rBackgroundModelPart.CreateNewNode( 13,  0.0 ,  3.0 , 0.0);
        auto p_node_14 = rBackgroundModelPart.CreateNewNode( 14,  1.0 ,  3.0 , 0.0);
        auto p_node_15 = rBackgroundModelPart.CreateNewNode( 15,  2.0 ,  3.0 , 0.0);
        auto p_node_16 = rBackgroundModelPart.CreateNewNode( 16,  3.0 ,  3.0 , 0.0);

        // Grid Elements
        rBackgroundModelPart.CreateNewElement("Element2D4N", 1, {  1,  2,  6,  5 }, nullptr);
        rBackgroundModelPart.CreateNewElement("Element2D4N", 2, {  2,  3,  7,  6 }, nullptr);
        rBackgroundModelPart.CreateNewElement("Element2D4N", 3, {  3,  4,  8,  7 }, nullptr);
        rBackgroundModelPart.CreateNewElement("Element2D4N", 4, {  5,  6, 10,  9 }, nullptr);
        rBackgroundModelPart.CreateNewElement("Element2D4N", 5, {  6,  7, 11, 10 }, nullptr);
        rBackgroundModelPart.CreateNewElement("Element2D4N", 6, {  7,  8, 12, 11 }, nullptr);
        rBackgroundModelPart.CreateNewElement("Element2D4N", 7, {  9, 10, 14, 13 }, nullptr);
        rBackgroundModelPart.CreateNewElement("Element2D4N", 8, { 10, 11, 15, 14 }, nullptr);
        rBackgroundModelPart.CreateNewElement("Element2D4N", 9, { 11, 12, 16, 15 }, nullptr);
    }

    void Prepare2D9EModelPart(ModelPart& rModelPart, ModelPart& rBackgroundModelPart)
    {
        // Properties
        Properties::Pointer p_elem_prop = rModelPart.CreateNewProperties(0);

        // MP Coordinates
        // Grid Element 1
        array_1d<double, 3> mp_coordinate1{0.211324865,0.211324865, 0.0};
        array_1d<double, 3> mp_coordinate2{0.788675135,0.211324865, 0.0};
        array_1d<double, 3> mp_coordinate3{0.788675135,0.788675135, 0.0};
        array_1d<double, 3> mp_coordinate4{0.211324865,0.788675135, 0.0};
        // Grid Element 2
        array_1d<double, 3> mp_coordinate5{1.211324865,0.211324865, 0.0};
        array_1d<double, 3> mp_coordinate6{1.788675135,0.211324865, 0.0};
        array_1d<double, 3> mp_coordinate7{1.788675135,0.788675135, 0.0};
        array_1d<double, 3> mp_coordinate8{1.211324865,0.788675135, 0.0};
        // Grid Element 3
        array_1d<double, 3> mp_coordinate9{2.211324865,0.211324865, 0.0};
        array_1d<double, 3> mp_coordinate10{2.788675135,0.211324865, 0.0};
        array_1d<double, 3> mp_coordinate11{2.788675135,0.788675135, 0.0};
        array_1d<double, 3> mp_coordinate12{2.211324865,0.788675135, 0.0};
        // Grid Element 4
        array_1d<double, 3> mp_coordinate13{0.211324865,1.211324865, 0.0};
        array_1d<double, 3> mp_coordinate14{0.788675135,1.211324865, 0.0};
        array_1d<double, 3> mp_coordinate15{0.788675135,1.788675135, 0.0};
        array_1d<double, 3> mp_coordinate16{0.211324865,1.788675135, 0.0};
        // Grid Element 5
        array_1d<double, 3> mp_coordinate17{1.211324865,1.211324865, 0.0};
        array_1d<double, 3> mp_coordinate18{1.788675135,1.211324865, 0.0};
        array_1d<double, 3> mp_coordinate19{1.788675135,1.788675135, 0.0};
        array_1d<double, 3> mp_coordinate20{1.211324865,1.788675135, 0.0};
        // Grid Element 6
        array_1d<double, 3> mp_coordinate21{2.211324865,1.211324865, 0.0};
        array_1d<double, 3> mp_coordinate22{2.788675135,1.211324865, 0.0};
        array_1d<double, 3> mp_coordinate23{2.788675135,1.788675135, 0.0};
        array_1d<double, 3> mp_coordinate24{2.211324865,1.788675135, 0.0};
        // Grid Element 7
        array_1d<double, 3> mp_coordinate25{0.211324865,2.211324865, 0.0};
        array_1d<double, 3> mp_coordinate26{0.788675135,2.211324865, 0.0};
        array_1d<double, 3> mp_coordinate27{0.788675135,2.788675135, 0.0};
        array_1d<double, 3> mp_coordinate28{0.211324865,2.788675135, 0.0};
        // Grid Element 8
        array_1d<double, 3> mp_coordinate29{1.211324865,2.211324865, 0.0};
        array_1d<double, 3> mp_coordinate30{1.788675135,2.211324865, 0.0};
        array_1d<double, 3> mp_coordinate31{1.788675135,2.788675135, 0.0};
        array_1d<double, 3> mp_coordinate32{1.211324865,2.788675135, 0.0};
        // Grid Element 9
        array_1d<double, 3> mp_coordinate33{2.211324865,2.211324865, 0.0};
        array_1d<double, 3> mp_coordinate34{2.788675135,2.211324865, 0.0};
        array_1d<double, 3> mp_coordinate35{2.788675135,2.788675135, 0.0};
        array_1d<double, 3> mp_coordinate36{2.211324865,2.788675135, 0.0};


        std::vector<array_1d<double, 3>> mp_coordinates = {mp_coordinate1 , mp_coordinate2 , mp_coordinate3 , mp_coordinate4 , // Grid Element 1
                                                           mp_coordinate5 , mp_coordinate6 , mp_coordinate7 , mp_coordinate8 , // Grid Element 2
                                                           mp_coordinate9 , mp_coordinate10, mp_coordinate11, mp_coordinate12, // Grid Element 3
                                                           mp_coordinate13, mp_coordinate14, mp_coordinate15, mp_coordinate16, // Grid Element 4
                                                           mp_coordinate17, mp_coordinate18, mp_coordinate19, mp_coordinate20, // Grid Element 5
                                                           mp_coordinate21, mp_coordinate22, mp_coordinate23, mp_coordinate24, // Grid Element 6
                                                           mp_coordinate25, mp_coordinate26, mp_coordinate27, mp_coordinate28, // Grid Element 7
                                                           mp_coordinate29, mp_coordinate30, mp_coordinate31, mp_coordinate32, // Grid Element 8
                                                           mp_coordinate29, mp_coordinate30, mp_coordinate31, mp_coordinate32};// Grid Element 9

        PrepareMP<2>(rModelPart, rBackgroundModelPart, mp_coordinates, 0.25);

        auto pElement1 = rModelPart.pGetElement(1);
        auto pElement2 = rModelPart.pGetElement(2);
        auto pElement3 = rModelPart.pGetElement(3);
        auto pElement4 = rModelPart.pGetElement(4);


        // mass
        SetUniformMPVolume(rModelPart, 2.5);

        // volume
        SetUniformMPVolume(rModelPart, 0.25);

        // MP Velocity is set to be coordinate + 1.0, MP Acceleration is set to be -(coordinate + 2.0)
        SetVelocityAndAccelerationToCoordinate(rModelPart);
    }

    /**
    * 
    */
    KRATOS_TEST_CASE_IN_SUITE(MPMFlipParticleMappingUtilityOneGridElement2D, KratosMPMFastSuite)
    {
        KRATOS_TRY;
        const unsigned int dimension = 2;
        Model current_model;
        ModelPart& r_mpm_model_part = current_model.CreateModelPart("MPMModelPart");
        ModelPart& r_background_model_part = current_model.CreateModelPart("MPMBackgroundModelPart");

        r_background_model_part.SetBufferSize(2);

        r_background_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_background_model_part.AddNodalSolutionStepVariable(REACTION);
        r_background_model_part.AddNodalSolutionStepVariable(PRESSURE);
        r_background_model_part.AddNodalSolutionStepVariable(POSITIVE_FACE_PRESSURE);
        r_background_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
        r_background_model_part.AddNodalSolutionStepVariable(NODAL_MASS);
        r_background_model_part.AddNodalSolutionStepVariable(NODAL_MOMENTUM);
        r_background_model_part.AddNodalSolutionStepVariable(NODAL_INERTIA);
        r_background_model_part.AddNodalSolutionStepVariable(VELOCITY);
        r_background_model_part.AddNodalSolutionStepVariable(ACCELERATION);

        Prepare2D1EBackgroundModelPart(r_background_model_part);
        r_mpm_model_part.SetNodes(r_background_model_part.pNodes());

        Prepare2D1EModelPart(r_mpm_model_part, r_background_model_part); // ------------------------------------------------------------

        MPMSearchElementUtility::SearchElement<dimension>(
            r_background_model_part, r_mpm_model_part, 1000, 1e-6);
        
        // const ProcessInfo& r_current_process_info = r_mpm_model_part.GetProcessInfo();
        

        // Run Flip P2G Mapping
        unsigned int echo_level = 0;
        MPMFlipParticleMappingUtility flip_mapping(r_mpm_model_part, r_background_model_part, echo_level);
        flip_mapping.RunP2GMapping();
        

        // Checking values at the nodes
        auto& r_node_1 = r_mpm_model_part.GetNode(1);
        auto& r_node_2 = r_mpm_model_part.GetNode(2);
        auto& r_node_3 = r_mpm_model_part.GetNode(3);
        auto& r_node_4 = r_mpm_model_part.GetNode(4);
        // Check mapped mass
        double& r_node_1_mass = r_node_1.FastGetSolutionStepValue(NODAL_MASS);
        double& r_node_2_mass = r_node_2.FastGetSolutionStepValue(NODAL_MASS);
        double& r_node_3_mass = r_node_3.FastGetSolutionStepValue(NODAL_MASS);
        double& r_node_4_mass = r_node_4.FastGetSolutionStepValue(NODAL_MASS);

        KRATOS_EXPECT_NEAR(r_node_1_mass, 2.5, 1e-6);
        KRATOS_EXPECT_NEAR(r_node_2_mass, 2.5, 1e-6);
        KRATOS_EXPECT_NEAR(r_node_3_mass, 2.5, 1e-6);
        KRATOS_EXPECT_NEAR(r_node_4_mass, 2.5, 1e-6);

        // Check mapped velocity_x
        array_1d<double,3>& r_node_1_velocity_previous = r_node_1.FastGetSolutionStepValue(VELOCITY, 1);
        array_1d<double,3>& r_node_2_velocity_previous = r_node_2.FastGetSolutionStepValue(VELOCITY, 1);
        array_1d<double,3>& r_node_3_velocity_previous = r_node_3.FastGetSolutionStepValue(VELOCITY, 1);
        array_1d<double,3>& r_node_4_velocity_previous = r_node_4.FastGetSolutionStepValue(VELOCITY, 1);

        KRATOS_EXPECT_NEAR(r_node_1_velocity_previous[0], 1.333333333, 1e-6);
        KRATOS_EXPECT_NEAR(r_node_2_velocity_previous[0], 1.666666666, 1e-6);
        KRATOS_EXPECT_NEAR(r_node_3_velocity_previous[0], 1.666666666, 1e-6);
        KRATOS_EXPECT_NEAR(r_node_4_velocity_previous[0], 1.333333333, 1e-6);

        // Check mapped acceleration_x
        array_1d<double,3>& r_node_1_acceleration_previous = r_node_1.FastGetSolutionStepValue(ACCELERATION, 1);
        array_1d<double,3>& r_node_2_acceleration_previous = r_node_2.FastGetSolutionStepValue(ACCELERATION, 1);
        array_1d<double,3>& r_node_3_acceleration_previous = r_node_3.FastGetSolutionStepValue(ACCELERATION, 1);
        array_1d<double,3>& r_node_4_acceleration_previous = r_node_4.FastGetSolutionStepValue(ACCELERATION, 1);

        KRATOS_EXPECT_NEAR(r_node_1_acceleration_previous[0], 2.333333333, 1e-6);
        KRATOS_EXPECT_NEAR(r_node_2_acceleration_previous[0], 2.666666666, 1e-6);
        KRATOS_EXPECT_NEAR(r_node_3_acceleration_previous[0], 2.666666666, 1e-6);
        KRATOS_EXPECT_NEAR(r_node_4_acceleration_previous[0], 2.333333333, 1e-6);

        KRATOS_CATCH("")
    } 
    // end of MPMFlipParticleMappingUtilityOneGridElement2D

    // KRATOS_TEST_CASE_IN_SUITE(MPMFlipParticleMappingUtilityNineGridElement2D, KratosMPMFastSuite)
    // {
    //     KRATOS_TRY;
    //     const unsigned int dimension = 2;
    //     Model current_model;
    //     ModelPart& r_mpm_model_part = current_model.CreateModelPart("MPMModelPart");
    //     ModelPart& r_background_model_part = current_model.CreateModelPart("MPMBackgroundModelPart");

    //     r_background_model_part.SetBufferSize(2);

    //     r_background_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    //     r_background_model_part.AddNodalSolutionStepVariable(REACTION);
    //     r_background_model_part.AddNodalSolutionStepVariable(PRESSURE);
    //     r_background_model_part.AddNodalSolutionStepVariable(POSITIVE_FACE_PRESSURE);
    //     r_background_model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
    //     r_background_model_part.AddNodalSolutionStepVariable(NODAL_MASS);
    //     r_background_model_part.AddNodalSolutionStepVariable(NODAL_MOMENTUM);
    //     r_background_model_part.AddNodalSolutionStepVariable(NODAL_INERTIA);
    //     r_background_model_part.AddNodalSolutionStepVariable(VELOCITY);
    //     r_background_model_part.AddNodalSolutionStepVariable(ACCELERATION);

    //     Prepare2D1EBackgroundModelPart(r_background_model_part);
    //     r_mpm_model_part.SetNodes(r_background_model_part.pNodes());

    //     Prepare2D1EModelPart(r_mpm_model_part, r_background_model_part); // ------------------------------------------------------------

    //     MPMSearchElementUtility::SearchElement<dimension>(r_background_model_part, r_mpm_model_part, 1000, 1e-6);
        
    //     // Run Flip P2G Mapping
    //     MPMFlipParticleMappingUtility flip_mapping(r_mpm_model_part);
    //     flip_mapping.RunP2GMapping();
        

    //     // Checking values at the nodes
    //     auto& r_node_1 = r_mpm_model_part.GetNode(1);
    //     auto& r_node_2 = r_mpm_model_part.GetNode(2);
    //     auto& r_node_3 = r_mpm_model_part.GetNode(3);
    //     auto& r_node_4 = r_mpm_model_part.GetNode(4);
    //     // Check mapped mass
    //     double& r_node_1_mass = r_node_1.FastGetSolutionStepValue(NODAL_MASS);
    //     double& r_node_2_mass = r_node_2.FastGetSolutionStepValue(NODAL_MASS);
    //     double& r_node_3_mass = r_node_3.FastGetSolutionStepValue(NODAL_MASS);
    //     double& r_node_4_mass = r_node_4.FastGetSolutionStepValue(NODAL_MASS);

    //     KRATOS_EXPECT_NEAR(r_node_1_mass, 2.5, 1e-6);
    //     KRATOS_EXPECT_NEAR(r_node_2_mass, 2.5, 1e-6);
    //     KRATOS_EXPECT_NEAR(r_node_3_mass, 2.5, 1e-6);
    //     KRATOS_EXPECT_NEAR(r_node_4_mass, 2.5, 1e-6);

    //     // Check mapped velocity_x
    //     array_1d<double,3>& r_node_1_velocity_previous = r_node_1.FastGetSolutionStepValue(VELOCITY, 1);
    //     array_1d<double,3>& r_node_2_velocity_previous = r_node_2.FastGetSolutionStepValue(VELOCITY, 1);
    //     array_1d<double,3>& r_node_3_velocity_previous = r_node_3.FastGetSolutionStepValue(VELOCITY, 1);
    //     array_1d<double,3>& r_node_4_velocity_previous = r_node_4.FastGetSolutionStepValue(VELOCITY, 1);

    //     KRATOS_EXPECT_NEAR(r_node_1_velocity_previous[0], 1.333333333, 1e-6);
    //     KRATOS_EXPECT_NEAR(r_node_2_velocity_previous[0], 1.666666666, 1e-6);
    //     KRATOS_EXPECT_NEAR(r_node_3_velocity_previous[0], 1.666666666, 1e-6);
    //     KRATOS_EXPECT_NEAR(r_node_4_velocity_previous[0], 1.333333333, 1e-6);

    //     // Check mapped acceleration_x
    //     array_1d<double,3>& r_node_1_acceleration_previous = r_node_1.FastGetSolutionStepValue(ACCELERATION, 1);
    //     array_1d<double,3>& r_node_2_acceleration_previous = r_node_2.FastGetSolutionStepValue(ACCELERATION, 1);
    //     array_1d<double,3>& r_node_3_acceleration_previous = r_node_3.FastGetSolutionStepValue(ACCELERATION, 1);
    //     array_1d<double,3>& r_node_4_acceleration_previous = r_node_4.FastGetSolutionStepValue(ACCELERATION, 1);

    //     KRATOS_EXPECT_NEAR(r_node_1_acceleration_previous[0], 2.333333333, 1e-6);
    //     KRATOS_EXPECT_NEAR(r_node_2_acceleration_previous[0], 2.666666666, 1e-6);
    //     KRATOS_EXPECT_NEAR(r_node_3_acceleration_previous[0], 2.666666666, 1e-6);
    //     KRATOS_EXPECT_NEAR(r_node_4_acceleration_previous[0], 2.333333333, 1e-6);

    //     KRATOS_CATCH("")
    // }
    // void CheckVelocity(ModelPart& rMPMModelPart)
    // {
    //     for (auto iter = rMPMModelPart.NodesBegin(); iter != rMPMModelPart.NodesEnd(); ++iter )
    //     {
    //         const auto& node = rMPMModelPart.Nodes()[iter];

    //         KRATOS_EXPECT_NEAR(node.FastGetSolutionStepValue(VELOCITY, 1), 1.333333333, 1e-6);
    //         rMPMModelPart
    //     }
    // }

} // namespace Kratos::Testing
