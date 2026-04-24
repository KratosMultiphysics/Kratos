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
#include "includes/debug_helpers.h"
#include "testing/testing.h"
#include "containers/model.h"
#include "mpm_application_variables.h"
#include "custom_utilities/material_point_search_utility.h"
#include "custom_utilities/mapping_utilities/mpm_flip_particle_mapping_utility.hpp"


namespace Kratos::Testing
{
    void CreateEmptyModels(Model& rCurrentModel)
    {
        ModelPart& rGridModelPart = rCurrentModel.CreateModelPart("GridModelPart");
        ModelPart& rMpmModelPart  = rCurrentModel.CreateModelPart("MpmModelPart");
        rMpmModelPart.SetProcessInfo(rGridModelPart.pGetProcessInfo());
    }

    void AddVariablesToModelPart(ModelPart& rMpmModelPart, ModelPart& rGridModelPart)
    {
        rGridModelPart.SetBufferSize(2);
        rGridModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
        rGridModelPart.AddNodalSolutionStepVariable(REACTION);
        rGridModelPart.AddNodalSolutionStepVariable(POSITIVE_FACE_PRESSURE);
        rGridModelPart.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);
        rGridModelPart.AddNodalSolutionStepVariable(NODAL_MASS);
        rGridModelPart.AddNodalSolutionStepVariable(NODAL_MOMENTUM);
        rGridModelPart.AddNodalSolutionStepVariable(NODAL_INERTIA);
        rGridModelPart.AddNodalSolutionStepVariable(VELOCITY);
        rGridModelPart.AddNodalSolutionStepVariable(ACCELERATION);

        if (rGridModelPart.GetProcessInfo().GetValue(IS_MIXED_FORMULATION))
        {
            rGridModelPart.AddNodalSolutionStepVariable(PRESSURE);
            rGridModelPart.AddNodalSolutionStepVariable(NODAL_MPRESSURE);
        }

        rMpmModelPart.SetBufferSize(2);
        rMpmModelPart.SetNodes(rGridModelPart.pNodes());
        rMpmModelPart.SetNodalSolutionStepVariablesList(rGridModelPart.pGetNodalSolutionStepVariablesList());
    }

    template<class TDataType>
    void SetValuesOnNodes(
        ModelPart& rGridModelPart,
        const Variable<TDataType>& rNodalVariableName,
        const int buffer_index,
        const std::vector<TDataType>& rValues)
    {
        IndexType nodal_index = 0;
        for (Node& r_node : rGridModelPart.Nodes())
        {
            if (r_node.SolutionStepsDataHas(rNodalVariableName)){
                r_node.FastGetSolutionStepValue(rNodalVariableName) = rValues[nodal_index];
            }
            ++nodal_index;
        }
    }

    template<class TDataType>
    void SetUniformValueOnNodes(
        ModelPart& rGridModelPart,
        const Variable<TDataType>& rNodalVariableName,
        const int buffer_index,
        TDataType& rValue)
    {
        IndexType nodal_index = 0;
        for (Node& r_node : rGridModelPart.Nodes())
        {
            if (r_node.SolutionStepsDataHas(rNodalVariableName)){
                r_node.FastGetSolutionStepValue(rNodalVariableName) = rValue;
            }
            ++nodal_index;
        }
    }

    template<class TDataType>
    void SetValuesOnMaterialPoints(
        ModelPart& rMpmModelPart,
        const Variable<TDataType>& rMPVariableName,
        const std::vector<TDataType>& rValues)
    {
        const ProcessInfo& rProcessInfo = rMpmModelPart.GetProcessInfo();
        IndexType mp_index = 0;
        for (auto& material_point_i : rMpmModelPart.Elements())
        {
            material_point_i.SetValuesOnIntegrationPoints(rMPVariableName, {rValues[mp_index]}, rProcessInfo);
            ++mp_index;
        }
    }

    template<class TDataType>
    void SetUniformValueOnMaterialPoints(
        ModelPart& rMpmModelPart,
        const Variable<TDataType>& rMPVariableName,
        const TDataType& rValue)
    {
        const ProcessInfo& rProcessInfo = rMpmModelPart.GetProcessInfo();
        IndexType mp_index = 0;
        for (auto& material_point_i : rMpmModelPart.Elements())
        {
            material_point_i.SetValuesOnIntegrationPoints(rMPVariableName, {rValue}, rProcessInfo);
            ++mp_index;
        }
    }

    void Prepare2D1ElementGridModelPart(ModelPart& rGridModelPart)
    {
        // Grid scheme:
        //  4----3
        //  |    |
        //  |    |
        //  |    |
        //  1----2

        // Nodes
        auto p_node_1 = rGridModelPart.CreateNewNode( 1,  0.0 ,  0.0 , 0.0);
        auto p_node_2 = rGridModelPart.CreateNewNode( 2,  1.0 ,  0.0 , 0.0);
        auto p_node_3 = rGridModelPart.CreateNewNode( 3,  1.0 ,  1.0 , 0.0);
        auto p_node_4 = rGridModelPart.CreateNewNode( 4,  0.0 ,  1.0 , 0.0);

        rGridModelPart.CreateNewElement(
            "Element2D4N", 1, { 1, 2, 3, 4 }, nullptr);
    }

    template <SizeType TDimension>
    void CreateMP(
        ModelPart& rMpmModelPart,
        ModelPart& rGridModelPart,
        Properties::Pointer pProperties,
        const array_1d<double, 3>& rMPCoordinate,
        const double& rMPVolume,
        const Element& rElementType)
    {
        // Create new material point element
        unsigned int new_element_id = rMpmModelPart.NumberOfElements() + 1;

        BinBasedFastPointLocator<TDimension> SearchStructure(rGridModelPart);
        SearchStructure.UpdateSearchDatabase();
        typename BinBasedFastPointLocator<TDimension>::ResultContainerType results(100);
        typename BinBasedFastPointLocator<TDimension>::ResultIteratorType result_begin = results.begin();
        Element::Pointer p_grid_element;
        Vector N;

        bool is_found = SearchStructure.FindPointOnMesh(rMPCoordinate, N, p_grid_element, result_begin);
        if (!is_found) KRATOS_WARNING("MPMTestHelper") << "::search failed." << std::endl;

        auto p_new_mp_geometry = CreateQuadraturePointsUtility<Node>::CreateFromCoordinates(
                            p_grid_element->pGetGeometry(), rMPCoordinate, rMPVolume);

        Element::Pointer p_new_element = rElementType.Create(new_element_id, p_new_mp_geometry, pProperties);

        p_new_element->SetValuesOnIntegrationPoints(MP_COORD, {rMPCoordinate}, rMpmModelPart.GetProcessInfo());
        rMpmModelPart.AddElement(p_new_element);
    }
    template <SizeType TDimension>
    void CreateMaterialPointsFromCoordinates(
        ModelPart& rMpmModelPart,
        ModelPart& rGridModelPart,
        const std::vector<array_1d<double, 3>>& rMPCoordinates,
        const double& rMPVolume,
        const Properties::Pointer pProperties
        )
    {
        // Elements
        std::string element_type_name;
        if (rMpmModelPart.GetProcessInfo().GetValue(IS_MIXED_FORMULATION))
            element_type_name = "MPMUpdatedLagrangianUP";
        else
            element_type_name = "MPMUpdatedLagrangian";

        const Element& r_element_type = KratosComponents<Element>::Get(element_type_name);

        for (auto& mp_coordinate : rMPCoordinates)
        {
            CreateMP<TDimension>(rMpmModelPart, rGridModelPart, pProperties, mp_coordinate, rMPVolume, r_element_type);
        }
    }

    void Prepare2D1ElementMpmModelPart(
        ModelPart& rMpmModelPart,
        ModelPart& rGridModelPart,
        const Properties::Pointer pProperties)
    {
        // MP scheme:
        //  4------3
        //  | X  X |
        //  |      |
        //  | X  X |
        //  1------2
        const unsigned int dimension = 2;

        // Elements
        array_1d<double, 3> mp_coordinate1{0.211324865,0.211324865, 0.0};
        array_1d<double, 3> mp_coordinate2{0.788675135,0.211324865, 0.0};
        array_1d<double, 3> mp_coordinate3{0.788675135,0.788675135, 0.0};
        array_1d<double, 3> mp_coordinate4{0.211324865,0.788675135, 0.0};
        std::vector<array_1d<double, 3>> mp_coordinates = {mp_coordinate1, mp_coordinate2, mp_coordinate3, mp_coordinate4};


        CreateMaterialPointsFromCoordinates<dimension>(rMpmModelPart, rGridModelPart, mp_coordinates, 0.25, pProperties);

        // MP Volume
        SetUniformValueOnMaterialPoints(rMpmModelPart, MP_VOLUME, 0.25);

        // MP Density
        SetUniformValueOnMaterialPoints(rMpmModelPart, MP_DENSITY, 1.0);

        // MP Mass
        SetUniformValueOnMaterialPoints(rMpmModelPart, MP_MASS, 0.25);

        // MP Pressure
        if (rMpmModelPart.GetProcessInfo().GetValue(IS_MIXED_FORMULATION))
            SetUniformValueOnMaterialPoints(rMpmModelPart, MP_PRESSURE, 0.0);

        // Initialize MP Vector Variables to zero
        const array_1d<double, 3> set_initial_mp_vector_variables = ZeroVector(3);

        // MP Displacement
        SetUniformValueOnMaterialPoints(rMpmModelPart, MP_DISPLACEMENT, set_initial_mp_vector_variables);

        // MP Velocity
        SetUniformValueOnMaterialPoints(rMpmModelPart, MP_VELOCITY, set_initial_mp_vector_variables);

        // MP Acceleration
        SetUniformValueOnMaterialPoints(rMpmModelPart, MP_ACCELERATION, set_initial_mp_vector_variables);

        // MP Volume Acceleration
        SetUniformValueOnMaterialPoints(rMpmModelPart, MP_VOLUME_ACCELERATION, set_initial_mp_vector_variables);

        // Search and update shape function values
        MPMSearchElementUtility::SearchElement<dimension>(rGridModelPart, rMpmModelPart, 1000, 1e-6);
    }

    void SetVelocityAndAccelerationToCoordinate(ModelPart& rMpmModelPart)
    {
        ProcessInfo& rProcessInfo = rMpmModelPart.GetProcessInfo();
        std::vector<array_1d<double, 3>> coordinate;
        for (auto& element : rMpmModelPart.Elements())
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

    // void Prepare2D9EGridModelPart(ModelPart& rGridModelPart)
    // {
    //     // Grid scheme:
    //     //  13---14---15---16
    //     //  |    |    |    |
    //     //  | 7  | 8  | 9  |
    //     //  |    |    |    |
    //     //  9----10---11---12
    //     //  |    |    |    |
    //     //  | 4  | 5  | 6  |
    //     //  |    |    |    |
    //     //  5----6----7----8
    //     //  |    |    |    |
    //     //  | 1  | 2  | 3  |
    //     //  |    |    |    |
    //     //  1----2----3----4

    //     // Nodes
    //     auto p_node_1  = rGridModelPart.CreateNewNode( 1,  0.0 ,  0.0 , 0.0);
    //     auto p_node_2  = rGridModelPart.CreateNewNode( 2,  1.0 ,  0.0 , 0.0);
    //     auto p_node_3  = rGridModelPart.CreateNewNode( 3,  2.0 ,  0.0 , 0.0);
    //     auto p_node_4  = rGridModelPart.CreateNewNode( 4,  3.0 ,  0.0 , 0.0);

    //     auto p_node_5  = rGridModelPart.CreateNewNode( 5,  0.0 ,  1.0 , 0.0);
    //     auto p_node_6  = rGridModelPart.CreateNewNode( 6,  1.0 ,  1.0 , 0.0);
    //     auto p_node_7  = rGridModelPart.CreateNewNode( 7,  2.0 ,  1.0 , 0.0);
    //     auto p_node_8  = rGridModelPart.CreateNewNode( 8,  3.0 ,  1.0 , 0.0);

    //     auto p_node_9  = rGridModelPart.CreateNewNode( 9 ,  0.0 ,  2.0 , 0.0);
    //     auto p_node_10 = rGridModelPart.CreateNewNode( 10,  1.0 ,  2.0 , 0.0);
    //     auto p_node_11 = rGridModelPart.CreateNewNode( 11,  2.0 ,  2.0 , 0.0);
    //     auto p_node_12 = rGridModelPart.CreateNewNode( 12,  3.0 ,  2.0 , 0.0);

    //     auto p_node_13 = rGridModelPart.CreateNewNode( 13,  0.0 ,  3.0 , 0.0);
    //     auto p_node_14 = rGridModelPart.CreateNewNode( 14,  1.0 ,  3.0 , 0.0);
    //     auto p_node_15 = rGridModelPart.CreateNewNode( 15,  2.0 ,  3.0 , 0.0);
    //     auto p_node_16 = rGridModelPart.CreateNewNode( 16,  3.0 ,  3.0 , 0.0);

    //     // Grid Elements
    //     rGridModelPart.CreateNewElement("Element2D4N", 1, {  1,  2,  6,  5 }, nullptr);
    //     rGridModelPart.CreateNewElement("Element2D4N", 2, {  2,  3,  7,  6 }, nullptr);
    //     rGridModelPart.CreateNewElement("Element2D4N", 3, {  3,  4,  8,  7 }, nullptr);
    //     rGridModelPart.CreateNewElement("Element2D4N", 4, {  5,  6, 10,  9 }, nullptr);
    //     rGridModelPart.CreateNewElement("Element2D4N", 5, {  6,  7, 11, 10 }, nullptr);
    //     rGridModelPart.CreateNewElement("Element2D4N", 6, {  7,  8, 12, 11 }, nullptr);
    //     rGridModelPart.CreateNewElement("Element2D4N", 7, {  9, 10, 14, 13 }, nullptr);
    //     rGridModelPart.CreateNewElement("Element2D4N", 8, { 10, 11, 15, 14 }, nullptr);
    //     rGridModelPart.CreateNewElement("Element2D4N", 9, { 11, 12, 16, 15 }, nullptr);
    // }

    // void Prepare2D9EModelPart(ModelPart& rMpmModelPart, ModelPart& rGridModelPart)
    // {
    //     // Properties
    //     Properties::Pointer p_elem_prop = rMpmModelPart.CreateNewProperties(0);

    //     // MP Coordinates
    //     // Grid Element 1
    //     array_1d<double, 3> mp_coordinate1{0.211324865,0.211324865, 0.0};
    //     array_1d<double, 3> mp_coordinate2{0.788675135,0.211324865, 0.0};
    //     array_1d<double, 3> mp_coordinate3{0.788675135,0.788675135, 0.0};
    //     array_1d<double, 3> mp_coordinate4{0.211324865,0.788675135, 0.0};
    //     // Grid Element 2
    //     array_1d<double, 3> mp_coordinate5{1.211324865,0.211324865, 0.0};
    //     array_1d<double, 3> mp_coordinate6{1.788675135,0.211324865, 0.0};
    //     array_1d<double, 3> mp_coordinate7{1.788675135,0.788675135, 0.0};
    //     array_1d<double, 3> mp_coordinate8{1.211324865,0.788675135, 0.0};
    //     // Grid Element 3
    //     array_1d<double, 3> mp_coordinate9{2.211324865,0.211324865, 0.0};
    //     array_1d<double, 3> mp_coordinate10{2.788675135,0.211324865, 0.0};
    //     array_1d<double, 3> mp_coordinate11{2.788675135,0.788675135, 0.0};
    //     array_1d<double, 3> mp_coordinate12{2.211324865,0.788675135, 0.0};
    //     // Grid Element 4
    //     array_1d<double, 3> mp_coordinate13{0.211324865,1.211324865, 0.0};
    //     array_1d<double, 3> mp_coordinate14{0.788675135,1.211324865, 0.0};
    //     array_1d<double, 3> mp_coordinate15{0.788675135,1.788675135, 0.0};
    //     array_1d<double, 3> mp_coordinate16{0.211324865,1.788675135, 0.0};
    //     // Grid Element 5
    //     array_1d<double, 3> mp_coordinate17{1.211324865,1.211324865, 0.0};
    //     array_1d<double, 3> mp_coordinate18{1.788675135,1.211324865, 0.0};
    //     array_1d<double, 3> mp_coordinate19{1.788675135,1.788675135, 0.0};
    //     array_1d<double, 3> mp_coordinate20{1.211324865,1.788675135, 0.0};
    //     // Grid Element 6
    //     array_1d<double, 3> mp_coordinate21{2.211324865,1.211324865, 0.0};
    //     array_1d<double, 3> mp_coordinate22{2.788675135,1.211324865, 0.0};
    //     array_1d<double, 3> mp_coordinate23{2.788675135,1.788675135, 0.0};
    //     array_1d<double, 3> mp_coordinate24{2.211324865,1.788675135, 0.0};
    //     // Grid Element 7
    //     array_1d<double, 3> mp_coordinate25{0.211324865,2.211324865, 0.0};
    //     array_1d<double, 3> mp_coordinate26{0.788675135,2.211324865, 0.0};
    //     array_1d<double, 3> mp_coordinate27{0.788675135,2.788675135, 0.0};
    //     array_1d<double, 3> mp_coordinate28{0.211324865,2.788675135, 0.0};
    //     // Grid Element 8
    //     array_1d<double, 3> mp_coordinate29{1.211324865,2.211324865, 0.0};
    //     array_1d<double, 3> mp_coordinate30{1.788675135,2.211324865, 0.0};
    //     array_1d<double, 3> mp_coordinate31{1.788675135,2.788675135, 0.0};
    //     array_1d<double, 3> mp_coordinate32{1.211324865,2.788675135, 0.0};
    //     // Grid Element 9
    //     array_1d<double, 3> mp_coordinate33{2.211324865,2.211324865, 0.0};
    //     array_1d<double, 3> mp_coordinate34{2.788675135,2.211324865, 0.0};
    //     array_1d<double, 3> mp_coordinate35{2.788675135,2.788675135, 0.0};
    //     array_1d<double, 3> mp_coordinate36{2.211324865,2.788675135, 0.0};


    //     std::vector<array_1d<double, 3>> mp_coordinates = {mp_coordinate1 , mp_coordinate2 , mp_coordinate3 , mp_coordinate4 , // Grid Element 1
    //                                                        mp_coordinate5 , mp_coordinate6 , mp_coordinate7 , mp_coordinate8 , // Grid Element 2
    //                                                        mp_coordinate9 , mp_coordinate10, mp_coordinate11, mp_coordinate12, // Grid Element 3
    //                                                        mp_coordinate13, mp_coordinate14, mp_coordinate15, mp_coordinate16, // Grid Element 4
    //                                                        mp_coordinate17, mp_coordinate18, mp_coordinate19, mp_coordinate20, // Grid Element 5
    //                                                        mp_coordinate21, mp_coordinate22, mp_coordinate23, mp_coordinate24, // Grid Element 6
    //                                                        mp_coordinate25, mp_coordinate26, mp_coordinate27, mp_coordinate28, // Grid Element 7
    //                                                        mp_coordinate29, mp_coordinate30, mp_coordinate31, mp_coordinate32, // Grid Element 8
    //                                                        mp_coordinate29, mp_coordinate30, mp_coordinate31, mp_coordinate32};// Grid Element 9

    //     PrepareMP<2>(rMpmModelPart, rGridModelPart, mp_coordinates, 0.25);

    //     // MP Velocity is set to be coordinate + 1.0, MP Acceleration is set to be -(coordinate + 2.0)
    //     SetVelocityAndAccelerationToCoordinate(rMpmModelPart);
    // }

    template<class TDataType>
    void AssertMPVariables(
        const ModelPart& rMpmModelPart,
        const Variable<TDataType>& rMPVariableName,
        const std::vector<TDataType>& rReferenceValues,
        const double tolerance)
    {
        const ProcessInfo& rProcessInfo = rMpmModelPart.GetProcessInfo();
        IndexType mp_index = 0;
        for (auto& material_point_i : rMpmModelPart.Elements())
        {
            std::vector<TDataType> mp_variable_value;
            material_point_i.CalculateOnIntegrationPoints(rMPVariableName, mp_variable_value, rProcessInfo);
            if constexpr(std::is_same_v<TDataType, double>) {
                KRATOS_EXPECT_RELATIVE_NEAR(mp_variable_value[0], rReferenceValues[mp_index], tolerance);
            } else {
                KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(mp_variable_value[0], rReferenceValues[mp_index], tolerance);
            }
            ++mp_index;
        }
    }

    template<class TDataType>
    void AssertNodalVariables(
        const ModelPart& rGridModelPart,
        const Variable<TDataType>& rNodalVariableName,
        const int buffer_index,
        const std::vector<TDataType>& rReferenceValues,
        const double tolerance)
    {
        IndexType nodal_index = 0;
        for (Node& r_node : rGridModelPart.Nodes())
        {
            TDataType nodal_variable_value = r_node.FastGetSolutionStepValue(rNodalVariableName, buffer_index);

            if constexpr(std::is_same_v<TDataType, double>) {
                KRATOS_EXPECT_RELATIVE_NEAR(nodal_variable_value, rReferenceValues[nodal_index], tolerance);
            } else {
                KRATOS_EXPECT_VECTOR_RELATIVE_NEAR(nodal_variable_value, rReferenceValues[nodal_index], tolerance);
            }
            ++nodal_index;
        }
    }
    /**
    *
    */
    KRATOS_TEST_CASE_IN_SUITE(MPMFlipParticleMappingUtilityOneGridElement2D, KratosMPMFastSuite) // #ToDo: To be renamed as MPMBaseParticleMappingUtilityOneGridElement2D
    {
        KRATOS_TRY;
        // --------------------------------------------------------------------------------------- Model Creation --------------------------------------------------------------------------------------- //
        Model current_model;
        CreateEmptyModels(current_model);
        ModelPart& r_grid_model_part = current_model.GetModelPart("GridModelPart");
        ModelPart& r_mpm_model_part  = current_model.GetModelPart("MpmModelPart");

        // Activate Pressure DoF
        r_mpm_model_part.GetProcessInfo().SetValue(IS_MIXED_FORMULATION, true);

        // Properties
        Properties::Pointer p_elem_prop = r_mpm_model_part.CreateNewProperties(0);

        AddVariablesToModelPart(r_mpm_model_part, r_grid_model_part);

        Prepare2D1ElementGridModelPart(r_grid_model_part);
        Prepare2D1ElementMpmModelPart(r_mpm_model_part, r_grid_model_part, p_elem_prop);

        // Set MP Mass
        SetUniformValueOnMaterialPoints(r_mpm_model_part, MP_MASS, 2.5);

        // Set MP Velocity
        const array_1d<double, 3> set_mp_1_velocity{1.211324865,0.0, 0.0};
        const array_1d<double, 3> set_mp_2_velocity{1.788675135,0.0, 0.0};
        const array_1d<double, 3> set_mp_3_velocity{1.788675135,0.0, 0.0};
        const array_1d<double, 3> set_mp_4_velocity{1.211324865,0.0, 0.0};

        const std::vector<array_1d<double,3>> set_mp_velocity_values{set_mp_1_velocity,
                                                                     set_mp_2_velocity,
                                                                     set_mp_3_velocity,
                                                                     set_mp_4_velocity};

        SetValuesOnMaterialPoints(r_mpm_model_part, MP_VELOCITY, set_mp_velocity_values);

        // Set MP Acceleration
        const array_1d<double, 3> set_mp_1_acceleration{2.211324865, 0.0, 0.0};
        const array_1d<double, 3> set_mp_2_acceleration{2.788675135, 0.0, 0.0};
        const array_1d<double, 3> set_mp_3_acceleration{2.788675135, 0.0, 0.0};
        const array_1d<double, 3> set_mp_4_acceleration{2.211324865, 0.0, 0.0};

        const std::vector<array_1d<double,3>> set_mp_acceleration_values{set_mp_1_acceleration,
                                                                         set_mp_2_acceleration,
                                                                         set_mp_3_acceleration,
                                                                         set_mp_4_acceleration};

        SetValuesOnMaterialPoints(r_mpm_model_part, MP_ACCELERATION, set_mp_acceleration_values);

        // Set MP Pressure
        const std::vector<double> set_mp_pressure_values{0.5, 1.0, 1.5, 2.0};
        SetValuesOnMaterialPoints(r_mpm_model_part, MP_PRESSURE, set_mp_pressure_values);

        // ------------------------------------------------------------------------------------------ P2G Test ------------------------------------------------------------------------------------------ //
        // Initialize and run FLIP mapping scheme
        unsigned int echo_level = 0;
        MPMFlipParticleMappingUtility flip_mapping(r_mpm_model_part, r_grid_model_part, echo_level);
        flip_mapping.RunP2GMapping();

        // Checking values at the nodes
        // Check mapped mass
        const std::vector<double> ref_nodal_mass{2.5, 2.5, 2.5, 2.5};
        AssertNodalVariables(r_mpm_model_part, NODAL_MASS, 0, ref_nodal_mass, 1e-6);

        // Check mapped velocity
        const array_1d<double, 3> ref_nodal_velocity_1 { 1.333333333, 0.0, 0.0};
        const array_1d<double, 3> ref_nodal_velocity_2 { 1.666666666, 0.0, 0.0};
        const array_1d<double, 3> ref_nodal_velocity_3 { 1.666666666, 0.0, 0.0};
        const array_1d<double, 3> ref_nodal_velocity_4 { 1.333333333, 0.0, 0.0};

        const std::vector<array_1d<double, 3>> ref_nodal_velocity{ref_nodal_velocity_1,
                                                                  ref_nodal_velocity_2,
                                                                  ref_nodal_velocity_3,
                                                                  ref_nodal_velocity_4};

        AssertNodalVariables(r_mpm_model_part, VELOCITY, 1, ref_nodal_velocity, 1e-6);

        // Check mapped acceleration
        const array_1d<double, 3> ref_nodal_acceleration_1 { 2.333333333, 0.0, 0.0};
        const array_1d<double, 3> ref_nodal_acceleration_2 { 2.666666666, 0.0, 0.0};
        const array_1d<double, 3> ref_nodal_acceleration_3 { 2.666666666, 0.0, 0.0};
        const array_1d<double, 3> ref_nodal_acceleration_4 { 2.333333333, 0.0, 0.0};

        const std::vector<array_1d<double, 3>> ref_nodal_acceleration{ref_nodal_acceleration_1,
                                                                      ref_nodal_acceleration_2,
                                                                      ref_nodal_acceleration_3,
                                                                      ref_nodal_acceleration_4};

        AssertNodalVariables(r_mpm_model_part, ACCELERATION, 1, ref_nodal_acceleration, 1e-6);

        // Check mapped pressure
        const std::vector<double> ref_nodal_pressure{0.877991531432732,
                                                     1.044658198567268,
                                                     1.455341801432732,
                                                     1.622008468567268};

        AssertNodalVariables(r_mpm_model_part, PRESSURE, 1, ref_nodal_pressure, 1e-6);

        // ------------------------------------------------------------------------------------------ G2P Test ------------------------------------------------------------------------------------------ //
        // Setting values for current nodal velocity, acceleration, and pressure to simulate solving
        // Velocity
        const array_1d<double,3> node_1_velocity{ 1.1,  0.1, 0.0};
        const array_1d<double,3> node_2_velocity{ 0.1, -0.1, 0.0};
        const array_1d<double,3> node_3_velocity{-0.1, -0.1, 0.0};
        const array_1d<double,3> node_4_velocity{ 0.1,  1.1, 0.0};

        const std::vector<array_1d<double,3>> node_velocity_values{node_1_velocity,
                                                                   node_2_velocity,
                                                                   node_3_velocity,
                                                                   node_4_velocity};

        SetValuesOnNodes(r_grid_model_part, VELOCITY, 0, node_velocity_values);

        // Acceleration
        const array_1d<double,3> node_1_acceleration{ 1.1,  0.1, 0.0};
        const array_1d<double,3> node_2_acceleration{ 0.1, -0.1, 0.0};
        const array_1d<double,3> node_3_acceleration{-0.1, -0.1, 0.0};
        const array_1d<double,3> node_4_acceleration{ 0.1,  1.1, 0.0};

        const std::vector<array_1d<double,3>> node_acceleration_values{node_1_acceleration,
                                                                       node_2_acceleration,
                                                                       node_3_acceleration,
                                                                       node_4_acceleration};

        SetValuesOnNodes(r_grid_model_part, ACCELERATION, 0, node_acceleration_values);

        // Pressure
        const std::vector<double> node_pressure_values{0.5, 1.0, 1.5, 2.0};
        SetValuesOnNodes(r_grid_model_part, PRESSURE, 0, node_pressure_values);

        // Adding current displacement at grid nodes to be mapped back to MP
        const array_1d<double,3> node_1_displacement{0.0 , 0.0 , 0.0};
        const array_1d<double,3> node_2_displacement{0.15, 0.0 , 0.0};
        const array_1d<double,3> node_3_displacement{0.1 , 0.1 , 0.0};
        const array_1d<double,3> node_4_displacement{0.0 , 0.05, 0.0};
        const std::vector<array_1d<double,3>> node_displacement_values{node_1_displacement,
                                                                       node_2_displacement,
                                                                       node_3_displacement,
                                                                       node_4_displacement};

        SetValuesOnNodes(r_grid_model_part, DISPLACEMENT, 0, node_displacement_values);

        // Adding initial MP displacement
        const array_1d<double, 3> set_mp_1_displacement{0.1, 0.2, 0.0};
        const array_1d<double, 3> set_mp_2_displacement{0.3, 0.4, 0.0};
        const array_1d<double, 3> set_mp_3_displacement{0.5, 0.6, 0.0};
        const array_1d<double, 3> set_mp_4_displacement{0.7, 0.8, 0.0};
        const std::vector<array_1d<double,3>> set_mp_displacement_values{set_mp_1_displacement,
                                                                         set_mp_2_displacement,
                                                                         set_mp_3_displacement,
                                                                         set_mp_4_displacement};

        SetValuesOnMaterialPoints(r_mpm_model_part, MP_DISPLACEMENT, set_mp_displacement_values);

        flip_mapping.RunG2PMapping();

        // Material points G2P Checks

        // Check MP acceleration
        array_1d<double, 3> ref_mp_acceleration_1 { 0.713076828853815,  0.224401693432732, 0.0};
        array_1d<double, 3> ref_mp_acceleration_2 { 0.233333333146185, -0.013076828432732, 0.0};
        array_1d<double, 3> ref_mp_acceleration_3 { 0.020256504853815,  0.108931639432732, 0.0};
        array_1d<double, 3> ref_mp_acceleration_4 { 0.233333333146185,  0.679743495567268, 0.0};

        const std::vector<array_1d<double, 3>> ref_mp_accelerations{ref_mp_acceleration_1,
                                                                    ref_mp_acceleration_2,
                                                                    ref_mp_acceleration_3,
                                                                    ref_mp_acceleration_4};

        AssertMPVariables(r_mpm_model_part, MP_ACCELERATION, ref_mp_accelerations, 1e-6);

        // Check MP displacement
        array_1d<double, 3> ref_mp_displacement_1 {0.129465819821637, 0.212799153178363, 0.0};
        array_1d<double, 3> ref_mp_displacement_2 {0.409967936928363, 0.418899576571637, 0.0};
        array_1d<double, 3> ref_mp_displacement_3 {0.587200846821637, 0.670534180178363, 0.0};
        array_1d<double, 3> ref_mp_displacement_4 {0.723365396428363, 0.847767090071637, 0.0};

        const std::vector<array_1d<double, 3>> ref_mp_displacements{ref_mp_displacement_1,
                                                                    ref_mp_displacement_2,
                                                                    ref_mp_displacement_3,
                                                                    ref_mp_displacement_4};

        AssertMPVariables(r_mpm_model_part, MP_DISPLACEMENT, ref_mp_displacements, 1e-6);

        // Check MP coordinate
        array_1d<double, 3> ref_mp_coordinate_1 {0.240790684821637, 0.224124018178363, 0.0};
        array_1d<double, 3> ref_mp_coordinate_2 {0.898643071928363, 0.230224441571637, 0.0};
        array_1d<double, 3> ref_mp_coordinate_3 {0.875875981821637, 0.859209315178363, 0.0};
        array_1d<double, 3> ref_mp_coordinate_4 {0.234690261428363, 0.836442225071637, 0.0};

        const std::vector<array_1d<double, 3>> ref_mp_coordinates{ref_mp_coordinate_1,
                                                                  ref_mp_coordinate_2,
                                                                  ref_mp_coordinate_3,
                                                                  ref_mp_coordinate_4};

        AssertMPVariables(r_mpm_model_part, MP_COORD, ref_mp_coordinates, 1e-6);

        // Check MP velocity // -------------------------------------------------------------------------- #ToDo: To be moved to FLIP specific test
        // std::vector<array_1d<double, 3>> mp_velocity_1{};
        // std::vector<array_1d<double, 3>> mp_velocity_2{};
        // std::vector<array_1d<double, 3>> mp_velocity_3{};
        // std::vector<array_1d<double, 3>> mp_velocity_4{};

        // r_element_1.CalculateOnIntegrationPoints(MP_VELOCITY, mp_velocity_1, rProcessInfo);
        // r_element_2.CalculateOnIntegrationPoints(MP_VELOCITY, mp_velocity_2, rProcessInfo);
        // r_element_3.CalculateOnIntegrationPoints(MP_VELOCITY, mp_velocity_3, rProcessInfo);
        // r_element_4.CalculateOnIntegrationPoints(MP_VELOCITY, mp_velocity_4, rProcessInfo);

        // array_1d<double, 3> ref_mp_velocity_1 {0.0, 0.0, 0.0};
        // array_1d<double, 3> ref_mp_velocity_2 {0.0, 0.0, 0.0};
        // array_1d<double, 3> ref_mp_velocity_3 {0.0, 0.0, 0.0};
        // array_1d<double, 3> ref_mp_velocity_4 {0.0, 0.0, 0.0};

        // KRATOS_WATCH(mp_velocity_1)
        // KRATOS_WATCH(mp_velocity_2)
        // KRATOS_WATCH(mp_velocity_3)
        // KRATOS_WATCH(mp_velocity_4)

        // Check MP pressure
        const double ref_mp_pressure_1 = 0.877991531432732;
        const double ref_mp_pressure_2 = 1.044658198567268;
        const double ref_mp_pressure_3 = 1.455341801432732;
        const double ref_mp_pressure_4 = 1.622008468567268;

        const std::vector<double> ref_mp_pressures{ref_mp_pressure_1,
                                                   ref_mp_pressure_2,
                                                   ref_mp_pressure_3,
                                                   ref_mp_pressure_4};

        AssertMPVariables(r_mpm_model_part, MP_PRESSURE, ref_mp_pressures, 1e-6);

        KRATOS_CATCH("")
    }
    // end of MPMFlipParticleMappingUtilityOneGridElement2D

    // FAIL: test_execution (mpm_test_factory.GravityTimeStepTableTest.test_execution)
    // reference for 1 thread
    // AssertionError: False is not true : -4.3750899682120523e-07 != -1.0246735978938274e-07, rel_tol = 1e-07, abs_tol = 1e-07 : Error checking material point 26 MP_ACCELERATION results.
    // reference for 2 threads
    // AssertionError: False is not true : -4.276550023090116e-07 != -1.0246735978938274e-07, rel_tol = 1e-07, abs_tol = 1e-07 : Error checking material point 26 MP_ACCELERATION results.
    // reference for 3 threads and above
    // AssertionError: False is not true : -1.0246735978938274e-07 != 5.829681202114759e-23, rel_tol = 1e-07, abs_tol = 1e-07 : Error checking material point 26 MP_ACCELERATION results.


} // namespace Kratos::Testing
