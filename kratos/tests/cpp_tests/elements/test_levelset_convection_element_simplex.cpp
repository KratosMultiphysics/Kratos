//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Daniel Diez
//

// System includes
#include <limits>

// External includes

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "elements/levelset_convection_element_simplex.h"

// Utility includes
#include "utilities/geometry_utilities.h"

namespace Kratos
{
namespace Testing
{
    KRATOS_TEST_CASE_IN_SUITE(LevelSetConvectionElement2D, KratosCoreFastSuite) {
        Model current_model;
        ModelPart& model_part = current_model.CreateModelPart("Main");
        model_part.SetBufferSize(3);

        // Variables addition
        model_part.AddNodalSolutionStepVariable(VELOCITY);
        model_part.AddNodalSolutionStepVariable(DISTANCE);

        // Process info creation
        double delta_time = 0.1;
        double cross_wind_diff = 0.5;
        model_part.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.001);
        model_part.GetProcessInfo().SetValue(DELTA_TIME, delta_time);

        ConvectionDiffusionSettings::Pointer p_conv_diff_settings = Kratos::make_unique<ConvectionDiffusionSettings>();
        model_part.GetProcessInfo().SetValue(CONVECTION_DIFFUSION_SETTINGS, p_conv_diff_settings);
        model_part.GetProcessInfo().SetValue(CROSS_WIND_STABILIZATION_FACTOR, cross_wind_diff);
        p_conv_diff_settings->SetUnknownVariable(DISTANCE);
        p_conv_diff_settings->SetConvectionVariable(VELOCITY);

        // Set the element properties
        Properties::Pointer pElemProp = model_part.CreateNewProperties(0);

        // Geometry creation
        model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
        model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
        std::vector<ModelPart::IndexType> elemNodes {1, 2, 3};
        auto p_element = model_part.CreateNewElement("LevelSetConvectionElementSimplex2D3N", 1, elemNodes, pElemProp);
        
        // Define the nodal values
        Matrix vel_original(3,2);
        vel_original(0,0) = 0.0; vel_original(0,1) = 0.1;
        vel_original(1,0) = 0.1; vel_original(1,1) = 0.2;
        vel_original(2,0) = 0.2; vel_original(2,1) = 0.3;

        for(unsigned int i=0; i<3; i++){
            for(unsigned int k=0; k<2; k++){
                p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[k]    = vel_original(i,k);
                p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.9*vel_original(i,k);
                p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = 0.75*vel_original(i,k);
            }
        }
        p_element->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -1.0;
        p_element->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = -1.0;
        p_element->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) =  0.5;

        p_element->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE,1) = -0.8;
        p_element->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE,1) = -0.8;
        p_element->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE,1) =  0.7;

        // Compute RHS and LHS
        Vector rhs = ZeroVector(4);
        Matrix lhs = ZeroMatrix(4,4);

        p_element->Initialize();
        const ProcessInfo& const_process_info = model_part.GetProcessInfo();
        p_element->CalculateLocalSystem(lhs, rhs, const_process_info);

        KRATOS_CHECK_EQUAL(rhs.size(), 3);
        KRATOS_CHECK_EQUAL(lhs.size1(), 3);
        KRATOS_CHECK_EQUAL(lhs.size2(), 3);
        KRATOS_CHECK_NEAR(rhs(0), 0.1186006876 , 1e-7);
        KRATOS_CHECK_NEAR(rhs(1), 0.4123679368 , 1e-7);
        KRATOS_CHECK_NEAR(rhs(2), 0.3265313756 , 1e-7);
        KRATOS_CHECK_NEAR(lhs(0,0), 0.367583692212 , 1e-7);
        KRATOS_CHECK_NEAR(lhs(0,1), -0.301736964939 , 1e-7);
        KRATOS_CHECK_NEAR(lhs(0,2), -0.336905789443 , 1e-7);
        KRATOS_CHECK_NEAR(lhs(1,0), 0.493383924505 , 1e-7);
        KRATOS_CHECK_NEAR(lhs(1,1), 1.1281300063 , 1e-7);
        KRATOS_CHECK_NEAR(lhs(1,2), 0.677825753415 , 1e-7);
        KRATOS_CHECK_NEAR(lhs(2,0), 0.73444904995 , 1e-7);
        KRATOS_CHECK_NEAR(lhs(2,1), 0.864023625302 , 1e-7);
        KRATOS_CHECK_NEAR(lhs(2,2), 1.37324670269 , 1e-7);        
    }

    KRATOS_TEST_CASE_IN_SUITE(LevelSetConvectionElement2DImplicit, KratosCoreFastSuite) {
        Model current_model;
        ModelPart& model_part = current_model.CreateModelPart("Main");
        model_part.SetBufferSize(3);

        // Variables addition
        model_part.AddNodalSolutionStepVariable(VELOCITY);
        model_part.AddNodalSolutionStepVariable(DISTANCE);

        // Process info creation
        double delta_time = 0.1;
        double cross_wind_diff = 0.5;
        model_part.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.001);
        model_part.GetProcessInfo().SetValue(DELTA_TIME, delta_time);

        ConvectionDiffusionSettings::Pointer p_conv_diff_settings = Kratos::make_unique<ConvectionDiffusionSettings>();
        model_part.GetProcessInfo().SetValue(CONVECTION_DIFFUSION_SETTINGS, p_conv_diff_settings);
        model_part.GetProcessInfo().SetValue(CROSS_WIND_STABILIZATION_FACTOR, cross_wind_diff);
        model_part.GetProcessInfo().SetValue(TIME_INTEGRATION_THETA, 1.0);
        p_conv_diff_settings->SetUnknownVariable(DISTANCE);
        p_conv_diff_settings->SetConvectionVariable(VELOCITY);

        // Set the element properties
        Properties::Pointer pElemProp = model_part.CreateNewProperties(0);

        // Geometry creation
        model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
        model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
        std::vector<ModelPart::IndexType> elemNodes {1, 2, 3};
        auto p_element = model_part.CreateNewElement("LevelSetConvectionElementSimplex2D3N", 1, elemNodes, pElemProp);
        
        // Define the nodal values
        Matrix vel_original(3,2);
        vel_original(0,0) = 0.0; vel_original(0,1) = 0.1;
        vel_original(1,0) = 0.1; vel_original(1,1) = 0.2;
        vel_original(2,0) = 0.2; vel_original(2,1) = 0.3;

        for(unsigned int i=0; i<3; i++){
            for(unsigned int k=0; k<2; k++){
                p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[k]    = vel_original(i,k);
                p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.9*vel_original(i,k);
                p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = 0.75*vel_original(i,k);
            }
        }
        p_element->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -1.0;
        p_element->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = -1.0;
        p_element->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) =  0.5;

        p_element->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE,1) = -0.8;
        p_element->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE,1) = -0.8;
        p_element->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE,1) =  0.7;

        // Compute RHS and LHS
        Vector rhs = ZeroVector(4);
        Matrix lhs = ZeroMatrix(4,4);

        p_element->Initialize();
        const ProcessInfo& const_process_info = model_part.GetProcessInfo();
        p_element->CalculateLocalSystem(lhs, rhs, const_process_info);
        
        KRATOS_CHECK_EQUAL(rhs.size(), 3);
        KRATOS_CHECK_EQUAL(lhs.size1(), 3);
        KRATOS_CHECK_EQUAL(lhs.size2(), 3);
        KRATOS_CHECK_NEAR(rhs(0), 0.1186006876 , 1e-7);
        KRATOS_CHECK_NEAR(rhs(1), 0.4123679368 , 1e-7);
        KRATOS_CHECK_NEAR(rhs(2), 0.3265313756 , 1e-7);
        KRATOS_CHECK_NEAR(lhs(0,0), 0.490708692212 , 1e-7);
        KRATOS_CHECK_NEAR(lhs(0,1), -0.367257798272 , 1e-7);
        KRATOS_CHECK_NEAR(lhs(0,2), -0.39450995611 , 1e-7);
        KRATOS_CHECK_NEAR(lhs(1,0), 0.398175591172 , 1e-7);
        KRATOS_CHECK_NEAR(lhs(1,1), 1.2075050063 , 1e-7);
        KRATOS_CHECK_NEAR(lhs(1,2), 0.693659086748 , 1e-7);
        KRATOS_CHECK_NEAR(lhs(2,0), 0.635282383283 , 1e-7);
        KRATOS_CHECK_NEAR(lhs(2,1), 0.873919458635 , 1e-7);
        KRATOS_CHECK_NEAR(lhs(2,2), 1.46251753603 , 1e-7);
    }

    KRATOS_TEST_CASE_IN_SUITE(LevelSetConvectionElement3D, KratosCoreFastSuite) {
        Model current_model;
        ModelPart& model_part = current_model.CreateModelPart("Main");
        model_part.SetBufferSize(3);

        // Variables addition
        model_part.AddNodalSolutionStepVariable(VELOCITY);
        model_part.AddNodalSolutionStepVariable(DISTANCE);

        // Process info creation
        double delta_time = 0.1;
        double cross_wind_diff = 0.5;
        model_part.GetProcessInfo().SetValue(DYNAMIC_TAU, 0.001);
        model_part.GetProcessInfo().SetValue(DELTA_TIME, delta_time);

        ConvectionDiffusionSettings::Pointer p_conv_diff_settings = Kratos::make_unique<ConvectionDiffusionSettings>();
        model_part.GetProcessInfo().SetValue(CONVECTION_DIFFUSION_SETTINGS, p_conv_diff_settings);
        model_part.GetProcessInfo().SetValue(CROSS_WIND_STABILIZATION_FACTOR, cross_wind_diff);
        p_conv_diff_settings->SetUnknownVariable(DISTANCE);
        p_conv_diff_settings->SetConvectionVariable(VELOCITY);

        // Set the element properties
        Properties::Pointer pElemProp = model_part.CreateNewProperties(0);

        // Geometry creation
        model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
        model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
        model_part.CreateNewNode(3, 0.0, 1.0, 0.0);
        model_part.CreateNewNode(4, 0.0, 0.0, 1.0);
        std::vector<ModelPart::IndexType> elemNodes {1, 2, 3, 4};
        auto p_element = model_part.CreateNewElement("LevelSetConvectionElementSimplex3D4N", 1, elemNodes, pElemProp);
        
        // Define the nodal values
        Matrix vel_original(4,3);
        vel_original(0,0) = 0.0; vel_original(0,1) = 0.1; vel_original(0,2) = 0.1;
        vel_original(1,0) = 0.1; vel_original(1,1) = 0.2; vel_original(1,2) = 0.2;
        vel_original(2,0) = 0.2; vel_original(2,1) = 0.3; vel_original(2,2) = 0.3;
        vel_original(3,0) = 0.2; vel_original(3,1) = 0.3; vel_original(3,2) = 0.3;

        for(unsigned int i=0; i<4; i++){
            for(unsigned int k=0; k<3; k++){
                p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[k]    = vel_original(i,k);
                p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 1)[k] = 0.9*vel_original(i,k);
                p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, 2)[k] = 0.75*vel_original(i,k);
            }
        }
        p_element->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) = -1.0;
        p_element->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) = -1.0;
        p_element->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE) =  0.5;
        p_element->GetGeometry()[3].FastGetSolutionStepValue(DISTANCE) =  0.5;

        p_element->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE,1) = -0.8;
        p_element->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE,1) = -0.8;
        p_element->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE,1) =  0.7;
        p_element->GetGeometry()[3].FastGetSolutionStepValue(DISTANCE,1) =  0.7;

        // Compute RHS and LHS
        Vector rhs = ZeroVector(5);
        Matrix lhs = ZeroMatrix(5,4);

        p_element->Initialize();
        const ProcessInfo& const_process_info = model_part.GetProcessInfo();
        p_element->CalculateLocalSystem(lhs, rhs, const_process_info);

        KRATOS_CHECK_EQUAL(rhs.size(), 4);
        KRATOS_CHECK_EQUAL(lhs.size1(), 4);
        KRATOS_CHECK_EQUAL(lhs.size2(), 4);
        KRATOS_CHECK_NEAR(rhs(0), -0.0191355037271 , 1e-7);
        KRATOS_CHECK_NEAR(rhs(1), 0.0913265708816 , 1e-7);
        KRATOS_CHECK_NEAR(rhs(2), 0.0771336330894 , 1e-7);
        KRATOS_CHECK_NEAR(rhs(3), 0.0771336330894 , 1e-7);
        KRATOS_CHECK_NEAR(lhs(0,0), 0.0118217128739 , 1e-7);
        KRATOS_CHECK_NEAR(lhs(0,1), -0.124973472054 , 1e-7);
        KRATOS_CHECK_NEAR(lhs(0,2), -0.132050658822 , 1e-7);
        KRATOS_CHECK_NEAR(lhs(0,3), -0.132050658822 , 1e-7);
        KRATOS_CHECK_NEAR(lhs(1,0), 0.0952497205202 , 1e-7);
        KRATOS_CHECK_NEAR(lhs(1,1), 0.224296426094 , 1e-7);
        KRATOS_CHECK_NEAR(lhs(1,2), 0.133855853867 , 1e-7);
        KRATOS_CHECK_NEAR(lhs(1,3), 0.133855853867 , 1e-7);
        KRATOS_CHECK_NEAR(lhs(2,0), 0.13203719997 , 1e-7);
        KRATOS_CHECK_NEAR(lhs(2,1), 0.16361977298 , 1e-7);
        KRATOS_CHECK_NEAR(lhs(2,2), 0.264675745948 , 1e-7);
        KRATOS_CHECK_NEAR(lhs(2,3), 0.167998225674 , 1e-7);
        KRATOS_CHECK_NEAR(lhs(3,0), 0.13203719997 , 1e-7);
        KRATOS_CHECK_NEAR(lhs(3,1), 0.16361977298 , 1e-7);
        KRATOS_CHECK_NEAR(lhs(3,2), 0.167998225674 , 1e-7);
        KRATOS_CHECK_NEAR(lhs(3,3), 0.264675745948 , 1e-7);
    }

    

} // namespace Testing.
} // namespace Kratos.
