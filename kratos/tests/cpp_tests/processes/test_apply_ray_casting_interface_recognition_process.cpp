//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Daniel Diez
//

// Project includes
#include "testing/testing.h"
#include "processes/apply_ray_casting_interface_recognition_process.h"

namespace Kratos::Testing {
    KRATOS_TEST_CASE_IN_SUITE(ApplyRayCastingInterfaceRecognitionProcess3D, KratosCoreFastSuite)
    {
        Model current_model;
        ModelPart& r_volume_model_part = current_model.CreateModelPart("volume");
        ModelPart& r_skin_model_part = current_model.CreateModelPart("skin");
        auto p_prop_0 = r_volume_model_part.CreateNewProperties(0);
        r_volume_model_part.CreateNewNode(1, 0.0,0.0,0.0);
        r_volume_model_part.CreateNewNode(2, 1.0,0.0,0.0);
        r_volume_model_part.CreateNewNode(3, 0.0,1.0,0.0);
        r_volume_model_part.CreateNewNode(4, 0.0,0.0,1.0);
        r_volume_model_part.CreateNewNode(5, 0.0,0.0,-1.0);
        r_volume_model_part.CreateNewElement("Element3D4N",1,{1,2,3,4},p_prop_0);
        r_volume_model_part.CreateNewElement("Element3D4N",2,{1,4,5,2},p_prop_0);

        r_skin_model_part.CreateNewNode(1,-5.0,5.0,0.0);
        r_skin_model_part.CreateNewNode(2,-5.0,-5.0,0.0);
        r_skin_model_part.CreateNewNode(3, 5.0,0.0,0.0);
        r_skin_model_part.CreateNewNode(4,-5.0,0.0,100.0);
        auto p_prop_skin_0 = r_skin_model_part.CreateNewProperties(0);
        r_skin_model_part.CreateNewCondition("SurfaceCondition3D3N",1,{1,2,3},p_prop_skin_0);
        r_skin_model_part.CreateNewCondition("SurfaceCondition3D3N",2,{1,2,4},p_prop_skin_0);
        r_skin_model_part.CreateNewCondition("SurfaceCondition3D3N",3,{1,4,3},p_prop_skin_0);
        for (auto& r_node: r_volume_model_part.Nodes()) {
            r_node.SetValue(DISTANCE, 1.0);
        }
        const Parameters settings(R"({
            "volume_model_part" : "volume",
            "skin_model_part" : "skin",
            "distance_database" : "nodal_non_historical"
        })");
        ApplyRayCastingInterfaceRecognitionProcess<3>(
            current_model,
            settings).Execute();
        for (auto& r_node : r_volume_model_part.Nodes()) {
            double theoretical_distance = 0.0;
            if (r_node.Id() == 4) {
                theoretical_distance = -1.0;
            } else if (r_node.Id() == 5) {
                theoretical_distance = 1.0;
            }
            KRATOS_CHECK_NEAR(r_node.GetValue(DISTANCE),theoretical_distance,1e-10);
        }
    }
    KRATOS_TEST_CASE_IN_SUITE(ApplyRayCastingInterfaceRecognitionProcess2D, KratosCoreFastSuite)
    {
        Model current_model;
        ModelPart& r_volume_model_part = current_model.CreateModelPart("volume");
        ModelPart& r_skin_model_part = current_model.CreateModelPart("skin");
        auto p_prop_0 = r_volume_model_part.CreateNewProperties(0);
        r_volume_model_part.CreateNewNode(1, 0.0,0.0,0.0);
        r_volume_model_part.CreateNewNode(2, 1.0,0.0,0.0);
        r_volume_model_part.CreateNewNode(3, 0.0,1.0,0.0);
        r_volume_model_part.CreateNewNode(4, 0.0,-1.0,0.0);
        r_volume_model_part.CreateNewElement("Element2D3N",1,{1,2,3},p_prop_0);
        r_volume_model_part.CreateNewElement("Element2D3N",2,{1,2,4},p_prop_0);

        r_skin_model_part.CreateNewNode(1,-5.0,0,0.0);
        r_skin_model_part.CreateNewNode(2, 5.0,0.0,0.0);
        r_skin_model_part.CreateNewNode(3,-5.0,100.0,0.0);
        auto p_prop_skin_0 = r_skin_model_part.CreateNewProperties(0);
        std::vector<ModelPart::IndexType> cond1{1, 2};
        std::vector<ModelPart::IndexType> cond2{1, 3};
        std::vector<ModelPart::IndexType> cond3{2, 3};
        r_skin_model_part.CreateNewCondition("LineCondition2D2N",1,cond1,p_prop_skin_0);
        r_skin_model_part.CreateNewCondition("LineCondition2D2N",2,cond2,p_prop_skin_0);
        r_skin_model_part.CreateNewCondition("LineCondition2D2N",3,cond3,p_prop_skin_0);
        for (auto& r_node: r_volume_model_part.Nodes()) {
            r_node.SetValue(DISTANCE, 1.0);
        }
        const Parameters settings(R"({
            "volume_model_part" : "volume",
            "skin_model_part" : "skin",
            "distance_database" : "nodal_non_historical"
        })");
        ApplyRayCastingInterfaceRecognitionProcess<2>(
            current_model,
            settings).Execute();
        for (auto& r_node : r_volume_model_part.Nodes()) {
            double theoretical_distance = 0.0;
            if (r_node.Id() == 3) {
                theoretical_distance = -1.0;
            } else if (r_node.Id() == 4) {
                theoretical_distance = 1.0;
            }
            KRATOS_CHECK_NEAR(r_node.GetValue(DISTANCE),theoretical_distance,1e-10);
        }
    }
}