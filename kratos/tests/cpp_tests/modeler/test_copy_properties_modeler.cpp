//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//
//

// Project includes
#include "testing/testing.h"
#include "modeler/connectivity_preserve_modeler.h"
#include "modeler/duplicate_mesh_modeler.h"
#include "modeler/copy_properties_modeler.h"

namespace Kratos
{

namespace Testing
{

    KRATOS_TEST_CASE_IN_SUITE(CopyPropertiesModelerConnectivityPreserve, KratosCoreFastSuite)
    {
        Model model;
        auto& model_part_1 = model.CreateModelPart("origin");
        auto& model_part_2 = model.CreateModelPart("destination");

        model_part_1.CreateNewProperties(0);
        model_part_1.CreateNewProperties(1);
        model_part_1.GetProperties(1).SetValue(DISTANCE, 1.1);

        model_part_2.CreateNewProperties(0);
        model_part_2.CreateNewProperties(5);
        model_part_2.CreateNewProperties(6);

        model_part_1.CreateNewNode(1, 0.0, 0.0, 0.0);
        model_part_1.CreateNewNode(2, 1.0, 0.0, 0.0);
        model_part_1.CreateNewNode(3, 0.0, 1.0, 0.0);
        model_part_1.CreateNewElement("Element2D3N", 1, {1, 2, 3}, model_part_1.pGetProperties(1));

        ConnectivityPreserveModeler().GenerateModelPart(model_part_1, model_part_2,
            KratosComponents<Element>::Get("Element2D3N"), KratosComponents<Condition>::Get("LineCondition2D2N"));

        Parameters parameters(R"({
            "origin_model_part_name"      : "origin",
            "destination_model_part_name" : "destination"
        })");
        CopyPropertiesModeler(model, parameters).SetupModelPart();

        KRATOS_EXPECT_EQ(model_part_1.NumberOfProperties(), 2);
        KRATOS_EXPECT_EQ(model_part_2.NumberOfProperties(), 2);
        KRATOS_EXPECT_EQ(model_part_2.GetElement(1).GetProperties().Id(), 1);

        model_part_2.GetProperties(1).SetValue(DISTANCE, 2.2);
        KRATOS_EXPECT_EQ(model_part_1.GetProperties(1).GetValue(DISTANCE), 1.1);
        KRATOS_EXPECT_EQ(model_part_2.GetProperties(1).GetValue(DISTANCE), 2.2);
    }

    KRATOS_TEST_CASE_IN_SUITE(CopyPropertiesModelerDuplicateMesh, KratosCoreFastSuite)
    {
        Model model;
        auto& model_part_1 = model.CreateModelPart("origin");
        auto& model_part_2 = model.CreateModelPart("destination");

        model_part_1.CreateNewProperties(0);
        model_part_1.CreateNewProperties(1);
        model_part_1.GetProperties(1).SetValue(DISTANCE, 1.1);

        model_part_2.CreateNewProperties(0);
        model_part_2.CreateNewProperties(5);
        model_part_2.CreateNewProperties(6);

        model_part_1.CreateNewNode(1, 0.0, 0.0, 0.0);
        model_part_1.CreateNewNode(2, 1.0, 0.0, 0.0);
        model_part_1.CreateNewNode(3, 0.0, 1.0, 0.0);
        model_part_1.CreateNewElement("Element2D3N", 1, {1, 2, 3}, model_part_1.pGetProperties(1));

        DuplicateMeshModeler(model_part_1).GenerateMesh(model_part_2,
            KratosComponents<Element>::Get("Element2D3N"), KratosComponents<Condition>::Get("LineCondition2D2N"));

        CopyPropertiesModeler(model_part_1, model_part_2).SetupModelPart();

        KRATOS_EXPECT_EQ(model_part_1.NumberOfProperties(), 2);
        KRATOS_EXPECT_EQ(model_part_2.NumberOfProperties(), 2);
        KRATOS_EXPECT_EQ(model_part_2.GetElement(1).GetProperties().Id(), 1);

        model_part_2.GetProperties(1).SetValue(DISTANCE, 2.2);
        KRATOS_EXPECT_EQ(model_part_1.GetProperties(1).GetValue(DISTANCE), 1.1);
        KRATOS_EXPECT_EQ(model_part_2.GetProperties(1).GetValue(DISTANCE), 2.2);
    }

}

}  // namespace Kratos.
