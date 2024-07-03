//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/variables.h"
#include "mpi/utilities/distributed_model_part_initializer.h"

namespace Kratos::Testing {

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(DistributedModelPartInitializer_NoSubModelParts_Empty, KratosMPICoreFastSuite)
{
    Model current_model;
    ModelPart& main_model_part = current_model.CreateModelPart("main");
    main_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

    DistributedModelPartInitializer(main_model_part, Testing::GetDefaultDataCommunicator(), 0).Execute();

    KRATOS_EXPECT_TRUE(main_model_part.IsDistributed());

    KRATOS_EXPECT_EQ(main_model_part.NumberOfSubModelParts(), 0);
    KRATOS_EXPECT_EQ(main_model_part.GetCommunicator().GlobalNumberOfNodes(), 0);
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(DistributedModelPartInitializer_NoSubModelParts_WithNodes, KratosMPICoreFastSuite)
{
    Model current_model;
    ModelPart& main_model_part = current_model.CreateModelPart("main");
    main_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

    const int num_nodes = 10;

    for (int i=0; i<num_nodes; ++i) {
        main_model_part.CreateNewNode(i+1, 0,0,0); // coords don't matter for this test
    }

    DistributedModelPartInitializer(main_model_part, Testing::GetDefaultDataCommunicator(), 0).Execute();

    KRATOS_EXPECT_TRUE(main_model_part.IsDistributed());

    KRATOS_EXPECT_EQ(main_model_part.NumberOfSubModelParts(), 0);
    KRATOS_EXPECT_EQ(main_model_part.GetCommunicator().GlobalNumberOfNodes(), num_nodes);
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(DistributedModelPartInitializer_1_SubModelPart_Empty, KratosMPICoreFastSuite)
{
    Model current_model;
    ModelPart& main_model_part = current_model.CreateModelPart("main");
    main_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

    const DataCommunicator& r_world = Testing::GetDefaultDataCommunicator();

    if (r_world.Rank() == 0) {
        main_model_part.CreateSubModelPart("sub");
    }

    DistributedModelPartInitializer(main_model_part, Testing::GetDefaultDataCommunicator(), 0).Execute();

    KRATOS_EXPECT_TRUE(main_model_part.IsDistributed());

    KRATOS_EXPECT_EQ(main_model_part.NumberOfSubModelParts(), 1);
    KRATOS_EXPECT_TRUE(main_model_part.HasSubModelPart("sub"));
    KRATOS_EXPECT_EQ(main_model_part.GetCommunicator().GlobalNumberOfNodes(), 0);
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(DistributedModelPartInitializer_1_SubModelPart_WithNodes, KratosMPICoreFastSuite)
{
    Model current_model;
    ModelPart& main_model_part = current_model.CreateModelPart("main");
    main_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

    const int num_nodes_main = 10;
    const int num_nodes_sub = 5;

    for (int i=0; i<num_nodes_main; ++i) {
        main_model_part.CreateNewNode(i+1, 0,0,0); // coords don't matter for this test
    }

    const DataCommunicator& r_world = Testing::GetDefaultDataCommunicator();

    if (r_world.Rank() == 0) {
        auto& r_sub_mp = main_model_part.CreateSubModelPart("sub");

        for (int i=0; i<num_nodes_sub; ++i) {
            r_sub_mp.CreateNewNode(i+1+num_nodes_main, 0,0,0); // coords don't matter for this test
        }
    }

    DistributedModelPartInitializer(main_model_part, Testing::GetDefaultDataCommunicator(), 0).Execute();

    KRATOS_EXPECT_TRUE(main_model_part.IsDistributed());

    KRATOS_EXPECT_EQ(main_model_part.NumberOfSubModelParts(), 1);
    KRATOS_EXPECT_TRUE(main_model_part.HasSubModelPart("sub"));
    KRATOS_EXPECT_EQ(main_model_part.GetCommunicator().GlobalNumberOfNodes(), num_nodes_main+num_nodes_sub);
    KRATOS_EXPECT_EQ(main_model_part.GetSubModelPart("sub").GetCommunicator().GlobalNumberOfNodes(), num_nodes_sub);
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(DistributedModelPartInitializer_2_SubModelParts_Empty, KratosMPICoreFastSuite)
{
    Model current_model;
    ModelPart& main_model_part = current_model.CreateModelPart("main");
    main_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

    const DataCommunicator& r_world = Testing::GetDefaultDataCommunicator();

    if (r_world.Rank() == 0) {
        main_model_part.CreateSubModelPart("sub");
        main_model_part.CreateSubModelPart("another_sub");
    }

    DistributedModelPartInitializer(main_model_part, Testing::GetDefaultDataCommunicator(), 0).Execute();

    KRATOS_EXPECT_TRUE(main_model_part.IsDistributed());

    KRATOS_EXPECT_EQ(main_model_part.NumberOfSubModelParts(), 2);
    KRATOS_EXPECT_TRUE(main_model_part.HasSubModelPart("sub"));
    KRATOS_EXPECT_TRUE(main_model_part.HasSubModelPart("another_sub"));
    KRATOS_EXPECT_EQ(main_model_part.GetCommunicator().GlobalNumberOfNodes(), 0);
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(DistributedModelPartInitializer_2_SubSubModelParts_Empty, KratosMPICoreFastSuite)
{
    Model current_model;
    ModelPart& main_model_part = current_model.CreateModelPart("main");
    main_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

    const DataCommunicator& r_world = Testing::GetDefaultDataCommunicator();

    if (r_world.Rank() == 0) {
        ModelPart& smp_1 = main_model_part.CreateSubModelPart("sub");
        smp_1.CreateSubModelPart("sub_sub");
        main_model_part.CreateSubModelPart("another_sub");
    }

    DistributedModelPartInitializer(main_model_part, Testing::GetDefaultDataCommunicator(), 0).Execute();

    KRATOS_EXPECT_TRUE(main_model_part.IsDistributed());

    KRATOS_EXPECT_EQ(main_model_part.NumberOfSubModelParts(), 2);
    KRATOS_EXPECT_TRUE(main_model_part.HasSubModelPart("sub"));
    KRATOS_EXPECT_TRUE(main_model_part.HasSubModelPart("another_sub"));
    KRATOS_EXPECT_TRUE(main_model_part.GetSubModelPart("sub").IsDistributed());
    KRATOS_EXPECT_TRUE(main_model_part.GetSubModelPart("another_sub").IsDistributed());

    ModelPart& sub_mp = main_model_part.GetSubModelPart("sub");
    KRATOS_EXPECT_EQ(sub_mp.NumberOfSubModelParts(), 1);
    KRATOS_EXPECT_TRUE(sub_mp.HasSubModelPart("sub_sub"));
    KRATOS_EXPECT_TRUE(sub_mp.GetSubModelPart("sub_sub").IsDistributed());
    KRATOS_EXPECT_FALSE(main_model_part.HasSubModelPart("sub_sub"));
    KRATOS_EXPECT_FALSE(main_model_part.GetSubModelPart("another_sub").HasSubModelPart("sub_sub"));

    KRATOS_EXPECT_EQ(main_model_part.GetCommunicator().GlobalNumberOfNodes(), 0);
}

}  // namespace Kratos::Testing
