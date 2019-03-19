//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "includes/model_part_io.h"
#include "includes/kratos_application.h"
#include "includes/kernel.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(
    ModelPartIOSubModelPartsDivision, KratosCoreFastSuite) {
    Kratos::shared_ptr<std::iostream> p_input(new std::stringstream(
        R"input(
				Begin ModelPartData
                                 DENSITY 2700.000000
				 BODY_FORCE [3] (0.000000,0.000000,0.000000)
				End ModelPartData


				 Begin Properties  1
				 DENSITY 2700.000000
				 YOUNG_MODULUS 7000000.000000
				 POISSON_RATIO 0.300000
				 BODY_FORCE [3] (0.000000,0.000000,0.000000)
				 THICKNESS 1.000000
				 End Properties

				Begin Nodes
				       1        0.0        0.0         0.0                               //node number, coord x, cord y, coord z
				       2        1.0        0.0         0.0                               //node number, coord x, cord y, coord z
				       3        1.0        1.0         0.0                               //node number, coord x, cord y, coord z
				       4        0.0        1.0         0.0                               //node number, coord x, cord y, coord z
				End Nodes

				Begin Elements Element2D3N
				  1 1 1 2 4  //the first column is the property
				  2 1 3 4 2
				End Elements

				Begin NodalData DISPLACEMENT_X          //be careful, variables are case sensitive!
				1 1 100.0                // pos1 is the node, pos2 (a 1) means that the DOF is fixed, then (position 3) we write the fixed displacement (in this case, temperature)
				End NodalData

				Begin NodalData DISPLACEMENT_Y
				1 1 100.0
				End NodalData

				Begin NodalData FORCE_Y
				3    0    5.0             //fixing it or not does not change anything since it is not a degree of freedom, it's just info that will be used by the condition
				End NodalData

				Begin Conditions Condition2D2N
				1 1 1 2
                                2 1 3 4
				End Conditions

				Begin SubModelPart BasePart // Note that this would be a sub sub modelpart
                                   Begin SubModelPartData
                                     DENSITY 1700.000000
                                     BODY_FORCE [3] (1.000000,1.000000,1.000000)
                                   End SubModelPartData
				   Begin SubModelPartNodes
				     1
				     2
				   End SubModelPartNodes
                                   Begin SubModelPart inner_part
                                        Begin SubModelPartNodes
                                            1
                                        End SubModelPartNodes
                                        Begin SubModelPartConditions
                                            1
                                        End SubModelPartConditions
                                    End SubModelPart
                                End SubModelPart
			)input"));

    Kernel kernel;
    KratosApplication application(std::string("Kratos"));
    application.Register();
    kernel.Initialize();

    Model current_model;

    ModelPartIO model_part_io(p_input);

    Kratos::shared_ptr<std::stringstream> p_output_0(new std::stringstream);
    Kratos::shared_ptr<std::stringstream> p_output_1(new std::stringstream);
    Kratos::shared_ptr<std::iostream> streams[2] = {p_output_0, p_output_1};

    std::size_t number_of_partitions = 2;
    std::size_t number_of_colors = 1;
    IO::GraphType domains_colored_graph(number_of_partitions, number_of_colors);
    domains_colored_graph(0, 0) = 1;
    domains_colored_graph(1, 0) = 0;
    IO::PartitionIndicesType nodes_partitions = {0, 0, 1, 1};
    IO::PartitionIndicesType elements_partitions = {0, 1};
    IO::PartitionIndicesType conditions_partitions = {0, 1};
    IO::PartitionIndicesContainerType nodes_all_partitions = {
        {0}, {0, 1}, {1}, {0, 1}};
    IO::PartitionIndicesContainerType elements_all_partitions = {{0}, {1}};
    IO::PartitionIndicesContainerType conditions_all_partitions = {{0}, {1}};

    model_part_io.DivideInputToPartitions(streams, number_of_partitions,
        domains_colored_graph, nodes_partitions, elements_partitions,
        conditions_partitions, nodes_all_partitions, elements_all_partitions,
        conditions_all_partitions);

    ModelPart& model_part_0 = current_model.CreateModelPart("Partition 0");
    model_part_0.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part_0.AddNodalSolutionStepVariable(FORCE);
    model_part_0.AddNodalSolutionStepVariable(PARTITION_INDEX);
    model_part_0.GetCommunicator().SetNumberOfColors(number_of_colors);

    ModelPartIO * model_part_io_0 = new ModelPartIO(p_output_0);
    model_part_io_0->ReadModelPart(model_part_0);

    KRATOS_CHECK_EQUAL(model_part_0.NumberOfNodes(), 3);
    KRATOS_CHECK_EQUAL(model_part_0.NumberOfElements(), 1);
    KRATOS_CHECK_EQUAL(model_part_0.NumberOfConditions(), 1);
    KRATOS_CHECK_EQUAL(model_part_0.NumberOfSubModelParts(), 1);
    KRATOS_CHECK(model_part_0.HasSubModelPart("BasePart"));
    KRATOS_CHECK_EQUAL(
        model_part_0.GetSubModelPart("BasePart").NumberOfNodes(), 2);
    KRATOS_CHECK_EQUAL(
        model_part_0.GetSubModelPart("BasePart").NumberOfElements(), 0);
    KRATOS_CHECK_EQUAL(
        model_part_0.GetSubModelPart("BasePart").NumberOfConditions(), 1);

    ModelPart& model_part_1 = current_model.CreateModelPart("Partition 1");
    model_part_1.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part_1.AddNodalSolutionStepVariable(FORCE);
    model_part_1.AddNodalSolutionStepVariable(PARTITION_INDEX);
    model_part_1.GetCommunicator().SetNumberOfColors(number_of_colors + 1);

    ModelPartIO * model_part_io_1 = new ModelPartIO(p_output_1);
    model_part_io_1->ReadModelPart(model_part_1);

    KRATOS_CHECK_EQUAL(model_part_1.NumberOfNodes(), 3);
    KRATOS_CHECK_EQUAL(model_part_1.NumberOfElements(), 1);
    KRATOS_CHECK_EQUAL(model_part_1.NumberOfConditions(), 1);
    KRATOS_CHECK_EQUAL(model_part_1.NumberOfSubModelParts(), 1);
    KRATOS_CHECK(model_part_1.HasSubModelPart("BasePart"));
    KRATOS_CHECK_EQUAL(
        model_part_1.GetSubModelPart("BasePart").NumberOfNodes(), 1);
    KRATOS_CHECK_EQUAL(
        model_part_1.GetSubModelPart("BasePart").NumberOfElements(), 0);
    KRATOS_CHECK_EQUAL(
        model_part_1.GetSubModelPart("BasePart").NumberOfConditions(), 0);

    // Free the modelparts to prevent files still being opened
    delete model_part_io_0;
    delete model_part_io_1;

    //KRATOS_CHECK_STRING_CONTAIN_SUB_STRING(p_output_1->str(), R"/(Begin SubModelPart BaseNodes
    //Begin SubModelPartNodes
    //2
    //End SubModelPartNodes
    //End Mesh)/");
}

KRATOS_TEST_CASE_IN_SUITE(ModelPartIOWriteModelPart, KratosCoreFastSuite) {

    Model current_model;
    
    // Create a model part to write
    ModelPart& main_model_part = current_model.CreateModelPart("MainModelPart");
    main_model_part.SetBufferSize(1);
    Properties::Pointer p_properties_1(new Properties(1));
    p_properties_1->SetValue(DENSITY, 1000.0);
    p_properties_1->SetValue(VISCOSITY, 1E-05);
    main_model_part.AddProperties(p_properties_1);

    main_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    main_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    main_model_part.CreateNewNode(3, 1.0, 1.0, 0.0);
    main_model_part.CreateNewNode(4, 0.0, 1.0, 0.0);

    std::vector<ModelPart::IndexType> elem_nodes_1 = {1,2,4};
    std::vector<ModelPart::IndexType> elem_nodes_2 = {3,4,2};
    main_model_part.CreateNewElement("Element2D3N", 1, elem_nodes_1, p_properties_1);
    main_model_part.CreateNewElement("Element2D3N", 2, elem_nodes_2, p_properties_1);

    //elemental data
    Matrix stress = ScalarMatrix(2,2, 1.00);
    main_model_part.GetMesh().GetElement(1).SetValue(CAUCHY_STRESS_TENSOR,stress);
    bool is_restarted = true;
    main_model_part.GetMesh().GetElement(1).SetValue(IS_RESTARTED,is_restarted);
    int domain_size =2;
    main_model_part.GetMesh().GetElement(1).SetValue(DOMAIN_SIZE,domain_size);
    double temperature = 34.7;
    main_model_part.GetMesh().GetElement(1).SetValue(TEMPERATURE, temperature);
    double displacement_x = 1.2;
    main_model_part.GetMesh().GetElement(1).SetValue(DISPLACEMENT_X, displacement_x);


    std::vector<ModelPart::IndexType> cond_nodes_1 = {1,2};
    std::vector<ModelPart::IndexType> cond_nodes_2 = {3,4};
    std::vector<ModelPart::IndexType> cond_nodes_3 = {4};
    main_model_part.CreateNewCondition("Condition2D2N", 1, cond_nodes_1, p_properties_1);
    main_model_part.CreateNewCondition("Condition2D2N", 2, cond_nodes_2, p_properties_1);
    main_model_part.CreateNewCondition("PointCondition2D1N", 3, cond_nodes_3, p_properties_1);

    //conditional data
    main_model_part.GetMesh().GetCondition(1).SetValue(CAUCHY_STRESS_TENSOR,stress);
    main_model_part.GetMesh().GetCondition(1).SetValue(IS_RESTARTED,is_restarted);
    main_model_part.GetMesh().GetCondition(1).SetValue(DOMAIN_SIZE,domain_size);
    main_model_part.GetMesh().GetCondition(1).SetValue(TEMPERATURE, temperature);
    main_model_part.GetMesh().GetCondition(1).SetValue(DISPLACEMENT_X, displacement_x);

    ModelPart* p_sub_model_part = &main_model_part.CreateSubModelPart("SubModelPart");
    std::vector<ModelPart::IndexType> sub_model_part_nodes = {1,2,4};
    std::vector<ModelPart::IndexType> sub_model_part_elems = {1};
    std::vector<ModelPart::IndexType> sub_model_part_conds = {1,3};
    p_sub_model_part->AddNodes(sub_model_part_nodes);
    p_sub_model_part->AddElements(sub_model_part_elems);
    p_sub_model_part->AddConditions(sub_model_part_conds);

    // Create the output .mdpa file
    std::string output_file_name = "main_model_part_output";
    std::fstream output_file;
    output_file.open(output_file_name + ".mdpa", std::fstream::out);
    output_file.close();

    // Fill the output .mdpa file
    ModelPartIO * model_part_io_write = new ModelPartIO(output_file_name, IO::WRITE);
    model_part_io_write->WriteModelPart(main_model_part);

    // Read and check the written .mdpa file
    ModelPartIO * model_part_io_output = new ModelPartIO(output_file_name);
    ModelPart& main_model_part_output = current_model.CreateModelPart("MainModelPartOutput");
    model_part_io_output->ReadModelPart(main_model_part_output);

    // Assert results
    KRATOS_CHECK_EQUAL(main_model_part_output.NumberOfProperties(), 1);
    KRATOS_CHECK_EQUAL(main_model_part_output.NumberOfSubModelParts() ,1);
    KRATOS_CHECK_EQUAL(main_model_part_output.NumberOfNodes(), 4);
    KRATOS_CHECK_EQUAL(main_model_part_output.NumberOfElements(), 2);
    KRATOS_CHECK_EQUAL(main_model_part_output.NumberOfConditions(), 3);
    KRATOS_CHECK_EQUAL(main_model_part_output.GetSubModelPart("SubModelPart").NumberOfNodes(), 3);
    KRATOS_CHECK_EQUAL(main_model_part_output.GetSubModelPart("SubModelPart").NumberOfElements(), 1);
    KRATOS_CHECK_EQUAL(main_model_part_output.GetSubModelPart("SubModelPart").NumberOfConditions(), 2);
    KRATOS_CHECK_EQUAL(main_model_part_output.GetMesh().GetElement(1).GetValue(CAUCHY_STRESS_TENSOR)(1,1),1);
    KRATOS_CHECK_EQUAL(main_model_part_output.GetMesh().GetElement(1).GetValue(IS_RESTARTED),is_restarted);
    KRATOS_CHECK_EQUAL(main_model_part_output.GetMesh().GetElement(1).GetValue(DOMAIN_SIZE),domain_size);
    KRATOS_CHECK_EQUAL(main_model_part_output.GetMesh().GetElement(1).GetValue(TEMPERATURE), temperature);
    KRATOS_CHECK_EQUAL(main_model_part_output.GetMesh().GetElement(1).GetValue(DISPLACEMENT_X), displacement_x);
    KRATOS_CHECK_EQUAL(main_model_part_output.GetMesh().GetCondition(1).GetValue(CAUCHY_STRESS_TENSOR)(1,1),1);
    KRATOS_CHECK_EQUAL(main_model_part_output.GetMesh().GetCondition(1).GetValue(IS_RESTARTED),is_restarted);
    KRATOS_CHECK_EQUAL(main_model_part_output.GetMesh().GetCondition(1).GetValue(DOMAIN_SIZE),domain_size);
    KRATOS_CHECK_EQUAL(main_model_part_output.GetMesh().GetCondition(1).GetValue(TEMPERATURE), temperature);
    KRATOS_CHECK_EQUAL(main_model_part_output.GetMesh().GetCondition(1).GetValue(DISPLACEMENT_X), displacement_x);

    // Free the modelparts to prevent files still being opened
    delete model_part_io_write;
    delete model_part_io_output;

    // Remove the generated files
    std::string aux_string_mdpa = output_file_name + ".mdpa";
    std::string aux_string_time = output_file_name + ".time";
    const char *mdpa_to_remove = aux_string_mdpa.c_str();
    const char *time_to_remove = aux_string_time.c_str();
    std::string error_msg = "Error deleting test output file: " + output_file_name;
    if (remove(mdpa_to_remove) != 0) {
        KRATOS_ERROR << error_msg + ".mdpa";
    }
    if (remove(time_to_remove) != 0) {
        KRATOS_ERROR << error_msg + ".time";
    }
}


KRATOS_TEST_CASE_IN_SUITE(ModelPartIOVariableNotInSolutionStepData, KratosCoreFastSuite) {
    Kratos::shared_ptr<std::iostream> p_input(new std::stringstream(
        R"input(
			    Begin Properties  0
                End Properties

				Begin Nodes
				       1        0.0        0.0         0.0
				       2        1.0        0.0         0.0
				       3        1.0        1.0         0.0
				       4        0.0        1.0         0.0
				End Nodes

				Begin Elements Element2D3N
				  1 0 1 2 4
				  2 0 3 4 2
				End Elements

				Begin NodalData DISPLACEMENT_X
				1 1 100.0
				End NodalData

				Begin NodalData DISPLACEMENT_Y
				1 1 100.0
				End NodalData

				Begin NodalData PRESSURE
				2 0 50.0
				End NodalData

				Begin NodalData TEMPERATURE
				4 0 33.0
				End NodalData

				Begin NodalData FORCE_Y
				3    0    5.0
				End NodalData
			)input"));

    Kernel kernel;
    KratosApplication application(std::string("Kratos"));
    application.Register();
    kernel.Initialize();

    Model current_model;

    // 1. Reading without IGNORE flag -> Error
    ModelPartIO default_model_part_io(p_input);

    ModelPart& model_part_0 = current_model.CreateModelPart("ErrorForce");
    model_part_0.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part_0.AddNodalSolutionStepVariable(PRESSURE);
    model_part_0.AddNodalSolutionStepVariable(TEMPERATURE);

    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        default_model_part_io.ReadModelPart(model_part_0),
        "The nodal solution step container does not have this variable: FORCE_Y");

    ModelPart& model_part_1 = current_model.CreateModelPart("ErrorTemperature");
    model_part_1.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part_1.AddNodalSolutionStepVariable(FORCE);
    model_part_1.AddNodalSolutionStepVariable(PRESSURE);

    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        default_model_part_io.ReadModelPart(model_part_1),
        "The nodal solution step container does not have this variable: TEMPERATURE");

    // 2. Reading with IGNORE flag -> Not set
    ModelPartIO ignore_model_part_io(p_input,ModelPartIO::READ|ModelPartIO::IGNORE_VARIABLES_ERROR);

    ModelPart& model_part_2 = current_model.CreateModelPart("IgnoreForce");
    model_part_0.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part_0.AddNodalSolutionStepVariable(PRESSURE);
    model_part_0.AddNodalSolutionStepVariable(TEMPERATURE);

    ignore_model_part_io.ReadModelPart(model_part_2);
    KRATOS_CHECK_EQUAL(model_part_2.NumberOfNodes(), 4);
    KRATOS_CHECK_EQUAL(model_part_2.NodesBegin()->Has(FORCE), false);

    ModelPart& model_part_3 = current_model.CreateModelPart("IgnoreTemperature");
    model_part_1.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part_1.AddNodalSolutionStepVariable(FORCE);
    model_part_1.AddNodalSolutionStepVariable(PRESSURE);

    ignore_model_part_io.ReadModelPart(model_part_3);
    KRATOS_CHECK_EQUAL(model_part_3.NumberOfNodes(), 4);
    KRATOS_CHECK_EQUAL(model_part_3.NodesBegin()->Has(TEMPERATURE), false);
}

}  // namespace Testing.
}  // namespace Kratos.
