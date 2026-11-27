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
#include "includes/kratos_components.h"
#include "tests/test_utilities/test_constitutive_law.h"

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

				Begin NodalData BOUNDARY
				    3
				End NodalData

				Begin NodalData DISPLACEMENT_X          //be careful, variables are case sensitive!
				    1 1 100.0                // pos1 is the node, pos2 (a 1) means that the DOF is fixed, then (position 3) we write the fixed displacement (in this case, temperature)
				End NodalData

				Begin NodalData DISPLACEMENT_Y
				    1 1 100.0
				End NodalData

				Begin NodalData FORCE_Y
				    3    0    5.0             //fixing it or not does not change anything since it is not a degree of freedom, it's just info that will be used by the condition
				End NodalData

				Begin Conditions LineCondition2D2N
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

    Model current_model;

    ModelPartIO model_part_io(p_input);

    Kratos::shared_ptr<std::stringstream> p_output_0(new std::stringstream);
    Kratos::shared_ptr<std::stringstream> p_output_1(new std::stringstream);
    Kratos::shared_ptr<std::iostream> streams[2] = {p_output_0, p_output_1};

    constexpr std::size_t number_of_partitions = 2;
    constexpr std::size_t number_of_colors = 1;
    IO::GraphType domains_colored_graph(number_of_partitions, number_of_colors);
    domains_colored_graph(0, 0) = 1;
    domains_colored_graph(1, 0) = 0;

    IO::PartitioningInfo part_info;
    part_info.Graph = domains_colored_graph;
    part_info.NodesPartitions = {0, 0, 1, 1};
    part_info.ElementsPartitions = {0, 1};
    part_info.ConditionsPartitions = {0, 1};
    part_info.NodesAllPartitions = {{0}, {0, 1}, {1}, {0, 1}};
    part_info.ElementsAllPartitions = {{0}, {1}};
    part_info.ConditionsAllPartitions = {{0}, {1}};

    model_part_io.DivideInputToPartitions(streams, number_of_partitions, part_info);

    ModelPart& model_part_0 = current_model.CreateModelPart("Partition 0");
    model_part_0.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part_0.AddNodalSolutionStepVariable(FORCE);
    model_part_0.AddNodalSolutionStepVariable(PARTITION_INDEX);
    model_part_0.GetCommunicator().SetNumberOfColors(number_of_colors);

    ModelPartIO * model_part_io_0 = new ModelPartIO(p_output_0);
    model_part_io_0->ReadModelPart(model_part_0);

    KRATOS_EXPECT_EQ(model_part_0.NumberOfNodes(), 3);
    KRATOS_EXPECT_EQ(model_part_0.NumberOfElements(), 1);
    KRATOS_EXPECT_EQ(model_part_0.NumberOfConditions(), 1);
    KRATOS_EXPECT_EQ(model_part_0.NumberOfSubModelParts(), 1);
    KRATOS_EXPECT_TRUE(model_part_0.HasSubModelPart("BasePart"));
    KRATOS_EXPECT_EQ(
        model_part_0.GetSubModelPart("BasePart").NumberOfNodes(), 2);
    KRATOS_EXPECT_EQ(
        model_part_0.GetSubModelPart("BasePart").NumberOfElements(), 0);
    KRATOS_EXPECT_EQ(
        model_part_0.GetSubModelPart("BasePart").NumberOfConditions(), 1);

    ModelPart& model_part_1 = current_model.CreateModelPart("Partition 1");
    model_part_1.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part_1.AddNodalSolutionStepVariable(FORCE);
    model_part_1.AddNodalSolutionStepVariable(PARTITION_INDEX);
    model_part_1.GetCommunicator().SetNumberOfColors(number_of_colors + 1);

    ModelPartIO * model_part_io_1 = new ModelPartIO(p_output_1);
    model_part_io_1->ReadModelPart(model_part_1);

    KRATOS_EXPECT_TRUE(model_part_1.GetNode(3).Is(BOUNDARY));
    KRATOS_EXPECT_EQ(model_part_1.NumberOfNodes(), 3);
    KRATOS_EXPECT_EQ(model_part_1.NumberOfElements(), 1);
    KRATOS_EXPECT_EQ(model_part_1.NumberOfConditions(), 1);
    KRATOS_EXPECT_EQ(model_part_1.NumberOfSubModelParts(), 1);
    KRATOS_EXPECT_TRUE(model_part_1.HasSubModelPart("BasePart"));
    KRATOS_EXPECT_EQ(
        model_part_1.GetSubModelPart("BasePart").NumberOfNodes(), 1);
    KRATOS_EXPECT_EQ(
        model_part_1.GetSubModelPart("BasePart").NumberOfElements(), 0);
    KRATOS_EXPECT_EQ(
        model_part_1.GetSubModelPart("BasePart").NumberOfConditions(), 0);

    // Free the modelparts to prevent files still being opened
    delete model_part_io_0;
    delete model_part_io_1;

    //KRATOS_EXPECT_HAS_SUBSTRING(p_output_1->str(), R"/(Begin SubModelPart BaseNodes
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
    main_model_part.CreateNewCondition("LineCondition2D2N", 1, cond_nodes_1, p_properties_1);
    main_model_part.CreateNewCondition("LineCondition2D2N", 2, cond_nodes_2, p_properties_1);
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
    KRATOS_EXPECT_EQ(main_model_part_output.NumberOfProperties(), 1);
    KRATOS_EXPECT_EQ(main_model_part_output.NumberOfSubModelParts() ,1);
    KRATOS_EXPECT_EQ(main_model_part_output.NumberOfNodes(), 4);
    KRATOS_EXPECT_EQ(main_model_part_output.NumberOfElements(), 2);
    KRATOS_EXPECT_EQ(main_model_part_output.NumberOfConditions(), 3);
    KRATOS_EXPECT_EQ(main_model_part_output.GetSubModelPart("SubModelPart").NumberOfNodes(), 3);
    KRATOS_EXPECT_EQ(main_model_part_output.GetSubModelPart("SubModelPart").NumberOfElements(), 1);
    KRATOS_EXPECT_EQ(main_model_part_output.GetSubModelPart("SubModelPart").NumberOfConditions(), 2);
    KRATOS_EXPECT_EQ(main_model_part_output.GetMesh().GetElement(1).GetValue(CAUCHY_STRESS_TENSOR)(1,1),1);
    KRATOS_EXPECT_EQ(main_model_part_output.GetMesh().GetElement(1).GetValue(IS_RESTARTED),is_restarted);
    KRATOS_EXPECT_EQ(main_model_part_output.GetMesh().GetElement(1).GetValue(DOMAIN_SIZE),domain_size);
    KRATOS_EXPECT_EQ(main_model_part_output.GetMesh().GetElement(1).GetValue(TEMPERATURE), temperature);
    KRATOS_EXPECT_EQ(main_model_part_output.GetMesh().GetElement(1).GetValue(DISPLACEMENT_X), displacement_x);
    KRATOS_EXPECT_EQ(main_model_part_output.GetMesh().GetCondition(1).GetValue(CAUCHY_STRESS_TENSOR)(1,1),1);
    KRATOS_EXPECT_EQ(main_model_part_output.GetMesh().GetCondition(1).GetValue(IS_RESTARTED),is_restarted);
    KRATOS_EXPECT_EQ(main_model_part_output.GetMesh().GetCondition(1).GetValue(DOMAIN_SIZE),domain_size);
    KRATOS_EXPECT_EQ(main_model_part_output.GetMesh().GetCondition(1).GetValue(TEMPERATURE), temperature);
    KRATOS_EXPECT_EQ(main_model_part_output.GetMesh().GetCondition(1).GetValue(DISPLACEMENT_X), displacement_x);

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

    Model current_model;

    // 1. Reading without IGNORE flag -> Error
    ModelPartIO default_model_part_io(p_input);

    ModelPart& model_part_0 = current_model.CreateModelPart("ErrorForce");
    model_part_0.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part_0.AddNodalSolutionStepVariable(PRESSURE);
    model_part_0.AddNodalSolutionStepVariable(TEMPERATURE);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        default_model_part_io.ReadModelPart(model_part_0),
        "The nodal solution step container does not have this variable: FORCE_Y");

    ModelPart& model_part_1 = current_model.CreateModelPart("ErrorTemperature");
    model_part_1.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part_1.AddNodalSolutionStepVariable(FORCE);
    model_part_1.AddNodalSolutionStepVariable(PRESSURE);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        default_model_part_io.ReadModelPart(model_part_1),
        "The nodal solution step container does not have this variable: TEMPERATURE");

    // 2. Reading with IGNORE flag -> Not set
    ModelPartIO ignore_model_part_io(p_input,ModelPartIO::READ|ModelPartIO::IGNORE_VARIABLES_ERROR);

    ModelPart& model_part_2 = current_model.CreateModelPart("IgnoreForce");
    model_part_0.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part_0.AddNodalSolutionStepVariable(PRESSURE);
    model_part_0.AddNodalSolutionStepVariable(TEMPERATURE);

    ignore_model_part_io.ReadModelPart(model_part_2);
    KRATOS_EXPECT_EQ(model_part_2.NumberOfNodes(), 4);
    KRATOS_EXPECT_EQ(model_part_2.NodesBegin()->Has(FORCE), false);

    ModelPart& model_part_3 = current_model.CreateModelPart("IgnoreTemperature");
    model_part_1.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part_1.AddNodalSolutionStepVariable(FORCE);
    model_part_1.AddNodalSolutionStepVariable(PRESSURE);

    ignore_model_part_io.ReadModelPart(model_part_3);
    KRATOS_EXPECT_EQ(model_part_3.NumberOfNodes(), 4);
    KRATOS_EXPECT_EQ(model_part_3.NodesBegin()->Has(TEMPERATURE), false);
}

KRATOS_TEST_CASE_IN_SUITE(ModelPartIOExtractValueFailsOnMalformedValue, KratosCoreFastSuite) {
    // A malformed value in a Nodes line ("bad" instead of a coordinate) used to be
    // silently ignored by ExtractValue (value >> rValue without any failure check),
    // leaving the coordinate untouched (indeterminate) and feeding it to CreateNewNode.
    Kratos::shared_ptr<std::iostream> p_input(new std::stringstream(
        R"input(
                Begin Nodes
                    1 bad 0.0 0.0
                    2 1.0 0.0 0.0
                End Nodes
        )input"));

    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");

    ModelPartIO model_part_io(p_input);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(model_part_io.ReadModelPart(r_model_part), "is not a valid value");
}

KRATOS_TEST_CASE_IN_SUITE(ModelPartIODivideInputDetectsMissingInterface, KratosCoreFastSuite) {
    // If the partition graph does not assign an interface between two partitions that
    // actually share nodes, interface_indices keeps its -1 initialization value. Using
    // that -1 as a SizeType index used to index out of bounds (huge index) instead of
    // raising an error about the inconsistent partition graph.
    Kratos::shared_ptr<std::iostream> p_input(new std::stringstream(
        R"input(
                Begin ModelPartData
                End ModelPartData

                Begin Nodes
                    1 0.0 0.0 0.0
                    2 1.0 0.0 0.0
                    3 1.0 1.0 0.0
                    4 0.0 1.0 0.0
                End Nodes
        )input"));

    ModelPartIO model_part_io(p_input);

    Kratos::shared_ptr<std::stringstream> p_output_0(new std::stringstream);
    Kratos::shared_ptr<std::stringstream> p_output_1(new std::stringstream);
    Kratos::shared_ptr<std::iostream> streams[2] = {p_output_0, p_output_1};

    constexpr std::size_t number_of_partitions = 2;
    constexpr std::size_t number_of_colors = 1;
    IO::GraphType domains_colored_graph(number_of_partitions, number_of_colors);
    // No interface assigned anywhere: entries stay -1 although node 2 is shared
    // between partitions 0 and 1 (NodesAllPartitions below).
    domains_colored_graph(0, 0) = -1;
    domains_colored_graph(1, 0) = -1;

    IO::PartitioningInfo part_info;
    part_info.Graph = domains_colored_graph;
    part_info.NodesPartitions = {0, 0, 1, 1};
    part_info.ElementsPartitions = {};
    part_info.ConditionsPartitions = {};
    part_info.NodesAllPartitions = {{0}, {0, 1}, {1}, {0, 1}};
    part_info.ElementsAllPartitions = {};
    part_info.ConditionsAllPartitions = {};

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        model_part_io.DivideInputToPartitions(streams, number_of_partitions, part_info),
        "is not connected");
}

KRATOS_TEST_CASE_IN_SUITE(ModelPartIOWritePropertiesKeepsFieldsAfterConstitutiveLaw, KratosCoreFastSuite) {
    // When a CONSTITUTIVE_LAW is present, WriteProperties removes the (unprintable)
    // constitutive law line coming from PrintData and re-appends it by hand with the
    // registered name. The erase used erase(begin, end_offset) as if it were
    // erase(begin, end): it removes end_offset characters starting at begin, i.e. the
    // CL line plus (begin - 1) characters of everything stored after it. The damage
    // grows with the position of the CL line inside the block (PrintData keeps
    // insertion order), so fields stored after a late CL are silently lost.
    static TestConstitutiveLaw s_dummy_constitutive_law;
    if (!KratosComponents<ConstitutiveLaw>::Has("TestConstitutiveLaw"))
        KratosComponents<ConstitutiveLaw>::Add("TestConstitutiveLaw", s_dummy_constitutive_law);

    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");

    Properties::Pointer p_properties(new Properties(1));
    // Fill the block BEFORE the constitutive law so that the "CONSTITUTIVE_LAW"
    // substring starts deep in the string: the buggy erase(begin, end_offset)
    // removes count = end_offset characters from begin, i.e. the CL line PLUS
    // begin-1 characters of everything that follows it. With CL first the extra
    // is only 3 characters and the damage stays hidden.
    p_properties->SetValue(POISSON_RATIO, 0.3);
    p_properties->SetValue(THICKNESS, 1.0);
    p_properties->SetValue(VISCOSITY, 1.0e-5);
    p_properties->SetValue(CONDUCTIVITY, 50.0);
    p_properties->SetValue(SPECIFIC_HEAT, 900.0);
    p_properties->SetValue(CONVECTION_COEFFICIENT, 10.0);
    p_properties->SetValue(AMBIENT_TEMPERATURE, 293.0);
    array_1d<double, 3> body_force;
    body_force[0] = 0.0;
    body_force[1] = 0.0;
    body_force[2] = 0.0;
    p_properties->SetValue(BODY_FORCE, body_force);
    p_properties->SetValue(CONSTITUTIVE_LAW, ConstitutiveLaw::Pointer(new TestConstitutiveLaw()));
    // Canary fields stored AFTER the constitutive law (PrintData keeps insertion order)
    p_properties->SetValue(DENSITY, 2700.0);
    p_properties->SetValue(YOUNG_MODULUS, 7000000.0);
    r_model_part.AddProperties(p_properties);

    Kratos::shared_ptr<std::stringstream> p_output(new std::stringstream);
    ModelPartIO model_part_io(p_output, IO::WRITE | IO::SKIP_TIMER);
    model_part_io.WriteModelPart(r_model_part);

    const std::string output = p_output->str();
    KRATOS_EXPECT_TRUE(output.find("Begin Properties 1") != std::string::npos);
    KRATOS_EXPECT_TRUE(output.find("DENSITY") != std::string::npos);
    KRATOS_EXPECT_TRUE(output.find("YOUNG_MODULUS") != std::string::npos);
    KRATOS_EXPECT_TRUE(output.find("CONSTITUTIVE_LAW TestConstitutiveLaw") != std::string::npos);

    // The constitutive law line must appear exactly once (removed + re-added by hand)
    std::size_t count = 0;
    for (std::size_t pos = output.find("CONSTITUTIVE_LAW"); pos != std::string::npos; pos = output.find("CONSTITUTIVE_LAW", pos + 1))
        count++;
    KRATOS_EXPECT_EQ(count, 1u);
}

KRATOS_TEST_CASE_IN_SUITE(ModelPartIOEmptyConstraintDeclarationThrows, KratosCoreFastSuite) {
    // The line after "Begin Constraints <Name>" carries the slave (first word) and
    // master variables. With an empty line, words.size()-1 underflows to ~SizeType(0)
    // (and words[0] reads out of bounds), used to read an absurd number of master
    // DOFs. Note that a single word is legal (1x0 constraints): only the empty case
    // must be rejected.
    Kratos::shared_ptr<std::iostream> p_input(new std::stringstream(
        R"input(
                Begin ModelPartData
                End ModelPartData

                Begin Constraints LinearMasterSlaveConstraint

                End Constraints
        )input"));

    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");

    ModelPartIO model_part_io(p_input);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        model_part_io.ReadModelPart(r_model_part),
        "Empty constraint declaration line");
}

KRATOS_TEST_CASE_IN_SUITE(ModelPartIOWriteTableIncludesFirstRow, KratosCoreFastSuite) {
    // WriteTableBlock iterated the table points from j = 1, silently dropping the
    // first point of every table, while ReadTableBlock reads every point.
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");

    ModelPart::TableType::Pointer p_table = Kratos::make_shared<ModelPart::TableType>();
    // NOTE: stream default precision prints only 6 significant digits, so the
    // table entries use short values that survive operator<< verbatim.
    p_table->insert(0.25, 8.75);
    p_table->insert(0.75, 9.5);
    r_model_part.AddTable(1, p_table);

    Kratos::shared_ptr<std::stringstream> p_output(new std::stringstream);
    ModelPartIO model_part_io(p_output, IO::WRITE | IO::SKIP_TIMER);
    model_part_io.WriteModelPart(r_model_part);

    const std::string output = p_output->str();
    KRATOS_EXPECT_TRUE(output.find("Begin Table") != std::string::npos);
    KRATOS_EXPECT_TRUE(output.find("0.25") != std::string::npos); // first row (used to be dropped)
    KRATOS_EXPECT_TRUE(output.find("8.75") != std::string::npos);
    KRATOS_EXPECT_TRUE(output.find("0.75") != std::string::npos); // second row
    KRATOS_EXPECT_TRUE(output.find("9.5") != std::string::npos);
}

}  // namespace Testing.
}  // namespace Kratos.
