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
#include "testing/testing.h"
#include "includes/model_part_io.h"
#include "includes/kratos_application.h"
#include "includes/kernel.h"


namespace Kratos {
	namespace Testing {

		KRATOS_TEST_CASE_IN_SUITE(ModelPartIOSubModelPartsDivision, KratosCoreFastSuite)
		{
			boost::shared_ptr<std::iostream> p_input(new std::stringstream( 
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
                        
                                Begin NodalData VELOCITY
                                1 0 [3] (0.1, 0.2, 0.3)
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
			)input"));

			Kernel kernel;
			KratosApplication application;
			application.Register();
			kernel.Initialize();

			ModelPartIO model_part_io(p_input);

			boost::shared_ptr<std::stringstream> p_output_0(new std::stringstream);
			boost::shared_ptr<std::stringstream> p_output_1(new std::stringstream);
			boost::shared_ptr<std::iostream> streams[2] = { p_output_0, p_output_1 };

			std::size_t number_of_partitions = 2;
			std::size_t number_of_colors = 1;
			IO::GraphType domains_colored_graph(number_of_partitions, number_of_colors);
			domains_colored_graph(0, 0) = 1;
			domains_colored_graph(1, 0) = 0;
			IO::PartitionIndicesType nodes_partitions = { 0,0,1,1 };
			IO::PartitionIndicesType elements_partitions = { 0,1 };
			IO::PartitionIndicesType conditions_partitions = { 0,1 };
			IO::PartitionIndicesContainerType nodes_all_partitions = { { 0 },{ 0,1 },{ 1 },{ 0,1 } };
			IO::PartitionIndicesContainerType elements_all_partitions = { { 0 },{ 1 } };
			IO::PartitionIndicesContainerType conditions_all_partitions = { { 0 },{ 1 } };

			model_part_io.DivideInputToPartitions(streams, number_of_partitions, domains_colored_graph,
				nodes_partitions,
				elements_partitions,
				conditions_partitions,
				nodes_all_partitions,
				elements_all_partitions,
				conditions_all_partitions);

			ModelPart model_part_0("Partition 0");
			model_part_0.AddNodalSolutionStepVariable(DISPLACEMENT);
			model_part_0.AddNodalSolutionStepVariable(FORCE);
			model_part_0.AddNodalSolutionStepVariable(PARTITION_INDEX);
			model_part_0.GetCommunicator().SetNumberOfColors(number_of_colors);

			ModelPartIO model_part_io_0(p_output_0);
			model_part_io_0.ReadModelPart(model_part_0);


			KRATOS_CHECK_EQUAL(model_part_0.NumberOfNodes(), 3);
			KRATOS_CHECK_EQUAL(model_part_0.NumberOfElements(), 1);
			KRATOS_CHECK_EQUAL(model_part_0.NumberOfConditions(), 1);
			KRATOS_CHECK_EQUAL(model_part_0.NumberOfSubModelParts(), 1);
			KRATOS_CHECK(model_part_0.HasSubModelPart("BasePart"));
			KRATOS_CHECK_EQUAL(model_part_0.GetSubModelPart("BasePart").NumberOfNodes(), 2);
			KRATOS_CHECK_EQUAL(model_part_0.GetSubModelPart("BasePart").NumberOfElements(), 0);
			KRATOS_CHECK_EQUAL(model_part_0.GetSubModelPart("BasePart").NumberOfConditions(), 1);

			ModelPart model_part_1("Partition 1");
			model_part_1.AddNodalSolutionStepVariable(DISPLACEMENT);
			model_part_1.AddNodalSolutionStepVariable(FORCE);
			model_part_1.AddNodalSolutionStepVariable(PARTITION_INDEX);
			model_part_1.GetCommunicator().SetNumberOfColors(number_of_colors + 1);

			ModelPartIO model_part_io_1(p_output_1);
			model_part_io_1.ReadModelPart(model_part_1);

			KRATOS_CHECK_EQUAL(model_part_1.NumberOfNodes(), 3);
			KRATOS_CHECK_EQUAL(model_part_1.NumberOfElements(), 1);
			KRATOS_CHECK_EQUAL(model_part_1.NumberOfConditions(), 1);
			KRATOS_CHECK_EQUAL(model_part_1.NumberOfSubModelParts(), 1);
			KRATOS_CHECK(model_part_1.HasSubModelPart("BasePart"));
			KRATOS_CHECK_EQUAL(model_part_1.GetSubModelPart("BasePart").NumberOfNodes(), 1);
			KRATOS_CHECK_EQUAL(model_part_1.GetSubModelPart("BasePart").NumberOfElements(), 0);
			KRATOS_CHECK_EQUAL(model_part_1.GetSubModelPart("BasePart").NumberOfConditions(), 0);


//KRATOS_CHECK_STRING_CONTAIN_SUB_STRING(p_output_1->str(), R"/(Begin SubModelPart BaseNodes
//Begin SubModelPartNodes
//2
//End SubModelPartNodes
//End Mesh)/");

		}
	}
}  // namespace Kratos.


