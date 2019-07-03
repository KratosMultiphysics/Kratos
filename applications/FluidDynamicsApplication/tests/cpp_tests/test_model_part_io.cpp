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
//                   Carlos A. Roig
//
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/model_part_io.h"

namespace Kratos {
namespace Testing {

// // Tests that two elements with the same geometry are properly distinguished 
// KRATOS_TEST_CASE_IN_SUITE(ModelPartIOWriteModelPartDiffElementSameGeometry, KratosCheckWriteElementsTestSuite) {
//     // Create a model part to write
//     ModelPart main_model_part("MainModelPart");
//     main_model_part.SetBufferSize(1);
//     Properties::Pointer p_properties_1(new Properties(1));
//     main_model_part.AddProperties(p_properties_1);

//     main_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
//     main_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
//     main_model_part.CreateNewNode(3, 1.0, 1.0, 0.0);

//     std::vector<ModelPart::IndexType> elem_nodes_1 = {1, 2, 3};
//     main_model_part.CreateNewElement("VMS2D3N", 1, elem_nodes_1, p_properties_1);
//     main_model_part.CreateNewElement("QSVMS2D3N", 2, elem_nodes_1, p_properties_1);

//     // Create the output .mdpa file
//     std::string output_file_name = "main_model_part_output";
//     std::fstream output_file;
//     output_file.open(output_file_name + ".mdpa", std::fstream::out);
//     output_file.close();

//     // Fill the output .mdpa file
//     ModelPartIO model_part_io_write(output_file_name, IO::WRITE);
//     model_part_io_write.WriteModelPart(main_model_part);

//     // Read and check the written .mdpa file
//     ModelPartIO model_part_io_output(output_file_name);
//     ModelPart main_model_part_output("MainModelPartOutput");
//     model_part_io_output.ReadModelPart(main_model_part_output);

//     // Remove the generated files
//     std::string aux_string_mdpa = output_file_name + ".mdpa"; 
//     std::string aux_string_time = output_file_name + ".time";
//     const char *mdpa_to_remove = aux_string_mdpa.c_str();
//     const char *time_to_remove = aux_string_time.c_str();
//     std::string error_msg = "Error deleting test output file: " + output_file_name;
//     if (remove(mdpa_to_remove) != 0) {
//         KRATOS_ERROR << error_msg + ".mdpa";
//     }
//     if (remove(time_to_remove) != 0) {
//         KRATOS_ERROR << error_msg + ".time";
//     }
// }

// // Tests that two elements with the same geometry are properly distinguished when added out of order and with repetitions
// KRATOS_TEST_CASE_IN_SUITE(ModelPartIOWriteModelPartDiffElementSameGeometryUnordered, KratosCheckWriteElementsTestSuite) {
//     // Create a model part to write
//     ModelPart main_model_part("MainModelPart");
//     main_model_part.SetBufferSize(1);
//     Properties::Pointer p_properties_1(new Properties(1));
//     main_model_part.AddProperties(p_properties_1);

//     main_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
//     main_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
//     main_model_part.CreateNewNode(3, 1.0, 1.0, 0.0);

//     std::vector<ModelPart::IndexType> elem_nodes_1 = {1, 2, 3};
//     main_model_part.CreateNewElement("VMS2D3N", 1, elem_nodes_1, p_properties_1);
//     main_model_part.CreateNewElement("QSVMS2D3N", 2, elem_nodes_1, p_properties_1);
//     main_model_part.CreateNewElement("VMS2D3N", 3, elem_nodes_1, p_properties_1);
//     main_model_part.CreateNewElement("QSVMS2D3N", 4, elem_nodes_1, p_properties_1);

//     // Create the output .mdpa file
//     std::string output_file_name = "main_model_part_output";
//     std::fstream output_file;
//     output_file.open(output_file_name + ".mdpa", std::fstream::out);
//     output_file.close();

//     // Fill the output .mdpa file
//     ModelPartIO model_part_io_write(output_file_name, IO::WRITE);
//     model_part_io_write.WriteModelPart(main_model_part);

//     // Read and check the written .mdpa file
//     ModelPartIO model_part_io_output(output_file_name);
//     ModelPart main_model_part_output("MainModelPartOutput");
//     model_part_io_output.ReadModelPart(main_model_part_output);

//     // Remove the generated files
//     std::string aux_string_mdpa = output_file_name + ".mdpa"; 
//     std::string aux_string_time = output_file_name + ".time";
//     const char *mdpa_to_remove = aux_string_mdpa.c_str();
//     const char *time_to_remove = aux_string_time.c_str();
//     std::string error_msg = "Error deleting test output file: " + output_file_name;
//     if (remove(mdpa_to_remove) != 0) {
//         KRATOS_ERROR << error_msg + ".mdpa";
//     }
//     if (remove(time_to_remove) != 0) {
//         KRATOS_ERROR << error_msg + ".time";
//     }
// }

// // Tests that the same element with different geometries are properly distinguished
// KRATOS_TEST_CASE_IN_SUITE(ModelPartIOWriteModelPartSameElementDiffGeometry, KratosCheckWriteElementsTestSuite) {
//     // Create a model part to write
//     ModelPart main_model_part("MainModelPart");
//     main_model_part.SetBufferSize(1);
//     Properties::Pointer p_properties_1(new Properties(1));
//     main_model_part.AddProperties(p_properties_1);

//     main_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
//     main_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
//     main_model_part.CreateNewNode(3, 1.0, 1.0, 0.0);
//     main_model_part.CreateNewNode(4, 1.0, 1.0, 1.0);

//     std::vector<ModelPart::IndexType> elem_nodes_1 = {1, 2, 3};
//     std::vector<ModelPart::IndexType> elem_nodes_2 = {1, 2, 3, 4};
//     main_model_part.CreateNewElement("SpalartAllmaras2D", 1, elem_nodes_1, p_properties_1);
//     main_model_part.CreateNewElement("SpalartAllmaras3D", 2, elem_nodes_2, p_properties_1);

//     // Create the output .mdpa file
//     std::string output_file_name = "main_model_part_output";
//     std::fstream output_file;
//     output_file.open(output_file_name + ".mdpa", std::fstream::out);
//     output_file.close();

//     // Fill the output .mdpa file
//     ModelPartIO model_part_io_write(output_file_name, IO::WRITE);
//     model_part_io_write.WriteModelPart(main_model_part);

//     // Read and check the written .mdpa file
//     ModelPartIO model_part_io_output(output_file_name);
//     ModelPart main_model_part_output("MainModelPartOutput");
//     model_part_io_output.ReadModelPart(main_model_part_output);

//     // Remove the generated files
//     std::string aux_string_mdpa = output_file_name + ".mdpa"; 
//     std::string aux_string_time = output_file_name + ".time";
//     const char *mdpa_to_remove = aux_string_mdpa.c_str();
//     const char *time_to_remove = aux_string_time.c_str();
//     std::string error_msg = "Error deleting test output file: " + output_file_name;
//     if (remove(mdpa_to_remove) != 0) {
//         KRATOS_ERROR << error_msg + ".mdpa";
//     }
//     if (remove(time_to_remove) != 0) {
//         KRATOS_ERROR << error_msg + ".time";
//     }
// }

// // Tests that the same element with different geometries are properly distinguished when added out of order and with repetitions
// KRATOS_TEST_CASE_IN_SUITE(ModelPartIOWriteModelPartSameElementDiffGeometryUordered, KratosCheckWriteElementsTestSuite) {
//     // Create a model part to write
//     ModelPart main_model_part("MainModelPart");
//     main_model_part.SetBufferSize(1);
//     Properties::Pointer p_properties_1(new Properties(1));
//     main_model_part.AddProperties(p_properties_1);

//     main_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
//     main_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
//     main_model_part.CreateNewNode(3, 1.0, 1.0, 0.0);
//     main_model_part.CreateNewNode(4, 1.0, 1.0, 1.0);

//     std::vector<ModelPart::IndexType> elem_nodes_1 = {1, 2, 3};
//     std::vector<ModelPart::IndexType> elem_nodes_2 = {1, 2, 3, 4};
//     main_model_part.CreateNewElement("SpalartAllmaras2D", 1, elem_nodes_1, p_properties_1);
//     main_model_part.CreateNewElement("SpalartAllmaras3D", 2, elem_nodes_2, p_properties_1);
//     main_model_part.CreateNewElement("SpalartAllmaras2D", 3, elem_nodes_1, p_properties_1);
//     main_model_part.CreateNewElement("SpalartAllmaras3D", 4, elem_nodes_2, p_properties_1);

//     // Create the output .mdpa file
//     std::string output_file_name = "main_model_part_output";
//     std::fstream output_file;
//     output_file.open(output_file_name + ".mdpa", std::fstream::out);
//     output_file.close();

//     // Fill the output .mdpa file
//     ModelPartIO model_part_io_write(output_file_name, IO::WRITE);
//     model_part_io_write.WriteModelPart(main_model_part);

//     // Read and check the written .mdpa file
//     ModelPartIO model_part_io_output(output_file_name);
//     ModelPart main_model_part_output("MainModelPartOutput");
//     model_part_io_output.ReadModelPart(main_model_part_output);

//     // Remove the generated files
//     std::string aux_string_mdpa = output_file_name + ".mdpa"; 
//     std::string aux_string_time = output_file_name + ".time";
//     const char *mdpa_to_remove = aux_string_mdpa.c_str();
//     const char *time_to_remove = aux_string_time.c_str();
//     std::string error_msg = "Error deleting test output file: " + output_file_name;
//     if (remove(mdpa_to_remove) != 0) {
//         KRATOS_ERROR << error_msg + ".mdpa";
//     }
//     if (remove(time_to_remove) != 0) {
//         KRATOS_ERROR << error_msg + ".time";
//     }
// }

} // Namespace Testing
} // Namespace Kratos