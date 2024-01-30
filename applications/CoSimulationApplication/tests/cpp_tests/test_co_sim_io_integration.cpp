// KRATOS  / ___|___/ ___|(_)_ __ ___  _   _| | __ _| |_(_) ___  _ ___
//        | |   / _ \___ \| | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_  |
//        | |__| (_) |__) | | | | | | | |_| | | (_| | |_| | (_) | | | |
//         \____\___/____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|
//
//  License:		 BSD License
//                   license: CoSimulationApplication/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/parallel_environment.h"
#include "containers/model.h"
#include "testing/testing.h"
#include "custom_utilities/co_sim_io_conversion_utilities.h"
#include "co_simulation_testing_utilities.h"

// Utilities

namespace Kratos {
namespace Testing {

using DataLocation = Globals::DataLocation;

KRATOS_TEST_CASE_IN_SUITE(CoSimIOModelPartToKratosModelPart_NodesOnly, KratosCosimulationFastSuite)
{
    Model model;
    auto& kratos_model_part = model.CreateModelPart("kratos_mp");

    CoSimIO::ModelPart co_sim_io_model_part("co_sim_io_mp");

    constexpr std::size_t num_nodes = 5;

    for (std::size_t i=0; i<num_nodes; ++i) {
        co_sim_io_model_part.CreateNewNode(i+1, i*1.5, i+3.5, i-8.6);
    }

    KRATOS_EXPECT_EQ(co_sim_io_model_part.NumberOfNodes(), num_nodes);
    KRATOS_EXPECT_EQ(co_sim_io_model_part.NumberOfLocalNodes(), num_nodes);
    KRATOS_EXPECT_EQ(co_sim_io_model_part.NumberOfGhostNodes(), 0);
    KRATOS_EXPECT_EQ(co_sim_io_model_part.NumberOfElements(), 0);

    const auto& r_serial_data_comm = ParallelEnvironment::GetDataCommunicator("Serial");

    CoSimIOConversionUtilities::CoSimIOModelPartToKratosModelPart(co_sim_io_model_part, kratos_model_part, r_serial_data_comm);

    CheckModelPartsAreEqual(kratos_model_part, co_sim_io_model_part);
}

KRATOS_TEST_CASE_IN_SUITE(CoSimIOModelPartToKratosModelPart_NodesOnly_Unordered, KratosCosimulationFastSuite)
{
    Model model;
    auto& kratos_model_part = model.CreateModelPart("kratos_mp");

    CoSimIO::ModelPart co_sim_io_model_part("co_sim_io_mp");

    constexpr std::size_t num_nodes = 5;

    constexpr std::array<CoSimIO::IdType, num_nodes> node_ids {1,2,3,6,4};

    for (std::size_t i=0; i<num_nodes; ++i) {
        co_sim_io_model_part.CreateNewNode(node_ids[i], i*1.5, i+3.5, i-8.6);
    }

    KRATOS_EXPECT_EQ(co_sim_io_model_part.NumberOfNodes(), num_nodes);
    KRATOS_EXPECT_EQ(co_sim_io_model_part.NumberOfLocalNodes(), num_nodes);
    KRATOS_EXPECT_EQ(co_sim_io_model_part.NumberOfGhostNodes(), 0);
    KRATOS_EXPECT_EQ(co_sim_io_model_part.NumberOfElements(), 0);

    const auto& r_serial_data_comm = ParallelEnvironment::GetDataCommunicator("Serial");

    CoSimIOConversionUtilities::CoSimIOModelPartToKratosModelPart(co_sim_io_model_part, kratos_model_part, r_serial_data_comm);

    CheckModelPartsAreEqualButEntitiesAreOrderedDifferently(kratos_model_part, co_sim_io_model_part);
}

KRATOS_TEST_CASE_IN_SUITE(KratosModelPartToCoSimIOModelPart_NodesOnly, KratosCosimulationFastSuite)
{
    Model model;
    auto& kratos_model_part = model.CreateModelPart("kratos_mp");

    CoSimIO::ModelPart co_sim_io_model_part("co_sim_io_mp");

    constexpr std::size_t num_nodes = 5;

    for (std::size_t i=0; i<num_nodes; ++i) {
        kratos_model_part.CreateNewNode(i+1, i*1.5, i+3.5, i-8.6);
    }

    KRATOS_EXPECT_EQ(kratos_model_part.NumberOfNodes(), num_nodes);
    KRATOS_EXPECT_EQ(kratos_model_part.NumberOfProperties(), 0);

    CoSimIOConversionUtilities::KratosModelPartToCoSimIOModelPart(kratos_model_part, co_sim_io_model_part);

    CheckModelPartsAreEqual(kratos_model_part, co_sim_io_model_part);
}

KRATOS_TEST_CASE_IN_SUITE(CoSimIOModelPartToKratosModelPart, KratosCosimulationFastSuite)
{
    Model model;
    auto& kratos_model_part = model.CreateModelPart("kratos_mp");

    CoSimIO::ModelPart co_sim_io_model_part("co_sim_io_mp");

    const int node_ids[] = {2,61,159};
    const std::array<double, 3> node_coords = {1.0, -2.7, 9.44};
    co_sim_io_model_part.CreateNewNode(node_ids[0], node_coords[0], node_coords[1], node_coords[2]);
    co_sim_io_model_part.CreateNewNode(node_ids[1], node_coords[1], node_coords[2], node_coords[0]);
    co_sim_io_model_part.CreateNewNode(node_ids[2], node_coords[2], node_coords[0], node_coords[1]);

    KRATOS_EXPECT_EQ(co_sim_io_model_part.NumberOfNodes(), 3);

    const int elem_ids[] = {1,19,21};
    const CoSimIO::ElementType elem_types[] = {
        CoSimIO::ElementType::Point2D,
        CoSimIO::ElementType::Point2D,
        CoSimIO::ElementType::Line2D2
    };

    co_sim_io_model_part.CreateNewElement(elem_ids[0], elem_types[0], {node_ids[0]});
    co_sim_io_model_part.CreateNewElement(elem_ids[1], elem_types[1], {node_ids[1]});
    co_sim_io_model_part.CreateNewElement(elem_ids[2], elem_types[2], {node_ids[1], node_ids[2]});

    KRATOS_EXPECT_EQ(co_sim_io_model_part.NumberOfElements(), 3);

    const auto& r_serial_data_comm = ParallelEnvironment::GetDataCommunicator("Serial");
    CoSimIOConversionUtilities::CoSimIOModelPartToKratosModelPart(co_sim_io_model_part, kratos_model_part, r_serial_data_comm);

    CheckModelPartsAreEqual(kratos_model_part, co_sim_io_model_part);
}

KRATOS_TEST_CASE_IN_SUITE(CoSimIOModelPartToKratosModelPart_Unordered, KratosCosimulationFastSuite)
{
    Model model;
    auto& kratos_model_part = model.CreateModelPart("kratos_mp");

    CoSimIO::ModelPart co_sim_io_model_part("co_sim_io_mp");

    const int node_ids[] = {159, 2,61};
    const std::array<double, 3> node_coords = {1.0, -2.7, 9.44};
    co_sim_io_model_part.CreateNewNode(node_ids[0], node_coords[0], node_coords[1], node_coords[2]);
    co_sim_io_model_part.CreateNewNode(node_ids[1], node_coords[1], node_coords[2], node_coords[0]);
    co_sim_io_model_part.CreateNewNode(node_ids[2], node_coords[2], node_coords[0], node_coords[1]);

    KRATOS_EXPECT_EQ(co_sim_io_model_part.NumberOfNodes(), 3);

    const int elem_ids[] = {104,19,21};
    const CoSimIO::ElementType elem_types[] = {
        CoSimIO::ElementType::Point2D,
        CoSimIO::ElementType::Point2D,
        CoSimIO::ElementType::Line2D2
    };

    co_sim_io_model_part.CreateNewElement(elem_ids[0], elem_types[0], {node_ids[0]});
    co_sim_io_model_part.CreateNewElement(elem_ids[1], elem_types[1], {node_ids[1]});
    co_sim_io_model_part.CreateNewElement(elem_ids[2], elem_types[2], {node_ids[1], node_ids[2]});

    KRATOS_EXPECT_EQ(co_sim_io_model_part.NumberOfElements(), 3);

    const auto& r_serial_data_comm = ParallelEnvironment::GetDataCommunicator("Serial");
    CoSimIOConversionUtilities::CoSimIOModelPartToKratosModelPart(co_sim_io_model_part, kratos_model_part, r_serial_data_comm);

    CheckModelPartsAreEqualButEntitiesAreOrderedDifferently(kratos_model_part, co_sim_io_model_part);
}

KRATOS_TEST_CASE_IN_SUITE(KratosModelPartToCoSimIOModelPart, KratosCosimulationFastSuite)
{
    Model model;
    auto& kratos_model_part = model.CreateModelPart("kratos_mp");

    CoSimIO::ModelPart co_sim_io_model_part("co_sim_io_mp");

    constexpr std::size_t num_nodes = 5;

    for (std::size_t i=0; i<num_nodes; ++i) {
        kratos_model_part.CreateNewNode(i+1, i*1.5, i+3.5, i-8.6);
    }

    auto p_props = kratos_model_part.CreateNewProperties(0);

    std::vector<IndexType> conn {1,2};
    kratos_model_part.CreateNewElement("Element2D2N", 1, conn, p_props);
    conn = {2,3};
    kratos_model_part.CreateNewElement("Element2D2N", 2, conn, p_props);
    conn = {3,4,5};
    kratos_model_part.CreateNewElement("Element2D3N", 3, conn, p_props);

    KRATOS_EXPECT_EQ(kratos_model_part.NumberOfNodes(), num_nodes);
    KRATOS_EXPECT_EQ(kratos_model_part.NumberOfElements(), 3);
    KRATOS_EXPECT_EQ(kratos_model_part.NumberOfProperties(), 1);

    CoSimIOConversionUtilities::KratosModelPartToCoSimIOModelPart(kratos_model_part, co_sim_io_model_part);

    CheckModelPartsAreEqual(kratos_model_part, co_sim_io_model_part);
}

KRATOS_TEST_CASE_IN_SUITE(CoSimIODataExport_direct_scalar, KratosCosimulationFastSuite)
{
    Model model;
    auto& kratos_model_part = model.CreateModelPart("kratos_mp");
    kratos_model_part.AddNodalSolutionStepVariable(AUX_INDEX);

    constexpr std::size_t num_nodes = 5;
    const std::vector<double> exp_values {2,6,8,4,1};

    auto p_props = kratos_model_part.CreateNewProperties(0);

    for (std::size_t i=0; i<num_nodes; ++i) {
        auto p_node = kratos_model_part.CreateNewNode(i+1, i*1.5, i+3.5, i-8.6);
        auto p_elem = kratos_model_part.CreateNewElement("Element2D1N", i+1, std::vector<std::size_t>{i+1ul}, p_props);
        p_node->FastGetSolutionStepValue(AUX_INDEX) = exp_values[i];
        p_node->GetValue(PRESSURE) = exp_values[i];
        p_elem->GetValue(TEMPERATURE) = exp_values[i];
    }

    KRATOS_EXPECT_EQ(kratos_model_part.NumberOfNodes(), num_nodes);
    KRATOS_EXPECT_EQ(kratos_model_part.NumberOfElements(), num_nodes);
    KRATOS_EXPECT_EQ(kratos_model_part.NumberOfProperties(), 1);

    {std::vector<double> values;
    CoSimIOConversionUtilities::GetData(kratos_model_part, values, AUX_INDEX, DataLocation::NodeHistorical);
    KRATOS_EXPECT_VECTOR_EQ(exp_values, values);}
    {std::vector<double> values;
    CoSimIOConversionUtilities::GetData(kratos_model_part, values, PRESSURE, DataLocation::NodeNonHistorical);
    KRATOS_EXPECT_VECTOR_EQ(exp_values, values);}
    {std::vector<double> values;
    CoSimIOConversionUtilities::GetData(kratos_model_part, values, TEMPERATURE, DataLocation::Element);
    KRATOS_EXPECT_VECTOR_EQ(exp_values, values);}
}

KRATOS_TEST_CASE_IN_SUITE(CoSimIODataExport_reordered_scalar, KratosCosimulationFastSuite)
{
    CoSimIO::ModelPart co_sim_io_model_part("co_sim_io_mp");

    Model model;
    auto& kratos_model_part = model.CreateModelPart("kratos_mp");
    kratos_model_part.AddNodalSolutionStepVariable(AUX_INDEX);

    constexpr std::size_t num_nodes = 5;
    const std::vector<double> exp_values {2,6,8,4,1};
    const std::vector<int> cosimio_ids {1,5,2,6,3};

    for (std::size_t i=0; i<num_nodes; ++i) {
        co_sim_io_model_part.CreateNewNode(cosimio_ids[i], i*1.5, i+3.5, i-8.6);
        co_sim_io_model_part.CreateNewElement(cosimio_ids[i], CoSimIO::ElementType::Point2D, {cosimio_ids[i]});
    }

    const auto& r_serial_data_comm = ParallelEnvironment::GetDataCommunicator("Serial");
    CoSimIOConversionUtilities::CoSimIOModelPartToKratosModelPart(co_sim_io_model_part, kratos_model_part, r_serial_data_comm);

    // force a sort which causes the issues in MPI (and sometimes happens randomly)
    kratos_model_part.Nodes().Sort();
    kratos_model_part.Elements().Sort();

    KRATOS_EXPECT_EQ(kratos_model_part.NumberOfNodes(), num_nodes);
    KRATOS_EXPECT_EQ(kratos_model_part.NumberOfElements(), num_nodes);
    KRATOS_EXPECT_EQ(kratos_model_part.NumberOfProperties(), 1);

    for (std::size_t i=0; i<num_nodes; ++i) {
        kratos_model_part.Nodes().find(cosimio_ids[i])->FastGetSolutionStepValue(AUX_INDEX) = exp_values[i];
        kratos_model_part.Nodes().find(cosimio_ids[i])->GetValue(PRESSURE) = exp_values[i];
        kratos_model_part.Elements().find(cosimio_ids[i])->GetValue(TEMPERATURE) = exp_values[i];
    }

    {std::vector<double> values;
    CoSimIOConversionUtilities::GetData(kratos_model_part, values, AUX_INDEX, DataLocation::NodeHistorical);
    KRATOS_EXPECT_VECTOR_EQ(exp_values, values);}
    {std::vector<double> values;
    CoSimIOConversionUtilities::GetData(kratos_model_part, values, PRESSURE, DataLocation::NodeNonHistorical);
    KRATOS_EXPECT_VECTOR_EQ(exp_values, values);}
    {std::vector<double> values;
    CoSimIOConversionUtilities::GetData(kratos_model_part, values, TEMPERATURE, DataLocation::Element);
    KRATOS_EXPECT_VECTOR_EQ(exp_values, values);}
}

KRATOS_TEST_CASE_IN_SUITE(CoSimIODataExport_direct_vector, KratosCosimulationFastSuite)
{
    Model model;
    auto& kratos_model_part = model.CreateModelPart("kratos_mp");
    kratos_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

    constexpr std::size_t num_nodes = 5;
    const std::vector<double> ref_values {2,6,8,4,1.14};
    std::vector<double> exp_values(ref_values.size()*3);
    for (std::size_t i=0; i<num_nodes; ++i) {
        for (std::size_t j=0; j<3; ++j) {
            exp_values[i*3+j] = ref_values[i];
        }
    }

    auto p_props = kratos_model_part.CreateNewProperties(0);

    for (std::size_t i=0; i<num_nodes; ++i) {
        auto p_node = kratos_model_part.CreateNewNode(i+1, i*1.5, i+3.5, i-8.6);
        auto p_elem = kratos_model_part.CreateNewElement("Element2D1N", i+1, std::vector<std::size_t>{i+1ul}, p_props);
        const array_1d<double, 3> vals(3, ref_values[i]);
        p_node->FastGetSolutionStepValue(DISPLACEMENT) = vals;
        p_node->GetValue(ROTATION) = vals;
        p_elem->GetValue(VELOCITY) = vals;
    }

    KRATOS_EXPECT_EQ(kratos_model_part.NumberOfNodes(), num_nodes);
    KRATOS_EXPECT_EQ(kratos_model_part.NumberOfElements(), num_nodes);
    KRATOS_EXPECT_EQ(kratos_model_part.NumberOfProperties(), 1);

    {std::vector<double> values;
    CoSimIOConversionUtilities::GetData(kratos_model_part, values, DISPLACEMENT, DataLocation::NodeHistorical);
    KRATOS_EXPECT_VECTOR_EQ(exp_values, values);}
    {std::vector<double> values;
    CoSimIOConversionUtilities::GetData(kratos_model_part, values, ROTATION, DataLocation::NodeNonHistorical);
    KRATOS_EXPECT_VECTOR_EQ(exp_values, values);}
    {std::vector<double> values;
    CoSimIOConversionUtilities::GetData(kratos_model_part, values, VELOCITY, DataLocation::Element);
    KRATOS_EXPECT_VECTOR_EQ(exp_values, values);}
}

KRATOS_TEST_CASE_IN_SUITE(CoSimIODataExport_reordered_vector, KratosCosimulationFastSuite)
{
    CoSimIO::ModelPart co_sim_io_model_part("co_sim_io_mp");

    Model model;
    auto& kratos_model_part = model.CreateModelPart("kratos_mp");
    kratos_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

    constexpr std::size_t num_nodes = 5;
    const std::vector<double> ref_values {2,6,8,4,1.14};
    std::vector<double> exp_values(ref_values.size()*3);
    for (std::size_t i=0; i<num_nodes; ++i) {
        for (std::size_t j=0; j<3; ++j) {
            exp_values[i*3+j] = ref_values[i];
        }
    }
    const std::vector<int> cosimio_ids {1,5,2,6,3};

    for (std::size_t i=0; i<num_nodes; ++i) {
        co_sim_io_model_part.CreateNewNode(cosimio_ids[i], i*1.5, i+3.5, i-8.6);
        co_sim_io_model_part.CreateNewElement(cosimio_ids[i], CoSimIO::ElementType::Point2D, {cosimio_ids[i]});
    }

    const auto& r_serial_data_comm = ParallelEnvironment::GetDataCommunicator("Serial");
    CoSimIOConversionUtilities::CoSimIOModelPartToKratosModelPart(co_sim_io_model_part, kratos_model_part, r_serial_data_comm);

    // force a sort which causes the issues in MPI (and sometimes happens randomly)
    kratos_model_part.Nodes().Sort();
    kratos_model_part.Elements().Sort();

    KRATOS_EXPECT_EQ(kratos_model_part.NumberOfNodes(), num_nodes);
    KRATOS_EXPECT_EQ(kratos_model_part.NumberOfElements(), num_nodes);
    KRATOS_EXPECT_EQ(kratos_model_part.NumberOfProperties(), 1);

    for (std::size_t i=0; i<num_nodes; ++i) {
        const array_1d<double, 3> vals(3, ref_values[i]);
        kratos_model_part.Nodes().find(cosimio_ids[i])->FastGetSolutionStepValue(DISPLACEMENT) = vals;
        kratos_model_part.Nodes().find(cosimio_ids[i])->GetValue(ROTATION) = vals;
        kratos_model_part.Elements().find(cosimio_ids[i])->GetValue(VELOCITY) = vals;
    }

    {std::vector<double> values;
    CoSimIOConversionUtilities::GetData(kratos_model_part, values, DISPLACEMENT, DataLocation::NodeHistorical);
    KRATOS_EXPECT_VECTOR_EQ(exp_values, values);}
    {std::vector<double> values;
    CoSimIOConversionUtilities::GetData(kratos_model_part, values, ROTATION, DataLocation::NodeNonHistorical);
    KRATOS_EXPECT_VECTOR_EQ(exp_values, values);}
    {std::vector<double> values;
    CoSimIOConversionUtilities::GetData(kratos_model_part, values, VELOCITY, DataLocation::Element);
    KRATOS_EXPECT_VECTOR_EQ(exp_values, values);}
}


KRATOS_TEST_CASE_IN_SUITE(CoSimIODataImport_direct_scalar, KratosCosimulationFastSuite)
{
    Model model;
    auto& kratos_model_part = model.CreateModelPart("kratos_mp");
    kratos_model_part.AddNodalSolutionStepVariable(AUX_INDEX);

    constexpr std::size_t num_nodes = 5;
    const std::vector<double> ref_values {2,6,8,4,1};

    auto p_props = kratos_model_part.CreateNewProperties(0);

    for (std::size_t i=0; i<num_nodes; ++i) {
        auto p_node = kratos_model_part.CreateNewNode(i+1, i*1.5, i+3.5, i-8.6);
        auto p_elem = kratos_model_part.CreateNewElement("Element2D1N", i+1, std::vector<std::size_t>{i+1ul}, p_props);
    }

    KRATOS_EXPECT_EQ(kratos_model_part.NumberOfNodes(), num_nodes);
    KRATOS_EXPECT_EQ(kratos_model_part.NumberOfElements(), num_nodes);
    KRATOS_EXPECT_EQ(kratos_model_part.NumberOfProperties(), 1);

    CoSimIOConversionUtilities::SetData(kratos_model_part, ref_values, AUX_INDEX, DataLocation::NodeHistorical);
    CoSimIOConversionUtilities::SetData(kratos_model_part, ref_values, PRESSURE, DataLocation::NodeNonHistorical);
    CoSimIOConversionUtilities::SetData(kratos_model_part, ref_values, TEMPERATURE, DataLocation::Element);

    for (std::size_t i=0; i<num_nodes; ++i) {
        auto p_node = kratos_model_part.NodesBegin()+i;
        auto p_elem = kratos_model_part.ElementsBegin()+i;
        KRATOS_EXPECT_DOUBLE_EQ(p_node->FastGetSolutionStepValue(AUX_INDEX), ref_values[i]);
        KRATOS_EXPECT_DOUBLE_EQ(p_node->GetValue(PRESSURE), ref_values[i]);
        KRATOS_EXPECT_DOUBLE_EQ(p_elem->GetValue(TEMPERATURE), ref_values[i]);
    }
}

KRATOS_TEST_CASE_IN_SUITE(CoSimIODataImport_reordered_scalar, KratosCosimulationFastSuite)
{
    CoSimIO::ModelPart co_sim_io_model_part("co_sim_io_mp");

    Model model;
    auto& kratos_model_part = model.CreateModelPart("kratos_mp");
    kratos_model_part.AddNodalSolutionStepVariable(AUX_INDEX);

    constexpr std::size_t num_nodes = 5;
    const std::vector<double> ref_values {2,6,8,4,1.14};
    const std::vector<int> cosimio_ids {1,5,2,6,3};

    for (std::size_t i=0; i<num_nodes; ++i) {
        co_sim_io_model_part.CreateNewNode(cosimio_ids[i], i*1.5, i+3.5, i-8.6);
        co_sim_io_model_part.CreateNewElement(cosimio_ids[i], CoSimIO::ElementType::Point2D, {cosimio_ids[i]});
    }

    const auto& r_serial_data_comm = ParallelEnvironment::GetDataCommunicator("Serial");
    CoSimIOConversionUtilities::CoSimIOModelPartToKratosModelPart(co_sim_io_model_part, kratos_model_part, r_serial_data_comm);

    // force a sort which causes the issues in MPI (and sometimes happens randomly)
    kratos_model_part.Nodes().Sort();
    kratos_model_part.Elements().Sort();

    KRATOS_EXPECT_EQ(kratos_model_part.NumberOfNodes(), num_nodes);
    KRATOS_EXPECT_EQ(kratos_model_part.NumberOfElements(), num_nodes);
    KRATOS_EXPECT_EQ(kratos_model_part.NumberOfProperties(), 1);

    CoSimIOConversionUtilities::SetData(kratos_model_part, ref_values, AUX_INDEX, DataLocation::NodeHistorical);
    CoSimIOConversionUtilities::SetData(kratos_model_part, ref_values, PRESSURE, DataLocation::NodeNonHistorical);
    CoSimIOConversionUtilities::SetData(kratos_model_part, ref_values, TEMPERATURE, DataLocation::Element);

    for (std::size_t i=0; i<num_nodes; ++i) {
        auto& r_node = *(kratos_model_part.Nodes().find(cosimio_ids[i]));
        auto& r_elem = *(kratos_model_part.Elements().find(cosimio_ids[i]));
        KRATOS_EXPECT_DOUBLE_EQ(r_node.FastGetSolutionStepValue(AUX_INDEX), ref_values[i]);
        KRATOS_EXPECT_DOUBLE_EQ(r_node.GetValue(PRESSURE), ref_values[i]);
        KRATOS_EXPECT_DOUBLE_EQ(r_elem.GetValue(TEMPERATURE), ref_values[i]);
    }
}

KRATOS_TEST_CASE_IN_SUITE(CoSimIODataImport_direct_vector, KratosCosimulationFastSuite)
{
    Model model;
    auto& kratos_model_part = model.CreateModelPart("kratos_mp");
    kratos_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

    constexpr std::size_t num_nodes = 5;
    const std::vector<double> ref_values {2,6,8,4,1};
    std::vector<double> set_values(ref_values.size()*3);
    for (std::size_t i=0; i<num_nodes; ++i) {
        for (std::size_t j=0; j<3; ++j) {
            set_values[i*3+j] = ref_values[i];
        }
    }

    auto p_props = kratos_model_part.CreateNewProperties(0);

    for (std::size_t i=0; i<num_nodes; ++i) {
        auto p_node = kratos_model_part.CreateNewNode(i+1, i*1.5, i+3.5, i-8.6);
        auto p_elem = kratos_model_part.CreateNewElement("Element2D1N", i+1, std::vector<std::size_t>{i+1ul}, p_props);
    }

    KRATOS_EXPECT_EQ(kratos_model_part.NumberOfNodes(), num_nodes);
    KRATOS_EXPECT_EQ(kratos_model_part.NumberOfElements(), num_nodes);
    KRATOS_EXPECT_EQ(kratos_model_part.NumberOfProperties(), 1);

    CoSimIOConversionUtilities::SetData(kratos_model_part, set_values, DISPLACEMENT, DataLocation::NodeHistorical);
    CoSimIOConversionUtilities::SetData(kratos_model_part, set_values, ROTATION, DataLocation::NodeNonHistorical);
    CoSimIOConversionUtilities::SetData(kratos_model_part, set_values, VELOCITY, DataLocation::Element);

    for (std::size_t i=0; i<num_nodes; ++i) {
        auto p_node = kratos_model_part.NodesBegin()+i;
        auto p_elem = kratos_model_part.ElementsBegin()+i;
        const array_1d<double, 3> vals(3, ref_values[i]);
        KRATOS_EXPECT_VECTOR_EQ(p_node->FastGetSolutionStepValue(DISPLACEMENT), vals);
        KRATOS_EXPECT_VECTOR_EQ(p_node->GetValue(ROTATION), vals);
        KRATOS_EXPECT_VECTOR_EQ(p_elem->GetValue(VELOCITY), vals);
    }
}

KRATOS_TEST_CASE_IN_SUITE(CoSimIODataImport_reordered_vector, KratosCosimulationFastSuite)
{
    CoSimIO::ModelPart co_sim_io_model_part("co_sim_io_mp");

    Model model;
    auto& kratos_model_part = model.CreateModelPart("kratos_mp");
    kratos_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

    constexpr std::size_t num_nodes = 5;
    const std::vector<double> ref_values {2,6,8,4,1.14};
    std::vector<double> set_values(ref_values.size()*3);
    for (std::size_t i=0; i<num_nodes; ++i) {
        for (std::size_t j=0; j<3; ++j) {
            set_values[i*3+j] = ref_values[i];
        }
    }
    const std::vector<int> cosimio_ids {1,5,2,6,3};

    for (std::size_t i=0; i<num_nodes; ++i) {
        co_sim_io_model_part.CreateNewNode(cosimio_ids[i], i*1.5, i+3.5, i-8.6);
        co_sim_io_model_part.CreateNewElement(cosimio_ids[i], CoSimIO::ElementType::Point2D, {cosimio_ids[i]});
    }

    const auto& r_serial_data_comm = ParallelEnvironment::GetDataCommunicator("Serial");
    CoSimIOConversionUtilities::CoSimIOModelPartToKratosModelPart(co_sim_io_model_part, kratos_model_part, r_serial_data_comm);

    // force a sort which causes the issues in MPI (and sometimes happens randomly)
    kratos_model_part.Nodes().Sort();
    kratos_model_part.Elements().Sort();

    KRATOS_EXPECT_EQ(kratos_model_part.NumberOfNodes(), num_nodes);
    KRATOS_EXPECT_EQ(kratos_model_part.NumberOfElements(), num_nodes);
    KRATOS_EXPECT_EQ(kratos_model_part.NumberOfProperties(), 1);

    CoSimIOConversionUtilities::SetData(kratos_model_part, set_values, DISPLACEMENT, DataLocation::NodeHistorical);
    CoSimIOConversionUtilities::SetData(kratos_model_part, set_values, ROTATION, DataLocation::NodeNonHistorical);
    CoSimIOConversionUtilities::SetData(kratos_model_part, set_values, VELOCITY, DataLocation::Element);

    for (std::size_t i=0; i<num_nodes; ++i) {
        auto& r_node = *(kratos_model_part.Nodes().find(cosimio_ids[i]));
        auto& r_elem = *(kratos_model_part.Elements().find(cosimio_ids[i]));
        const array_1d<double, 3> vals(3, ref_values[i]);
        KRATOS_EXPECT_VECTOR_EQ(r_node.FastGetSolutionStepValue(DISPLACEMENT), vals);
        KRATOS_EXPECT_VECTOR_EQ(r_node.GetValue(ROTATION), vals);
        KRATOS_EXPECT_VECTOR_EQ(r_elem.GetValue(VELOCITY), vals);
    }
}

KRATOS_TEST_CASE_IN_SUITE(CoSimIO_SetGetData_scalar, KratosCosimulationFastSuite)
{
    CoSimIO::ModelPart co_sim_io_model_part("co_sim_io_mp");

    Model model;
    auto& kratos_model_part = model.CreateModelPart("kratos_mp");
    kratos_model_part.AddNodalSolutionStepVariable(AUX_INDEX);

    constexpr std::size_t num_nodes = 5;
    const std::vector<double> exp_values {2,6,8,4,1};
    const std::vector<int> cosimio_ids {1,5,2,6,3};

    for (std::size_t i=0; i<num_nodes; ++i) {
        co_sim_io_model_part.CreateNewNode(cosimio_ids[i], i*1.5, i+3.5, i-8.6);
        co_sim_io_model_part.CreateNewElement(cosimio_ids[i], CoSimIO::ElementType::Point2D, {cosimio_ids[i]});
    }

    const auto& r_serial_data_comm = ParallelEnvironment::GetDataCommunicator("Serial");
    CoSimIOConversionUtilities::CoSimIOModelPartToKratosModelPart(co_sim_io_model_part, kratos_model_part, r_serial_data_comm);

    // force a sort which causes the issues in MPI (and sometimes happens randomly)
    kratos_model_part.Nodes().Sort();
    kratos_model_part.Elements().Sort();

    KRATOS_EXPECT_EQ(kratos_model_part.NumberOfNodes(), num_nodes);
    KRATOS_EXPECT_EQ(kratos_model_part.NumberOfElements(), num_nodes);
    KRATOS_EXPECT_EQ(kratos_model_part.NumberOfProperties(), 1);

    CoSimIOConversionUtilities::SetData(kratos_model_part, exp_values, AUX_INDEX, DataLocation::NodeHistorical);
    CoSimIOConversionUtilities::SetData(kratos_model_part, exp_values, PRESSURE, DataLocation::NodeNonHistorical);
    CoSimIOConversionUtilities::SetData(kratos_model_part, exp_values, TEMPERATURE, DataLocation::Element);

    {std::vector<double> values;
    CoSimIOConversionUtilities::GetData(kratos_model_part, values, AUX_INDEX, DataLocation::NodeHistorical);
    KRATOS_EXPECT_VECTOR_EQ(exp_values, values);}
    {std::vector<double> values;
    CoSimIOConversionUtilities::GetData(kratos_model_part, values, PRESSURE, DataLocation::NodeNonHistorical);
    KRATOS_EXPECT_VECTOR_EQ(exp_values, values);}
    {std::vector<double> values;
    CoSimIOConversionUtilities::GetData(kratos_model_part, values, TEMPERATURE, DataLocation::Element);
    KRATOS_EXPECT_VECTOR_EQ(exp_values, values);}
}

KRATOS_TEST_CASE_IN_SUITE(CoSimIO_SetGetData_vector, KratosCosimulationFastSuite)
{
    CoSimIO::ModelPart co_sim_io_model_part("co_sim_io_mp");

    Model model;
    auto& kratos_model_part = model.CreateModelPart("kratos_mp");
    kratos_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

    constexpr std::size_t num_nodes = 5;
    const std::vector<double> exp_values {2,6,8,4,1,
                                          6,3,7,4,5,
                                          8,3,6,1,4};
    const std::vector<int> cosimio_ids {1,5,2,6,3};

    for (std::size_t i=0; i<num_nodes; ++i) {
        co_sim_io_model_part.CreateNewNode(cosimio_ids[i], i*1.5, i+3.5, i-8.6);
        co_sim_io_model_part.CreateNewElement(cosimio_ids[i], CoSimIO::ElementType::Point2D, {cosimio_ids[i]});
    }

    const auto& r_serial_data_comm = ParallelEnvironment::GetDataCommunicator("Serial");
    CoSimIOConversionUtilities::CoSimIOModelPartToKratosModelPart(co_sim_io_model_part, kratos_model_part, r_serial_data_comm);

    // force a sort which causes the issues in MPI (and sometimes happens randomly)
    kratos_model_part.Nodes().Sort();
    kratos_model_part.Elements().Sort();

    KRATOS_EXPECT_EQ(kratos_model_part.NumberOfNodes(), num_nodes);
    KRATOS_EXPECT_EQ(kratos_model_part.NumberOfElements(), num_nodes);
    KRATOS_EXPECT_EQ(kratos_model_part.NumberOfProperties(), 1);

    CoSimIOConversionUtilities::SetData(kratos_model_part, exp_values, DISPLACEMENT, DataLocation::NodeHistorical);
    CoSimIOConversionUtilities::SetData(kratos_model_part, exp_values, ROTATION, DataLocation::NodeNonHistorical);
    CoSimIOConversionUtilities::SetData(kratos_model_part, exp_values, VELOCITY, DataLocation::Element);

    {std::vector<double> values;
    CoSimIOConversionUtilities::GetData(kratos_model_part, values, DISPLACEMENT, DataLocation::NodeHistorical);
    KRATOS_EXPECT_VECTOR_EQ(exp_values, values);}
    {std::vector<double> values;
    CoSimIOConversionUtilities::GetData(kratos_model_part, values, ROTATION, DataLocation::NodeNonHistorical);
    KRATOS_EXPECT_VECTOR_EQ(exp_values, values);}
    {std::vector<double> values;
    CoSimIOConversionUtilities::GetData(kratos_model_part, values, VELOCITY, DataLocation::Element);
    KRATOS_EXPECT_VECTOR_EQ(exp_values, values);}
}

namespace DistributedTestHelpers {

std::size_t GetPartnerRank()
{
    const auto& r_world_data_comm = ParallelEnvironment::GetDataCommunicator("World");
    const std::size_t my_rank = r_world_data_comm.Rank();
    const std::size_t world_size = r_world_data_comm.Size();

    return (my_rank+1) % world_size;
}

std::size_t GetId(const std::size_t NumLocalNodesPerRank, const std::size_t LocalNodeIndex)
{
    const std::size_t my_rank = ParallelEnvironment::GetDataCommunicator("World").Rank();
    return my_rank*NumLocalNodesPerRank + LocalNodeIndex+1;
}

std::size_t GetGhostId(const std::size_t NumLocalNodesPerRank, const std::size_t LocalNodeIndex)
{
    return GetPartnerRank()*NumLocalNodesPerRank + LocalNodeIndex+1;
}

void CreateDistributedNodes(
    CoSimIO::ModelPart& rCoSimIOModelPart,
    const std::size_t NumLocalNodesPerRank,
    const std::size_t NumGhostNodesPerRank)
{
    KRATOS_EXPECT_GT(NumLocalNodesPerRank, NumGhostNodesPerRank);

    const auto& r_world_data_comm = ParallelEnvironment::GetDataCommunicator("World");
    const std::size_t my_rank = r_world_data_comm.Rank();
    const std::size_t world_size = r_world_data_comm.Size();

    auto create_ghost_nodes = [&](){
        for (std::size_t i=0; i<NumGhostNodesPerRank; ++i) {
            rCoSimIOModelPart.CreateNewGhostNode(GetGhostId(NumLocalNodesPerRank, i), 0,0,0, GetPartnerRank());
        }
    };

    if (GetPartnerRank() < my_rank) {
        // Kratos requires a consecutively ordered input as it sorts the nodes internally and then the matching would be lost
        // in case the ghost nodes have smaller ids than the local nodes they have to be added first as currently
        // the CoSimIO::ModelPart does not allow sorting
        create_ghost_nodes();
    }

    // create local nodes
    for (std::size_t i=0; i<NumLocalNodesPerRank; ++i) {
        rCoSimIOModelPart.CreateNewNode(GetId(NumLocalNodesPerRank, i), 0,0,0);
    }

    if (GetPartnerRank() > my_rank) {
        // in this case the ghost nodes have larger Ids than the local nodes
        create_ghost_nodes();
    }

    KRATOS_EXPECT_EQ(rCoSimIOModelPart.NumberOfNodes(), NumLocalNodesPerRank+NumGhostNodesPerRank);
    KRATOS_EXPECT_EQ(rCoSimIOModelPart.NumberOfLocalNodes(), NumLocalNodesPerRank);
    KRATOS_EXPECT_EQ(rCoSimIOModelPart.NumberOfGhostNodes(), NumGhostNodesPerRank);
    KRATOS_EXPECT_EQ(r_world_data_comm.SumAll(static_cast<int>(rCoSimIOModelPart.NumberOfLocalNodes())), static_cast<int>(NumLocalNodesPerRank*world_size));
    KRATOS_EXPECT_EQ(rCoSimIOModelPart.NumberOfElements(), 0);
}

void CreateDistributedNodesUnordered(
    CoSimIO::ModelPart& rCoSimIOModelPart,
    const std::size_t NumLocalNodesPerRank,
    const std::size_t NumGhostNodesPerRank)
{
    KRATOS_EXPECT_GT(NumLocalNodesPerRank, NumGhostNodesPerRank);

    const auto& r_world_data_comm = ParallelEnvironment::GetDataCommunicator("World");
    const std::size_t world_size = r_world_data_comm.Size();

    for (std::size_t i=0; i<NumLocalNodesPerRank; ++i) {
        rCoSimIOModelPart.CreateNewNode(GetId(NumLocalNodesPerRank, i), 0,0,0);
    }

    for (int i=static_cast<int>(NumGhostNodesPerRank); --i>=0;) { // must use signed type here!
        rCoSimIOModelPart.CreateNewGhostNode(GetGhostId(NumLocalNodesPerRank, i), 0,0,0, GetPartnerRank());
    }

    KRATOS_EXPECT_EQ(rCoSimIOModelPart.NumberOfNodes(), NumLocalNodesPerRank+NumGhostNodesPerRank);
    KRATOS_EXPECT_EQ(rCoSimIOModelPart.NumberOfLocalNodes(), NumLocalNodesPerRank);
    KRATOS_EXPECT_EQ(rCoSimIOModelPart.NumberOfGhostNodes(), NumGhostNodesPerRank);
    KRATOS_EXPECT_EQ(r_world_data_comm.SumAll(static_cast<int>(rCoSimIOModelPart.NumberOfLocalNodes())), static_cast<int>(NumLocalNodesPerRank*world_size));
    KRATOS_EXPECT_EQ(rCoSimIOModelPart.NumberOfElements(), 0);
}

void CreateDistributedNodes(
    Kratos::ModelPart& rKratosModelPart,
    const std::size_t NumLocalNodesPerRank,
    const std::size_t NumGhostNodesPerRank)
{
    KRATOS_EXPECT_GT(NumLocalNodesPerRank, NumGhostNodesPerRank);

    const auto& r_world_data_comm = ParallelEnvironment::GetDataCommunicator("World");
    const std::size_t my_rank = r_world_data_comm.Rank();
    const std::size_t world_size = r_world_data_comm.Size();

    // create local nodes
    for (std::size_t i=0; i<NumLocalNodesPerRank; ++i) {
        auto node = rKratosModelPart.CreateNewNode(GetId(NumLocalNodesPerRank, i), i*10,i/3.58,i-99.14);
        node->FastGetSolutionStepValue(PARTITION_INDEX) = my_rank;
    }

    // create ghost nodes
    for (std::size_t i=0; i<NumGhostNodesPerRank; ++i) {
        auto node = rKratosModelPart.CreateNewNode(GetGhostId(NumLocalNodesPerRank, i), i-10,i*3,i+99.14);
        node->FastGetSolutionStepValue(PARTITION_INDEX) = GetPartnerRank();
    }

    // this calls the ParallelFillCommunicator
    ParallelEnvironment::CreateFillCommunicatorFromGlobalParallelism(rKratosModelPart, r_world_data_comm)->Execute();

    KRATOS_EXPECT_EQ(rKratosModelPart.NumberOfNodes(), NumLocalNodesPerRank+NumGhostNodesPerRank);
    KRATOS_EXPECT_EQ(rKratosModelPart.GetCommunicator().LocalMesh().NumberOfNodes(), NumLocalNodesPerRank);
    KRATOS_EXPECT_EQ(rKratosModelPart.GetCommunicator().GhostMesh().NumberOfNodes(), NumGhostNodesPerRank);
    KRATOS_EXPECT_EQ(rKratosModelPart.GetCommunicator().GlobalNumberOfNodes(), (NumLocalNodesPerRank)*world_size);
    KRATOS_EXPECT_EQ(rKratosModelPart.NumberOfProperties(), 0);
}

} // DistributedTestHelpers

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(DistributedCoSimIOModelPartToKratosModelPart_NodesOnly, KratosCosimulationMPIFastSuite)
{
    const auto& r_world_data_comm = ParallelEnvironment::GetDataCommunicator("World");

    Model model;
    auto& kratos_model_part = model.CreateModelPart("kratos_mp");

    CoSimIO::ModelPart co_sim_io_model_part("co_sim_io_mp");

    constexpr std::size_t num_local_nodes_per_rank = 5;
    constexpr std::size_t num_ghost_nodes_per_rank = 3;

    DistributedTestHelpers::CreateDistributedNodes(co_sim_io_model_part, num_local_nodes_per_rank, num_ghost_nodes_per_rank);

    CoSimIOConversionUtilities::CoSimIOModelPartToKratosModelPart(co_sim_io_model_part, kratos_model_part, r_world_data_comm);

    CheckDistributedModelPartsAreEqual(kratos_model_part, co_sim_io_model_part);
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(DistributedCoSimIOModelPartToKratosModelPart_NodesOnly_Unordered, KratosCosimulationMPIFastSuite)
{
    const auto& r_world_data_comm = ParallelEnvironment::GetDataCommunicator("World");

    Model model;
    auto& kratos_model_part = model.CreateModelPart("kratos_mp");

    CoSimIO::ModelPart co_sim_io_model_part("co_sim_io_mp");

    constexpr std::size_t num_local_nodes_per_rank = 5;
    constexpr std::size_t num_ghost_nodes_per_rank = 3;

    DistributedTestHelpers::CreateDistributedNodesUnordered(co_sim_io_model_part, num_local_nodes_per_rank, num_ghost_nodes_per_rank);

    CoSimIOConversionUtilities::CoSimIOModelPartToKratosModelPart(co_sim_io_model_part, kratos_model_part, r_world_data_comm);

    CheckDistributedModelPartsAreEqualButEntitiesAreOrderedDifferently(kratos_model_part, co_sim_io_model_part);
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(DistributedKratosModelPartToCoSimIOModelPart_NodesOnly, KratosCosimulationMPIFastSuite)
{
    Model model;
    auto& kratos_model_part = model.CreateModelPart("kratos_mp");
    kratos_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

    CoSimIO::ModelPart co_sim_io_model_part("co_sim_io_mp");

    constexpr std::size_t num_local_nodes_per_rank = 5;
    constexpr std::size_t num_ghost_nodes_per_rank = 3;

    DistributedTestHelpers::CreateDistributedNodes(kratos_model_part, num_local_nodes_per_rank, num_ghost_nodes_per_rank);

    CoSimIOConversionUtilities::KratosModelPartToCoSimIOModelPart(kratos_model_part, co_sim_io_model_part);

    CheckDistributedModelPartsAreEqual(kratos_model_part, co_sim_io_model_part);
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(DistributedCoSimIOModelPartToKratosModelPart, KratosCosimulationMPIFastSuite)
{
    const auto& r_world_data_comm = ParallelEnvironment::GetDataCommunicator("World");

    Model model;
    auto& kratos_model_part = model.CreateModelPart("kratos_mp");

    CoSimIO::ModelPart co_sim_io_model_part("co_sim_io_mp");

    constexpr std::size_t num_local_nodes_per_rank = 5;
    constexpr std::size_t num_ghost_nodes_per_rank = 3;

    DistributedTestHelpers::CreateDistributedNodes(co_sim_io_model_part, num_local_nodes_per_rank, num_ghost_nodes_per_rank);

    // elements that use two local nodes
    for (std::size_t i=0; i<num_ghost_nodes_per_rank; ++i) {
        CoSimIO::ConnectivitiesType conn {
            static_cast<int>(DistributedTestHelpers::GetId(num_local_nodes_per_rank, i)),
            static_cast<int>(DistributedTestHelpers::GetId(num_local_nodes_per_rank, i+1))};
        co_sim_io_model_part.CreateNewElement(DistributedTestHelpers::GetId(num_local_nodes_per_rank, i), CoSimIO::ElementType::Line2D2, conn);
    }

    // elements that use one local and one ghost node
    for (std::size_t i=0; i<num_ghost_nodes_per_rank; ++i) {
        CoSimIO::ConnectivitiesType conn {
            static_cast<int>(DistributedTestHelpers::GetId(num_local_nodes_per_rank, i)),
            static_cast<int>(DistributedTestHelpers::GetGhostId(num_local_nodes_per_rank, i))};
        co_sim_io_model_part.CreateNewElement(DistributedTestHelpers::GetId(num_local_nodes_per_rank, i)+num_ghost_nodes_per_rank, CoSimIO::ElementType::Line2D2, conn);
    }

    CoSimIOConversionUtilities::CoSimIOModelPartToKratosModelPart(co_sim_io_model_part, kratos_model_part, r_world_data_comm);

    CheckDistributedModelPartsAreEqual(kratos_model_part, co_sim_io_model_part);
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(DistributedCoSimIOModelPartToKratosModelPart_Unordered, KratosCosimulationMPIFastSuite)
{
    const auto& r_world_data_comm = ParallelEnvironment::GetDataCommunicator("World");

    Model model;
    auto& kratos_model_part = model.CreateModelPart("kratos_mp");

    CoSimIO::ModelPart co_sim_io_model_part("co_sim_io_mp");

    constexpr std::size_t num_local_nodes_per_rank = 5;
    constexpr std::size_t num_ghost_nodes_per_rank = 3;

    DistributedTestHelpers::CreateDistributedNodesUnordered(co_sim_io_model_part, num_local_nodes_per_rank, num_ghost_nodes_per_rank);

    // elements that use two local nodes
    for (std::size_t i=0; i<num_ghost_nodes_per_rank; ++i) {
        CoSimIO::ConnectivitiesType conn {
            static_cast<int>(DistributedTestHelpers::GetId(num_local_nodes_per_rank, i)),
            static_cast<int>(DistributedTestHelpers::GetId(num_local_nodes_per_rank, i+1))};
        co_sim_io_model_part.CreateNewElement(DistributedTestHelpers::GetId(num_local_nodes_per_rank, i), CoSimIO::ElementType::Line2D2, conn);
    }

    // elements that use one local and one ghost node
    for (std::size_t i=0; i<num_ghost_nodes_per_rank; ++i) {
        CoSimIO::ConnectivitiesType conn {
            static_cast<int>(DistributedTestHelpers::GetId(num_local_nodes_per_rank, i)),
            static_cast<int>(DistributedTestHelpers::GetGhostId(num_local_nodes_per_rank, i))};
        co_sim_io_model_part.CreateNewElement(DistributedTestHelpers::GetId(num_local_nodes_per_rank, i)+num_ghost_nodes_per_rank, CoSimIO::ElementType::Line2D2, conn);
    }

    CoSimIOConversionUtilities::CoSimIOModelPartToKratosModelPart(co_sim_io_model_part, kratos_model_part, r_world_data_comm);

    CheckDistributedModelPartsAreEqualButEntitiesAreOrderedDifferently(kratos_model_part, co_sim_io_model_part);
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(KratosDistributedModelPartToCoSimIOModelPart, KratosCosimulationMPIFastSuite)
{
    const auto& r_world_data_comm = ParallelEnvironment::GetDataCommunicator("World");

    Model model;
    auto& kratos_model_part = model.CreateModelPart("kratos_mp");
    kratos_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

    CoSimIO::ModelPart co_sim_io_model_part("co_sim_io_mp");

    constexpr std::size_t num_local_nodes_per_rank = 5;
    constexpr std::size_t num_ghost_nodes_per_rank = 3;

    DistributedTestHelpers::CreateDistributedNodes(kratos_model_part, num_local_nodes_per_rank, num_ghost_nodes_per_rank);

    Properties::Pointer p_props(kratos_model_part.CreateNewProperties(0));

    // elements that use two local nodes
    for (std::size_t i=0; i<num_ghost_nodes_per_rank; ++i) {
        std::vector<ModelPart::IndexType> conn {
            static_cast<ModelPart::IndexType>(DistributedTestHelpers::GetId(num_local_nodes_per_rank, i)),
            static_cast<ModelPart::IndexType>(DistributedTestHelpers::GetId(num_local_nodes_per_rank, i+1))};
        const std::size_t elem_id = DistributedTestHelpers::GetId(num_local_nodes_per_rank, i);
        kratos_model_part.CreateNewElement("Element2D2N", elem_id, conn, p_props);
    }

    // elements that use one local and one ghost node
    for (std::size_t i=0; i<num_ghost_nodes_per_rank; ++i) {
        std::vector<ModelPart::IndexType> conn {
            static_cast<ModelPart::IndexType>(DistributedTestHelpers::GetId(num_local_nodes_per_rank, i)),
            static_cast<ModelPart::IndexType>(DistributedTestHelpers::GetGhostId(num_local_nodes_per_rank, i))};
        const std::size_t elem_id = DistributedTestHelpers::GetId(num_local_nodes_per_rank, i)+num_ghost_nodes_per_rank;
        kratos_model_part.CreateNewElement("Element2D2N", elem_id, conn, p_props);
    }

    // this calls the ParallelFillCommunicator
    ParallelEnvironment::CreateFillCommunicatorFromGlobalParallelism(kratos_model_part, r_world_data_comm)->Execute();

    CoSimIOConversionUtilities::KratosModelPartToCoSimIOModelPart(kratos_model_part, co_sim_io_model_part);

    CheckDistributedModelPartsAreEqual(kratos_model_part, co_sim_io_model_part);
}

} // namespace Testing
} // namespace Kratos
