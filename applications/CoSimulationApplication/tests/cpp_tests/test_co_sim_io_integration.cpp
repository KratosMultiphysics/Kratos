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
#include "containers/model.h"
#include "testing/testing.h"
#include "custom_utilities/co_sim_io_conversion_utilities.h"

// Utilities

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(CoSimIOModelPartToKratosModelPart_NodesOnly, KratosCosimulationFastSuite)
{
    Model model;
    auto& kratos_model_part = model.CreateModelPart("kratos_mp");

    CoSimIO::ModelPart co_sim_io_model_part("co_sim_io_mp");

    constexpr std::size_t num_nodes = 5;

    for (std::size_t i=0; i<num_nodes; ++i) {
        co_sim_io_model_part.CreateNewNode(i+1, i*1.5, i+3.5, i-8.6);
    }

    KRATOS_CHECK_EQUAL(co_sim_io_model_part.NumberOfNodes(), num_nodes);
    KRATOS_CHECK_EQUAL(co_sim_io_model_part.NumberOfLocalNodes(), num_nodes);
    KRATOS_CHECK_EQUAL(co_sim_io_model_part.NumberOfGhostNodes(), 0);
    KRATOS_CHECK_EQUAL(co_sim_io_model_part.NumberOfElements(), 0);

    CoSimIOConversionUtilities::CoSimIOModelPartToKratosModelPart(co_sim_io_model_part, kratos_model_part);

    KRATOS_CHECK_EQUAL(kratos_model_part.NumberOfNodes(), num_nodes);
    KRATOS_CHECK_EQUAL(kratos_model_part.NumberOfElements(), 0);
    KRATOS_CHECK_EQUAL(kratos_model_part.NumberOfProperties(), 0);

    const auto nodes_begin = kratos_model_part.NodesBegin();

    for (std::size_t i=0; i<num_nodes; ++i) {
        const std::size_t exp_id = i+1;

        const auto& r_node = *(nodes_begin+i);

        const auto& r_co_sim_io_node = co_sim_io_model_part.GetNode(exp_id);

        KRATOS_CHECK_EQUAL(r_node.Id(), exp_id);
        KRATOS_CHECK_DOUBLE_EQUAL(r_node.X(), r_co_sim_io_node.X());
        KRATOS_CHECK_DOUBLE_EQUAL(r_node.X0(), r_co_sim_io_node.X());
        KRATOS_CHECK_DOUBLE_EQUAL(r_node.Y(), r_co_sim_io_node.Y());
        KRATOS_CHECK_DOUBLE_EQUAL(r_node.Y0(), r_co_sim_io_node.Y());
        KRATOS_CHECK_DOUBLE_EQUAL(r_node.Z(), r_co_sim_io_node.Z());
        KRATOS_CHECK_DOUBLE_EQUAL(r_node.Z0(), r_co_sim_io_node.Z());
    }
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

    KRATOS_CHECK_EQUAL(kratos_model_part.NumberOfNodes(), num_nodes);
    KRATOS_CHECK_EQUAL(kratos_model_part.NumberOfProperties(), 0);

    CoSimIOConversionUtilities::KratosModelPartToCoSimIOModelPart(kratos_model_part, co_sim_io_model_part);

    KRATOS_CHECK_EQUAL(co_sim_io_model_part.NumberOfNodes(), num_nodes);
    KRATOS_CHECK_EQUAL(co_sim_io_model_part.NumberOfLocalNodes(), num_nodes);
    KRATOS_CHECK_EQUAL(co_sim_io_model_part.NumberOfGhostNodes(), 0);
    KRATOS_CHECK_EQUAL(co_sim_io_model_part.NumberOfElements(), 0);

    const auto nodes_begin = co_sim_io_model_part.NodesBegin();

    for (std::size_t i=0; i<num_nodes; ++i) {
        const int exp_id = i+1;

        const auto& r_co_sim_io_node = **(nodes_begin+i);

        const auto& r_kratos_node = kratos_model_part.GetNode(exp_id);

        KRATOS_CHECK_EQUAL(r_co_sim_io_node.Id(), exp_id);
        KRATOS_CHECK_DOUBLE_EQUAL(r_kratos_node.X(), r_co_sim_io_node.X());
        KRATOS_CHECK_DOUBLE_EQUAL(r_kratos_node.X0(), r_co_sim_io_node.X());
        KRATOS_CHECK_DOUBLE_EQUAL(r_kratos_node.Y(), r_co_sim_io_node.Y());
        KRATOS_CHECK_DOUBLE_EQUAL(r_kratos_node.Y0(), r_co_sim_io_node.Y());
        KRATOS_CHECK_DOUBLE_EQUAL(r_kratos_node.Z(), r_co_sim_io_node.Z());
        KRATOS_CHECK_DOUBLE_EQUAL(r_kratos_node.Z0(), r_co_sim_io_node.Z());
    }
}

KRATOS_TEST_CASE_IN_SUITE(CoSimIOModelPartToKratosModelPart, KratosCosimulationFastSuite)
{

}

KRATOS_TEST_CASE_IN_SUITE(KratosModelPartToCoSimIOModelPart, KratosCosimulationFastSuite)
{

}



} // namespace Testing
} // namespace Kratos
