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
#include "co_simulation_testing_utilities.h"


namespace Kratos {
namespace Testing {

void CheckNodesAreEqual(
    const Kratos::Node<3>& rKratosNode,
    const CoSimIO::Node& rCoSimIONode)
{
    KRATOS_TRY

    KRATOS_CHECK_EQUAL(rKratosNode.Id(), static_cast<std::size_t>(rCoSimIONode.Id()));

    KRATOS_CHECK_DOUBLE_EQUAL(rKratosNode.X(),  rCoSimIONode.X());
    KRATOS_CHECK_DOUBLE_EQUAL(rKratosNode.X0(), rCoSimIONode.X());

    KRATOS_CHECK_DOUBLE_EQUAL(rKratosNode.Y(),  rCoSimIONode.Y());
    KRATOS_CHECK_DOUBLE_EQUAL(rKratosNode.Y0(), rCoSimIONode.Y());

    KRATOS_CHECK_DOUBLE_EQUAL(rKratosNode.Z(),  rCoSimIONode.Z());
    KRATOS_CHECK_DOUBLE_EQUAL(rKratosNode.Z0(), rCoSimIONode.Z());

    KRATOS_CATCH("")
}

void CheckElementsAreEqual(
    const Kratos::Element& rKratosElement,
    const CoSimIO::Element& rCoSimIOElement)
{
    KRATOS_TRY

    // basic checks
    KRATOS_CHECK_EQUAL(rKratosElement.Id(), static_cast<std::size_t>(rCoSimIOElement.Id()));
    // KRATOS_CHECK_EQUAL(rKratosElement.Type(), rCoSimIOElement.Type());
    KRATOS_CHECK_EQUAL(rKratosElement.GetGeometry().PointsNumber(), rCoSimIOElement.NumberOfNodes());

    // check nodes
    for (std::size_t i=0; i<rCoSimIOElement.NumberOfNodes(); ++i) {
        CheckNodesAreEqual(*(rKratosElement.GetGeometry().begin()+i), **(rCoSimIOElement.NodesBegin()+i));
    }

    KRATOS_CATCH("")
}

void CheckModelPartsAreEqual(
    const Kratos::ModelPart& rKratosModelPart,
    const CoSimIO::ModelPart& rCoSimIOModelPart)
{
    KRATOS_TRY

    // basic checks
    KRATOS_CHECK_EQUAL(rKratosModelPart.NumberOfNodes(),      rCoSimIOModelPart.NumberOfNodes());
    KRATOS_CHECK_EQUAL(rKratosModelPart.NumberOfElements(),   rCoSimIOModelPart.NumberOfElements());

    // check nodes
    for (std::size_t i=0; i<rCoSimIOModelPart.NumberOfNodes(); ++i) {
        CheckNodesAreEqual(*(rKratosModelPart.NodesBegin()+i), **(rCoSimIOModelPart.NodesBegin()+i));
    }

    // check elements
    for (std::size_t i=0; i<rCoSimIOModelPart.NumberOfElements(); ++i) {
        CheckElementsAreEqual(*(rKratosModelPart.ElementsBegin()+i), **(rCoSimIOModelPart.ElementsBegin()+i));
    }

    KRATOS_CATCH("")
}

} // namespace Kratos
} // namespace Testing
