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
#include "custom_external_libraries/CoSimIO/co_sim_io/co_sim_io.hpp"

// Project includes
#include "testing/testing.h"
#include "includes/model_part.h"


namespace Kratos {
namespace Testing {

void CheckNodesAreEqual(
    const Kratos::Node<3>& rKratosNode,
    const CoSimIO::Node& rCoSimIONode);

void CheckElementsAreEqual(
    const Kratos::Element& rKratosElement,
    const CoSimIO::Element& rCoSimIOElement);

void CheckModelPartsAreEqual(
    const Kratos::ModelPart& rKratosModelPart,
    const CoSimIO::ModelPart& rCoSimIOModelPart);

void CheckDistributedModelPartsAreEqual(
    const Kratos::ModelPart& rKratosModelPart,
    const CoSimIO::ModelPart& rCoSimIOModelPart);

} // namespace Kratos
} // namespace Testing