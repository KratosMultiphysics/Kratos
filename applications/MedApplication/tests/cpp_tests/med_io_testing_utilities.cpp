// KRATOS  __  __          _    _                _ _           _   _
//        |  \/  | ___  __| |  / \   _ __  _ __ | (_) ___ __ _| |_(_) ___  _ ___
//        | |\/| |/ _ \/ _` | / _ \ | '_ \| '_ \| | |/ __/ _` | __| |/ _ \| '_  |
//        | |  | |  __/ (_| |/ ___ \| |_) | |_) | | | (_| (_| | |_| | (_) | | | |
//        |_|  |_|\___|\__,_/_/   \_\ .__/| .__/|_|_|\___\__,_|\__|_|\___/|_| |_|
//                                  |_|   |_|
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

// System includes
#include <algorithm>

// External includes

// Project includes
#include "med_io_testing_utilities.h"


namespace Kratos::Testing {

namespace { // helpers namespace

void CheckEntitiesAreEqual(
    const Kratos::Node<3>& rNode1,
    const Kratos::Node<3>& rNode2)
{
    KRATOS_TRY

    KRATOS_CHECK_EQUAL(rNode1.Id(), rNode2.Id());

    KRATOS_CHECK_DOUBLE_EQUAL(rNode1.X(),  rNode2.X());
    KRATOS_CHECK_DOUBLE_EQUAL(rNode1.X0(), rNode2.X0());

    KRATOS_CHECK_DOUBLE_EQUAL(rNode1.Y(),  rNode2.Y());
    KRATOS_CHECK_DOUBLE_EQUAL(rNode1.Y0(), rNode2.Y0());

    KRATOS_CHECK_DOUBLE_EQUAL(rNode1.Z(),  rNode2.Z());
    KRATOS_CHECK_DOUBLE_EQUAL(rNode1.Z0(), rNode2.Z0());

    KRATOS_CATCH("")
}

template<class T>
bool contains(
    const std::vector<T> rVec,
    const T& rValue)
{
    return std::find(rVec.begin(), rVec.end(), rValue) != rVec.end();
}

} // helpers namespace

void CheckNodesAreEqual(
    const Kratos::Node<3>& rNode1,
    const Kratos::Node<3>& rNode2)
{
    KRATOS_TRY

    CheckEntitiesAreEqual(rNode1, rNode2);

    KRATOS_CATCH("")
}

void CheckModelPartsAreEqual(
    const Kratos::ModelPart& rModelPart1,
    const Kratos::ModelPart& rModelPart2)
{
    KRATOS_TRY

    // basic checks
    KRATOS_CHECK_EQUAL(rModelPart1.NumberOfNodes(), rModelPart2.NumberOfNodes());

    // check nodes
    for (std::size_t i=0; i<rModelPart1.NumberOfNodes(); ++i) {
        CheckNodesAreEqual(*(rModelPart1.NodesBegin()+i), *(rModelPart2.NodesBegin()+i));
    }

    KRATOS_CHECK_EQUAL(rModelPart1.NumberOfSubModelParts(), rModelPart2.NumberOfSubModelParts());

    const auto& r_smp2_names = rModelPart2.GetSubModelPartNames();

    for (const auto& r_smp_name : rModelPart1.GetSubModelPartNames()) {
        KRATOS_CHECK(contains(r_smp2_names, r_smp_name));
        // CheckModelPartsAreEqual(rModelPart1.GetSubModelPart(r_smp_name), rModelPart1.GetSubModelPart(r_smp_name)); // TODO requires const version of "GetSubModelPart"
    }

    KRATOS_CATCH("")
}

} // namespace Kratos::Testing
