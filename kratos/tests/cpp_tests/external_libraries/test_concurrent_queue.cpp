//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Suneth Warnakulasuriya
//
//

// System includes

// External includes
#include "concurrentqueue/concurrentqueue.h"

// Project includes
#include "testing/testing.h"
#include "utilities/parallel_utilities.h"
#include "includes/checks.h"


namespace Kratos {

namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(ConcurrentQueue, KratosExternalLibrariesFastSuite)
{
    moodycamel::ConcurrentQueue<int> concurrent_queue;

    int n = 1e+4;
    IndexPartition<int>(n).for_each([&](const int Index) {
        concurrent_queue.enqueue(Index);
    });

    bool found_value = true;
    int sum = 0;
    while (found_value) {
        int value = 0;
        found_value = concurrent_queue.try_dequeue(value);
        sum += value;
    }

    KRATOS_CHECK_EQUAL(n * (n - 1) / 2, sum);
}

} // namespace Testing.
} // namespace Kratos.
