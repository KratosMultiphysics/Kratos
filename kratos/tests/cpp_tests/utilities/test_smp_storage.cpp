//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <vector>

// Project includes
#include "utilities/reduction_utilities.h"
#include "utilities/parallel_utilities.h"
#include "utilities/smp_storage.h"
#include "testing/testing.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(SMPStorage, KratosCoreFastSuite)
{
    class Scheme
    {
    public:
        void CalculateLocalSystem(int& LHS, int ProcessInfo) {
            auto& tls = mSMPStorage.GetThreadLocalStorage();
            tls.b = ProcessInfo * 3;
            LHS = tls.b;
        }
    private:
        struct TLS
        {
            int b;
        };

        SMPStorage<TLS> mSMPStorage;
    };

    Scheme scheme;

    int n = 1e+4;
    const int result = IndexPartition<int>(n).for_each<SumReduction<int>>([&](const int ProcessInfo) -> int {
        int lhs;
        scheme.CalculateLocalSystem(lhs, ProcessInfo);
        return lhs;
    });

    KRATOS_CHECK_EQUAL(result, n * (n-1) * 3 / 2);
}

}
}