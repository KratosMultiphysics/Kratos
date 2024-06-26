//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos A. Roig
//                   Jordi Cotela Dalmau
//
//

// System includes

// External includes

// Project includes
#include "tests/test_utilities/test_suite.h"

/*
 * Suite for the mpi testing environment (mKernel(true))
*/
class KratosMPICoreFastSuite : public KratosCoreFastSuite
{
    public:
        void TearDown() override;

    protected:
        KratosMPICoreFastSuite(): mKernel(true) {}
        ~KratosMPICoreFastSuite() {}

    	Kratos::Kernel mKernel;
};