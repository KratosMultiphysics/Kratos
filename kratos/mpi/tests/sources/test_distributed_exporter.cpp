//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi

#include <utility>
#include <iostream>
#include <random>

#include "testing/testing.h"
#include "mpi/includes/mpi_data_communicator.h"
#include "containers/distributed_system_vector.h"
#include "containers/distributed_vector_exporter.h"



namespace Kratos {
namespace Testing {

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(DistributedVectorExporter, KratosCoreFastSuite)
{
    typedef std::size_t IndexType;
    DataCommunicator& rComm=ParallelEnvironment::GetDefaultDataCommunicator();

    IndexType local_size = 4;
    DistributedNumbering<IndexType> numbering(rComm,local_size);

    IndexType total_size =numbering.Size();
    DistributedSystemVector<double,IndexType> x(numbering);

    //test exporting
    std::vector<IndexType> indices_to_export;
    std::vector<double> data_to_export;
    if(rComm.Rank() == 0)
    {
        indices_to_export.push_back(0);
        data_to_export.push_back(5.0);

        indices_to_export.push_back(total_size/2);
        data_to_export.push_back(9.0);

        indices_to_export.push_back(total_size-1);
        data_to_export.push_back(15.0);
    }
    DistributedVectorExporter<IndexType> exporter(rComm,indices_to_export,x.GetNumbering());
    exporter.Apply(x,data_to_export);
    if(x.GetNumbering().IsLocal(0))
        KRATOS_CHECK_NEAR(x[x.GetNumbering().LocalId(0)], 5.0, 1e-14);
    if(x.GetNumbering().IsLocal(total_size/2))
        KRATOS_CHECK_NEAR(x[x.GetNumbering().LocalId(total_size/2)], 9.0, 1e-14);
    if(x.GetNumbering().IsLocal(total_size-1))
        KRATOS_CHECK_NEAR(x[x.GetNumbering().LocalId(total_size-1)], 15.0, 1e-14);
}

} // namespace Testing
} // namespace Kratos
