//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//

// Project includes
#include "testing/testing.h"
#include "includes/model_part.h"
#include "includes/mpi_communicator.h"
#include "mapping_application_variables.h"
#include "custom_utilities/mapper_utilities.h"

namespace Kratos {
namespace Testing {

void CreateNodesForMapping(ModelPart& rModelPart, const int NumNodes)
{
    const int rank = rModelPart.GetCommunicator().MyPID();
    const int size = rModelPart.GetCommunicator().TotalProcesses();

    const int start_id = NumNodes * rank + 1;

    // creating nodes with random coordinates
    for (int i=0; i< NumNodes; ++i)
        rModelPart.CreateNewNode(i+start_id, i*0.1*rank*size+0.134,
                                             i*0.2+rank*3.48*size,
                                             i*0.3*rank*6.13*size);
}

// Function to check if the INTERFACE_EQUATION_IDs (EquationIds of the MappingMatrix) are being assigned correctly
// Note that this test also exists for MPI
KRATOS_TEST_CASE_IN_SUITE(MapperUtilities_AssignInterfaceEquationIds_InMPI, KratosMappingApplicationGeneralTestSuite)
{
    const int num_nodes = 11;
    ModelPart model_part("ForTest");

// In MPI we replace the Comunicator with the MPICommunicator
#ifdef KRATOS_USING_MPI // mpi-parallel compilation
    int mpi_initialized;
    MPI_Initialized(&mpi_initialized);
    if (mpi_initialized)   // parallel execution, i.e. mpi imported in python
    {
        int comm_size;
        MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
        if (comm_size > 1)
        {
            // Note that this MPICommunicator
            VariablesList* var_list = &(model_part.GetNodalSolutionStepVariablesList()) ;
            model_part.SetCommunicator(Communicator::Pointer(new MPICommunicator(var_list)));
        }
    }
#endif

    CreateNodesForMapping(model_part, num_nodes);

    MapperUtilities::AssignInterfaceEquationIds(model_part.GetCommunicator());

    const int rank = model_part.GetCommunicator().MyPID();

    int idx = num_nodes * rank; // this simulates the ScanSum

    for (const auto& r_node : model_part/*.GetCommunicator().LocalMesh()*/.Nodes())
    {
        KRATOS_INFO("Rank") << rank << " , idx = " << idx << " , EdId: " << r_node.GetValue(INTERFACE_EQUATION_ID) << std::endl;
        KRATOS_CHECK_EQUAL(idx, r_node.GetValue(INTERFACE_EQUATION_ID));
        idx += 1;
    }
}

}  // namespace Testing
}  // namespace Kratos