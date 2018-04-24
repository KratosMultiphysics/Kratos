//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

// System includes

// External includes

// Project includes
#include "interface_preprocessor.h"

#ifdef KRATOS_USING_MPI // mpi-parallel compilation
#include "custom_utilities/parallel_fill_communicator.h" // in trilinos-application
#endif

namespace Kratos
{
    /***********************************************************************************/
    /* PUBLIC Methods */
    /***********************************************************************************/
    InterfacePreprocessor::InterfacePreprocessor(ModelPart& rModelPartDestination,
                                                 ModelPartPointerType pInterfaceModelPart)
        : mrModelPartDestination(rModelPartDestination),
          mpInterfaceModelPart(pInterfaceModelPart)
    {

    }

    void InterfacePreprocessor::GenerateInterfaceModelPart(Parameters InterfaceParameters)
    {
        CheckAndValidateParameters(InterfaceParameters);

        // here the InterfaceModelPart is resetted to start from a clear state
        mpInterfaceModelPart = Kratos::make_shared<ModelPart>("Mapper-Interface");

        // Adding the Nodes
        mpInterfaceModelPart->Nodes() = mrModelPartDestination.Nodes();


// #ifdef KRATOS_USING_MPI // mpi-parallel compilation
//         // if (MapperUtilities::TotalProcesses() > 1)
//         // {
//             // TODO check if this is actually necessary
//             // Set the MPICommunicator
//             std::cout << "Doing the ParallelFillCommunicator stuff" << std::endl;
//             ParallelFillCommunicator parallel_fill_communicator(*mpInterfaceModelPart);
//             parallel_fill_communicator.Execute();
//         // }
// #endif

    }
    /***********************************************************************************/
    /* PROTECTED Methods */
    /***********************************************************************************/


    /***********************************************************************************/
    /* PRIVATE Methods */
    /***********************************************************************************/
    void InterfacePreprocessor::CheckAndValidateParameters(Parameters InterfaceParameters)
    {
        InterfaceParameters.RecursivelyValidateAndAssignDefaults(mDefaultParameters);

        KRATOS_ERROR_IF(InterfaceParameters["mapper_condition_name"].GetString() == "")
            << "Condition name for Interface-ModelPart not specified" << std::endl;
    }

    void InterfacePreprocessor::CreateMapperConditions()
    {

    }


}  // namespace Kratos.
