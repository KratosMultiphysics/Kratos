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
        mpInterfaceModelPart->AddNodalSolutionStepVariable(PARTITION_INDEX); // I think this is needed in the ParallelFillCommunicator
    }

    void InterfacePreprocessor::GenerateInterfaceModelPart(Parameters InterfaceParameters)
    {
        CheckAndValidateParameters(InterfaceParameters);

        mpInterfaceModelPart->GetMesh().Clear();

        // create dummy properties for the mapper conditions
        Properties::Pointer dummy_properties = Kratos::make_shared<Properties>();
        mpInterfaceModelPart->AddProperties(dummy_properties);

        // Adding the Nodes
        mpInterfaceModelPart->Nodes() = mrModelPartDestination.Nodes();

        CreateMapperConditions(InterfaceParameters);


#ifdef KRATOS_USING_MPI // mpi-parallel compilation
        if (mrModelPartDestination.GetCommunicator().TotalProcesses() > 1)
        {
            // TODO check if this is actually necessary
            // Set the MPICommunicator
            std::cout << "Doing the ParallelFillCommunicator stuff" << std::endl;

            ParallelFillCommunicator parallel_fill_communicator(*mpInterfaceModelPart);
            // parallel_fill_communicator.Execute(); // TODO does not work atm
        }
#endif

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

    void InterfacePreprocessor::CreateMapperConditions(Parameters InterfaceParameters)
    {
        if (InterfaceParameters["use_nodes"].GetBool())
            CreateMapperConditionsFromNodes(InterfaceParameters);
        else
            CreateMapperConditionsFromGeometries(InterfaceParameters);
    }

    void InterfacePreprocessor::CreateMapperConditionsFromNodes(Parameters InterfaceParameters)
    {
        const std::string condition_name = InterfaceParameters["mapper_condition_name"].GetString();

        const Condition& rReferenceCondition = KratosComponents<Condition>::Get(condition_name);

        Condition::NodesArrayType condition_node;

        const auto p_props = mpInterfaceModelPart->pGetProperties(0); // dummy properties needed for the Create-fct

        for (ModelPart::NodesContainerType::const_iterator it_node = mpInterfaceModelPart->NodesBegin();
             it_node != mpInterfaceModelPart->NodesEnd(); ++it_node)
        {
            // Creating the mapper conditions
            // Using the ID of the Node which also works in MPI
            condition_node.clear();
            condition_node.push_back(*(it_node).base()); // TODO ask if this way is the most efficient one
            Condition::Pointer p_new_cond = Condition::Pointer(rReferenceCondition.Create(it_node->Id(),
                                                                                          condition_node,
                                                                                          p_props));
            mpInterfaceModelPart->AddCondition(p_new_cond);
        }
    }

    void InterfacePreprocessor::CreateMapperConditionsFromGeometries(Parameters InterfaceParameters)
    {
        KRATOS_ERROR << "This function is not implemented yet" << std::endl;
    }


}  // namespace Kratos.
