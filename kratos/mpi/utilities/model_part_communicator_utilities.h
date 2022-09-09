//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Jordi Cotela
//

#ifndef KRATOS_MODEL_PART_COMMUNICATOR_UTILITIES_H_INCLUDED
#define	KRATOS_MODEL_PART_COMMUNICATOR_UTILITIES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "mpi/includes/mpi_communicator.h"

namespace Kratos
{

///@addtogroup MPICore
///@{

///@name Kratos Classes
///@{

/// Utilitiy class for ModelPart::Comunicator management in an MPI context.
class ModelPartCommunicatorUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ModelPartCommunicatorUtilities
    KRATOS_CLASS_POINTER_DEFINITION(ModelPartCommunicatorUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Deleted default constructor.
    ModelPartCommunicatorUtilities() = delete;

    /// Deleted copy constructor
    ModelPartCommunicatorUtilities(const ModelPartCommunicatorUtilities& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /// Create and assign an MPICommunicator for a ModelPart instance.
    /** Note that this does not initialize the Communicator, since the ModelPart
     *  may not yet contain nodes (or anything else) when this is called.
     *  @param rThisModelPart The ModelPart that will get an MPICommunicator.
     *  @param rDataCommunicator The DataCommunicator that will be used for the MPICommunicator.
     */
    static inline void SetMPICommunicator(ModelPart& rThisModelPart, const DataCommunicator& rDataCommunicator)
    {
        KRATOS_ERROR_IF_NOT(rDataCommunicator.IsDistributed()) << "Only distributed DataCommunicators can be used!" << std::endl;
        VariablesList * p_variables_list = &rThisModelPart.GetNodalSolutionStepVariablesList();
        rThisModelPart.SetCommunicator(Kratos::make_shared<MPICommunicator>(p_variables_list, rDataCommunicator));
    }

    /// Create and assign an MPICommunicator for a ModelPart instance and its SubModelParts.
    /** Note that this does not initialize the Communicator, since the ModelPart
     *  may not yet contain nodes (or anything else) when this is called.
     *  @param rThisModelPart The ModelPart that will get an MPICommunicator.
     *  @param rDataCommunicator The DataCommunicator that will be used for the MPICommunicator.
     */
    static inline void SetMPICommunicatorRecursively(ModelPart& rThisModelPart, const DataCommunicator& rDataCommunicator)
    {
        SetMPICommunicator(rThisModelPart, rDataCommunicator);
        for (auto& r_smp : rThisModelPart.SubModelParts()) {
            SetMPICommunicatorRecursively(r_smp, rDataCommunicator);
        }
    }

    KRATOS_DEPRECATED_MESSAGE("This function is deprecated, please use the one that accepts a DataCommunicator") static inline void SetMPICommunicator(ModelPart& rThisModelPart)
    {
        VariablesList * p_variables_list = &rThisModelPart.GetNodalSolutionStepVariablesList();
        rThisModelPart.SetCommunicator(Kratos::make_shared<MPICommunicator>(p_variables_list, DataCommunicator::GetDefault()));
    }

    ///@}

}; // Class ModelPartCommunicatorUtilities

///@}

///@}

}  // namespace Kratos.

#endif	/* KRATOS_MODEL_PART_COMMUNICATOR_UTILITIES_H_INCLUDED */

