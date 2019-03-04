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

    /// Create and assign an MPI Communicator for a ModelPart instance.
    /** Note that this does not initialize the communicator, since the ModelPart
     *  may not yet contain nodes (or anything else) when this is called.
     *  @param rThisModelPart The model part that will get an MPI communicator.
     */
    static inline void SetMPICommunicator(ModelPart& rThisModelPart)
    {
        VariablesList * p_variables_list = &rThisModelPart.GetNodalSolutionStepVariablesList();
        rThisModelPart.SetCommunicator(Kratos::make_shared<MPICommunicator>(p_variables_list));
    }

    ///@}

}; // Class ModelPartCommunicatorUtilities

///@}

///@}

}  // namespace Kratos.

#endif	/* KRATOS_MODEL_PART_COMMUNICATOR_UTILITIES_H_INCLUDED */

