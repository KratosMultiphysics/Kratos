//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//


#if !defined(KRATOS_DISTRIBUTED_MODEL_PART_INITIALIZER_H_INCLUDED )
#define  KRATOS_DISTRIBUTED_MODEL_PART_INITIALIZER_H_INCLUDED


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"


namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/// Initialize a distributed ModelPart from a serial one.
/** This class initializes a distributed ModelPart from a serial one.
 * It creates the ModelPart hierarchy that exists on the source rank
 * also on the other ranks.
 * Furthermore it initializes the (MPI-)Communicators.
 * Note that all the entities are still only on the source rank,
 * no partitioning is done!
*/
class DistributedModelPartInitializer
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DistributedModelPartInitializer
    KRATOS_CLASS_POINTER_DEFINITION(DistributedModelPartInitializer);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    DistributedModelPartInitializer(
        ModelPart& rModelPart,
        const DataCommunicator& rDataComm,
        int SourceRank)
        : mrModelPart(rModelPart),
          mrDataComm(rDataComm),
          mSourceRank(SourceRank)
        { }

    /// Destructor.
    virtual ~DistributedModelPartInitializer() = default;

    ///@}
    ///@name Operations
    ///@{

    void CopySubModelPartStructure();

    void Execute();

    ///@}

private:
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;
    const DataCommunicator& mrDataComm;
    int mSourceRank;

    ///@}

}; // Class DistributedModelPartInitializer

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_DISTRIBUTED_MODEL_PART_INITIALIZER_H_INCLUDED defined
