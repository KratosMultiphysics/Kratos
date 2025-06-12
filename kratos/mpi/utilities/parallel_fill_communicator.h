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
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/fill_communicator.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/** 
 * @class ParallelFillCommunicator
 * @ingroup KratosMPI
 * @brief This function recomputes the communication plan for MPI
 * @details The objective of this class is to read the mesh owned by each node in a distributed context
 * and to fill the communication plan (coloring) so to allow the communication to be performed correctly
 * It fills the Ghost and Local lists and performs the coloring, then it updates the MPI communicator
 * @author Riccardo Rossi
 */
class KRATOS_API(KRATOS_MPI_CORE) ParallelFillCommunicator 
    : public FillCommunicator
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ParallelFillCommunicator
    KRATOS_CLASS_POINTER_DEFINITION(ParallelFillCommunicator);

    ///@}
    ///@name Life Cycle
    ///@{

    /** 
     * @brief Constructor (deprecated)
     * @param rModelPart The model part to recompute the communication plan for MPI
     */
    KRATOS_DEPRECATED_MESSAGE("This constructor is deprecated, please use the one that accepts a DataCommunicator")
    ParallelFillCommunicator(ModelPart& rModelPart);

    /** 
     * @brief Constructor.
     * @param rModelPart The model part to recompute the communication plan for MPI
     * @param rDataCommunicator The communicator to recompute the communication plan for MPI
     */
    ParallelFillCommunicator(
        ModelPart& rModelPart,
        const DataCommunicator& rDataCommunicator
        );

    /// Destructor.
    virtual ~ParallelFillCommunicator() = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Execute the communicator fill
     * @details This method is intended to perform the communicator filling
     */
    void Execute() override;

    /**
     * @brief Function to print mesh information of the provided model part
     * @details This function is intended to check and print some mesh information
     * @param rModelPart Reference to the model part to be checked
     */
    void PrintModelPartDebugInfo(const ModelPart& rModelPart) override;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

    ///@}
    ///@name Friends
    ///@{

    ///@}
protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief This function computes the communication plan
     * @param rModelPart The model part to compute the communication plan
    */
    void ComputeCommunicationPlan(ModelPart& rModelPart);

    /// Initialize the communicator's ghost, local and interface meshes for all communication pairs (colors).
    void InitializeParallelCommunicationMeshes(
        ModelPart& rModelPart,
        const std::vector<int>& rColors,
        const int MyRank
        );

    /// Generate the ghost, local and interface meshes for processes of a communication pair (color).
    void GenerateMeshes(
        const int NeighbourPID, 
        const int MyPID, 
        const unsigned int Color, 
        ModelPart& rModelPart
        );

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    bool mPartitionIndexCheckPerformed = false; /// If the partition index check is performed

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.

    ParallelFillCommunicator& operator=(ParallelFillCommunicator const& rOther) = delete;

    /// Copy constructor.

    ParallelFillCommunicator(ParallelFillCommunicator const& rOther) = delete;

    ///@}

}; // Class ParallelFillCommunicator

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream & operator >>(std::istream& rIStream,
                                  ParallelFillCommunicator& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream & operator <<(std::ostream& rOStream,
                                  const ParallelFillCommunicator& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

} // namespace Kratos.
