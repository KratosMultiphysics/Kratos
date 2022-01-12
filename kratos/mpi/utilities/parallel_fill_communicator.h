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
//

#if !defined(KRATOS_PARALLEL_FILL_COMMUNICATOR_H_INCLUDED )
#define  KRATOS_PARALLEL_FILL_COMMUNICATOR_H_INCLUDED

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

class ModelPart;

/// This function recomputes the communication plan for MPI
/** The objective of this class is to read the mesh owned by each node in a distributed context
 * and to fill the communication plan (coloring) so to allow the communication to be performed correctly
 * It fills the Ghost and Local lists and performs the coloring, then it updates the MPI communicator
 */
class KRATOS_API(KRATOS_MPI_CORE) ParallelFillCommunicator : public FillCommunicator
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ParallelFillCommunicator
    KRATOS_CLASS_POINTER_DEFINITION(ParallelFillCommunicator);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    KRATOS_DEPRECATED_MESSAGE("This constructor is deprecated, please use the one that accepts a DataCommunicator")
    ParallelFillCommunicator(ModelPart& rModelPart);

    ParallelFillCommunicator(
        ModelPart& rModelPart,
        const DataCommunicator& rDataComm);

    /// Destructor.
    virtual ~ParallelFillCommunicator() = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void Execute() override;

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

    void ComputeCommunicationPlan(ModelPart& rModelPart);

    /// Initialize the communicator's ghost, local and interface meshes for all communication pairs (colors).
    void InitializeParallelCommunicationMeshes(ModelPart& rModelPart,
                                               const std::vector<int>& rColors,
                                               int MyRank);

    /// Generate the ghost, local and interface meshes for processes of a communication pair (color).
    void GenerateMeshes(int NeighbourPID, int MyPID, unsigned Color, ModelPart& rModelPart);

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

    bool mPartitionIndexCheckPerformed = false;

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

#endif // KRATOS_PARALLEL_FILL_COMMUNICATOR_H_INCLUDED  defined


