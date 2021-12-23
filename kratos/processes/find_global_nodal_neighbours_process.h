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
//                   Suneth Warnakulasuriya
//

#if !defined(KRATOS_FIND_GLOBAL_NODAL_NEIGHBOURS_PROCESS_H_INCLUDED)
#define KRATOS_FIND_GLOBAL_NODAL_NEIGHBOURS_PROCESS_H_INCLUDED

// System includes
#include <iostream>
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/global_pointer_variables.h"
#include "includes/model_part.h"
#include "processes/find_global_nodal_neighbours_for_entities_process.h"
#include "utilities/pointer_communicator.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
 */
class FindGlobalNodalNeighboursProcess
    : public FindNodalNeighboursForEntitiesProcess<ModelPart::ElementsContainerType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FindGlobalNodalNeighboursProcess
    KRATOS_CLASS_POINTER_DEFINITION(FindGlobalNodalNeighboursProcess);

    using BaseType = FindNodalNeighboursForEntitiesProcess<ModelPart::ElementsContainerType>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    FindGlobalNodalNeighboursProcess(
        ModelPart& rModelPart)
        : BaseType(rModelPart, NEIGHBOUR_NODES)
    {
    }

    /// avg_elems ------ expected number of neighbour elements per node.,
    /// avg_nodes ------ expected number of neighbour Nodes
    /// the better the guess for the quantities above the less memory occupied and the fastest the algorithm
    KRATOS_DEPRECATED_MESSAGE(
        "Use of avg_nodes is deprecated. Please use constructor without it.")
    FindGlobalNodalNeighboursProcess(
        const DataCommunicator& rDataCommunicator,
        ModelPart& rModelPart,
        unsigned int avg_nodes)
        : FindGlobalNodalNeighboursProcess(rModelPart)
    {
    }

    KRATOS_DEPRECATED_MESSAGE("Use of DataCommunicator is deprecated. Please use constructor without it.")
    FindGlobalNodalNeighboursProcess(
        const DataCommunicator& rDataCommunicator,
        ModelPart& rModelPart)
        : FindGlobalNodalNeighboursProcess(rModelPart)
    {
    }

    /// Destructor.
    ~FindGlobalNodalNeighboursProcess() override = default;

    /// Assignment operator.
    FindGlobalNodalNeighboursProcess& operator=(
        FindGlobalNodalNeighboursProcess const& rOther) = delete;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "FindGlobalNodalNeighboursProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "FindGlobalNodalNeighboursProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}

}; // Class FindGlobalNodalNeighboursProcess

///@}

///@name Input and output
///@{

/// input stream function
inline std::istream& operator>>(
    std::istream& rIStream,
    FindGlobalNodalNeighboursProcess& rThis);

/// output stream function
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const FindGlobalNodalNeighboursProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

} // namespace Kratos.

#endif // KRATOS_FIND_GLOBAL_NODAL_NEIGHBOURS_PROCESS_H_INCLUDED  defined
