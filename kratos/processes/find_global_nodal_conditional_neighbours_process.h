//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

#if !defined(KRATOS_FIND_GLOBAL_NODAL_CONDITION_NEIGHBOURS_PROCESS_H_INCLUDED)
#define KRATOS_FIND_GLOBAL_NODAL_CONDITION_NEIGHBOURS_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/global_pointer_variables.h"
#include "processes/find_global_nodal_neighbours_for_entities_process.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
 */
class FindGlobalNodalConditionalNeighboursProcess
    : public FindNodalNeighboursForEntitiesProcess<ModelPart::ConditionsContainerType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FindGlobalNodalConditionalNeighboursProcess
    KRATOS_CLASS_POINTER_DEFINITION(FindGlobalNodalConditionalNeighboursProcess);

    using BaseType = FindNodalNeighboursForEntitiesProcess<ModelPart::ConditionsContainerType>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    FindGlobalNodalConditionalNeighboursProcess(
        ModelPart& rModelPart)
        : BaseType(rModelPart, NEIGHBOUR_CONDITION_NODES)
    {
    }

    /// Destructor.
    ~FindGlobalNodalConditionalNeighboursProcess() override = default;

    /// Assignment operator.
    FindGlobalNodalConditionalNeighboursProcess& operator=(
        FindGlobalNodalConditionalNeighboursProcess const& rOther) = delete;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "FindGlobalNodalConditionalNeighboursProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "FindGlobalNodalConditionalNeighboursProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}

}; // Class FindGlobalNodalConditionalNeighboursProcess

///@}

///@name Input and output
///@{

/// input stream function
inline std::istream& operator>>(
    std::istream& rIStream,
    FindGlobalNodalConditionalNeighboursProcess& rThis);

/// output stream function
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const FindGlobalNodalConditionalNeighboursProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

} // namespace Kratos.

#endif // KRATOS_FIND_GLOBAL_NODAL_CONDITION_NEIGHBOURS_PROCESS_H_INCLUDED defined
