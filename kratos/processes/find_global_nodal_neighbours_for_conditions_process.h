//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#if !defined(KRATOS_FIND_GLOBAL_NODAL_NEIGHBOURS_FOR_CONDITIONS_PROCESS_H_INCLUDED)
#define KRATOS_FIND_GLOBAL_NODAL_NEIGHBOURS_FOR_CONDITIONS_PROCESS_H_INCLUDED

// System includes
#include <iostream>
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/global_pointer_variables.h"
#include "includes/model_part.h"
#include "processes/find_global_nodal_neighbours_for_entities_process.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
 */
class FindGlobalNodalNeighboursForConditionsProcess
    : public FindNodalNeighboursForEntitiesProcess<ModelPart::ConditionsContainerType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FindGlobalNodalNeighboursForConditionsProcess
    KRATOS_CLASS_POINTER_DEFINITION(FindGlobalNodalNeighboursForConditionsProcess);

    using BaseType =
        FindNodalNeighboursForEntitiesProcess<ModelPart::ConditionsContainerType>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// constructor.
    FindGlobalNodalNeighboursForConditionsProcess(
        const DataCommunicator& rDataCommunicator,
        ModelPart& rModelPart)
        : BaseType(rDataCommunicator, rModelPart, NEIGHBOUR_CONDITION_NODES)
    {
    }

    /// Destructor.
    ~FindGlobalNodalNeighboursForConditionsProcess() override = default;

    /// Assignment operator.
    FindGlobalNodalNeighboursForConditionsProcess& operator=(
        FindGlobalNodalNeighboursForConditionsProcess const& rOther) = delete;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "FindGlobalNodalNeighboursForConditionsProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "FindGlobalNodalNeighboursForConditionsProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}

}; // Class FindGlobalNodalNeighboursForConditionsProcess

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator>>(
    std::istream& rIStream,
    FindGlobalNodalNeighboursForConditionsProcess& rThis);

/// output stream function
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const FindGlobalNodalNeighboursForConditionsProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

} // namespace Kratos.

#endif // KRATOS_FIND_GLOBAL_NODAL_NEIGHBOURS_FOR_CONDITIONS_PROCESS_H_INCLUDED defined
