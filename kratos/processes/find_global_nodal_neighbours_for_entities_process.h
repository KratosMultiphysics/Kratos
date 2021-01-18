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

#if !defined(KRATOS_FIND_NODAL_NEIGHBOURS_FOR_ENTITIES_PROCESS_H_INCLUDED)
#define KRATOS_FIND_NODAL_NEIGHBOURS_FOR_ENTITIES_PROCESS_H_INCLUDED

// System includes
#include <iostream>
#include <string>
#include <unordered_map>

// External includes

// Project includes
#include "containers/global_pointers_vector.h"
#include "includes/data_communicator.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
 */
template <class TContainerType>
class KRATOS_API(KRATOS_CORE) FindNodalNeighboursForEntitiesProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FindNodalNeighboursForEntitiesProcess
    KRATOS_CLASS_POINTER_DEFINITION(FindNodalNeighboursForEntitiesProcess);

    using NodeType = ModelPart::NodeType;
    using NodesContainerType = ModelPart::NodesContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    FindNodalNeighboursForEntitiesProcess(
        const DataCommunicator& rDataCommunicator,
        ModelPart& rModelPart,
        const Variable<GlobalPointersVector<NodeType>>& rOutputVariable)
        : mrModelPart(rModelPart),
          mrDataCommunicator(rDataCommunicator),
          mrOutputVariable(rOutputVariable)
    {
    }

    /// Destructor.
    ~FindNodalNeighboursForEntitiesProcess() override = default;

    /// Assignment operator.
    FindNodalNeighboursForEntitiesProcess& operator=(
        FindNodalNeighboursForEntitiesProcess const& rOther) = delete;

    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    void Execute() override;

    void ClearNeighbours();

    std::unordered_map<int, std::vector<int>> GetNeighbourIds(
        NodesContainerType& rNodes) const
    {
        return FindNodalNeighboursForEntitiesProcess<TContainerType>::GetNodalNeighbourIdsMap(
            rNodes, this->mrDataCommunicator, this->mrOutputVariable);
    }

    static std::unordered_map<int, std::vector<int>> GetNodalNeighbourIdsMap(
        NodesContainerType& rNodes,
        const DataCommunicator& rDataCommunicator,
        const Variable<GlobalPointersVector<NodeType>>& rNodalGlobalPointerVariable);

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "FindNodalNeighboursForEntitiesProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "FindNodalNeighboursForEntitiesProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;
    const DataCommunicator& mrDataCommunicator;
    const Variable<GlobalPointersVector<NodeType>>& mrOutputVariable;

    ///@}
    ///@name Private Operators
    ///@{

    TContainerType& GetContainer();

    void AddHangingNodeIds(
        std::unordered_map<int, std::unordered_map<int, std::vector<int>>>& rNeighbourIds) const;

    void AddUnique(
        std::vector<int>& rContainer,
        const int Item) const;

    template <class TDataType>
    void AddUniqueGlobalPointer(
        GlobalPointersVector<TDataType>& rGPVector,
        const GlobalPointer<TDataType>& rCandidate) const
    {
        bool found = false;
        for (const auto& gp : rGPVector.GetContainer()) {
            if (&(*gp) == &(*rCandidate) && gp.GetRank() == rCandidate.GetRank()) {
                found = true;
                break;
            }
        }
        if (!found) {
            rGPVector.push_back(rCandidate);
        }
    }

    ///@}
}; // Class FindNodalNeighboursForEntitiesProcess

///@}
///@name Input and output
///@{

/// input stream function
template <class TContainerType>
inline std::istream& operator>>(
    std::istream& rIStream,
    FindNodalNeighboursForEntitiesProcess<TContainerType>& rThis);

/// output stream function
template <class TContainerType>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const FindNodalNeighboursForEntitiesProcess<TContainerType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

} // namespace Kratos.

#endif // KRATOS_FIND_NODAL_NEIGHBOURS_FOR_ENTITIES_PROCESS_H_INCLUDED defined
