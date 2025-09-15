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

#pragma once

// System includes
#include <string>
#include <iostream>
#include <unordered_map>

// External includes

// Project includes
#include "includes/define.h"
#include "containers/model.h"
#include "processes/process.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "includes/global_pointer_variables.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/

template<class TContainerType>
class KRATOS_API(KRATOS_CORE) FindGlobalNodalEntityNeighboursProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    using NodeType = ModelPart::NodeType;

    using EntityType = typename TContainerType::value_type;

    using GlobalEntityPointersVectorType = GlobalPointersVector<EntityType>;

    //contains id vs vector_of_neighbours
    using NeighbourMapType = std::unordered_map<int, GlobalPointersVector<EntityType>>;

    using NonLocalMapType =  std::unordered_map<int, NeighbourMapType>;

    using IdMapType = std::unordered_map<int, std::vector<int>>;

    /// Pointer definition of FindGlobalNodalEntityNeighboursProcess
    KRATOS_CLASS_POINTER_DEFINITION(FindGlobalNodalEntityNeighboursProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    FindGlobalNodalEntityNeighboursProcess(
        Model& rModel,
        Parameters Params);

    FindGlobalNodalEntityNeighboursProcess(
        ModelPart& rModelPart);

    FindGlobalNodalEntityNeighboursProcess(
        ModelPart& rModelPart,
        const Variable<GlobalEntityPointersVectorType>& rOutputVariable);

    /// Destructor.
    ~FindGlobalNodalEntityNeighboursProcess() override = default;

    ///@}
    ///@name Operators
    ///@{

    FindGlobalNodalEntityNeighboursProcess<EntityType>& operator=(FindGlobalNodalEntityNeighboursProcess<EntityType> const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    void Execute() override;

    void Clear() override;

    KRATOS_DEPRECATED_MESSAGE("This is legacy version (use Clear)") void ClearNeighbours() { Clear(); }

    IdMapType GetNeighbourIds(const ModelPart::NodesContainerType& rNodes);

    const Parameters GetDefaultParameters() const override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "FindGlobalNodalEntityNeighboursProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "FindGlobalNodalEntityNeighboursProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }


    ///@}
private:
    ///@name Member Variables
    ///@{

    Model& mrModel;

    std::string mModelPartName;

    const Variable<GlobalEntityPointersVectorType>& mrOutputVariable;

    ///@}
    ///@name Private Operations
    ///@{

    static TContainerType& GetContainer(ModelPart& rModelPart);

    static const Variable<GlobalEntityPointersVectorType>& GetDefaultOutputVariable();

    ///@}

}; // Class FindGlobalNodalEntityNeighboursProcess

///@}
///@name Input and output
///@{


/// input stream function
template<class TContainerType>
inline std::istream& operator >> (std::istream& rIStream,
                                  FindGlobalNodalEntityNeighboursProcess<TContainerType>& rThis);

/// output stream function
template<class TContainerType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const FindGlobalNodalEntityNeighboursProcess<TContainerType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.
