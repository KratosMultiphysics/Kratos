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

#if !defined(KRATOS_FIND_GLOBAL_NODAL_ELEMENTAL_NEIGHBOURS_PROCESS_H_INCLUDED )
#define  KRATOS_FIND_GLOBAL_NODAL_ELEMENTAL_NEIGHBOURS_PROCESS_H_INCLUDED


// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "processes/find_global_nodal_entity_neighbours_process.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class FindGlobalNodalElementalNeighboursProcess
    : public FindGlobalNodalEntityNeighboursProcess<ModelPart::ElementsContainerType>
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = FindGlobalNodalEntityNeighboursProcess<ModelPart::ElementsContainerType>;

    /// Pointer definition of FindGlobalNodalElementalNeighboursProcess
    KRATOS_CLASS_POINTER_DEFINITION(FindGlobalNodalElementalNeighboursProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /// the better the guess for the quantities above the less memory occupied and the fastest the algorithm
    FindGlobalNodalElementalNeighboursProcess(ModelPart& rModelPart)
        : BaseType(rModelPart, NEIGHBOUR_ELEMENTS)
    {
    }

    KRATOS_DEPRECATED_MESSAGE("Use of DataCommunicator is deprecated. Please use constructor without it.")
    FindGlobalNodalElementalNeighboursProcess(const DataCommunicator& rComm,
                                              ModelPart& rModelPart)
        : FindGlobalNodalElementalNeighboursProcess(rModelPart)
    {
    }

    /// Destructor.
    ~FindGlobalNodalElementalNeighboursProcess() override = default;

    ///}

}; // Class FindGlobalNodalElementalNeighboursProcess

///@}

///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  FindGlobalNodalElementalNeighboursProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const FindGlobalNodalElementalNeighboursProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_FIND_GLOBAL_NODAL_ELEMENTAL_NEIGHBOURS_PROCESS_H_INCLUDED  defined


