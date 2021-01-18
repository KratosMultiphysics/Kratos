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
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/global_pointer_variables.h"
#include "utilities/communication_coloring_utilities.h"
#include "utilities/pointer_communicator.h"
#include "includes/mpi_serializer.h"
#include "utilities/compute_neighbour_list_functor.h"
#include "utilities/global_pointer_utilities.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
typedef  ModelPart::NodesContainerType NodesContainerType;
typedef  ModelPart::ElementsContainerType ElementsContainerType;


///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class KRATOS_API(KRATOS_CORE) FindGlobalNodalElementalNeighboursProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FindGlobalNodalElementalNeighboursProcess
    KRATOS_CLASS_POINTER_DEFINITION(FindGlobalNodalElementalNeighboursProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /// the better the guess for the quantities above the less memory occupied and the fastest the algorithm
    FindGlobalNodalElementalNeighboursProcess(const DataCommunicator& rComm, 
                                              ModelPart& model_part)
        : mrComm(rComm),mr_model_part(model_part)
    {
    }

    /// Destructor.
    ~FindGlobalNodalElementalNeighboursProcess() override
    {
    }


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

    std::unordered_map<int, std::vector<int> > GetNeighbourIds(
                ModelPart::NodesContainerType& rNodes
                );

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
    std::string Info() const override
    {
        return "FindGlobalNodalElementalNeighboursProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "FindGlobalNodalElementalNeighboursProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }


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
    const DataCommunicator& mrComm;
    ModelPart& mr_model_part;
    unsigned int mavg_nodes;


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
    FindGlobalNodalElementalNeighboursProcess& operator=(FindGlobalNodalElementalNeighboursProcess const& rOther);

    /// Copy constructor.
    //FindGlobalNodalElementalNeighboursProcess(FindGlobalNodalElementalNeighboursProcess const& rOther);


    ///@}

}; // Class FindGlobalNodalElementalNeighboursProcess

///@}

///@name Type Definitions
///@{


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


