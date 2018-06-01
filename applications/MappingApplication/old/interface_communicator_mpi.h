//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

#if !defined(KRATOS_INTERFACE_COMMUNICATOR_MPI_H)
#define  KRATOS_INTERFACE_COMMUNICATOR_MPI_H

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "interface_communicator.h"


namespace Kratos
{
///@addtogroup ApplicationNameApplication
///@{

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

/// Short class definition.
/** Detail class definition.
*/
class InterfaceCommunicatorMPI : public InterfaceCommunicator
{
    public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of InterfaceCommunicatorMPI
    KRATOS_CLASS_POINTER_DEFINITION(InterfaceCommunicatorMPI);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    InterfaceCommunicatorMPI(ModelPart& rModelPartOrigin,
                             ModelPart::Pointer mpInterfaceModelPart);

    /// Destructor.
    virtual ~InterfaceCommunicatorMPI() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


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
    virtual std::string Info() const override
    {
        return "InterfaceCommunicatorMPI";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override {}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override {}


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
    // InterfaceCommunicatorMPI& operator=(InterfaceCommunicatorMPI const& rOther) {}

    /// Copy constructor.
    // InterfaceCommunicatorMPI(InterfaceCommunicatorMPI const& rOther) {}


    ///@}

}; // Class InterfaceCommunicatorMPI

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_INTERFACE_COMMUNICATOR_MPI_H  defined
