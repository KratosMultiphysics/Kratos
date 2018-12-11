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

#if !defined(KRATOS_MPI_COLORING_UTILITIES_H_INCLUDED )
#define  KRATOS_MPI_COLORING_UTILITIES_H_INCLUDED

// System includes
#include <iostream>
#include <sstream>
#include <map>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"


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
class MPIColoringUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MPIColoringUtilities
    KRATOS_CLASS_POINTER_DEFINITION(MPIColoringUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MPIColoringUtilities() {}

    /// Destructor.
    virtual ~MPIColoringUtilities() {}

    ///@}
    ///@name Operators
    ///@{
    static std::vector<int> ComputeRecvList(
        const std::vector<int>& rLocalDestinationIds,
        MPIDataCommunicator& rComm
    );

    static std::vector<int> ComputeCommunicationScheduling(
        const std::vector<int>& rLocalDestinationIds,
        MPIDataCommunicator& rComm
    );




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
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "MPIColoringUtilities" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MPIColoringUtilities";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}


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
    static bool HasEdge(std::map<int, std::map<int, int> >& graph, int i, int j);


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

    friend class Serializer;

    void save(Serializer& rSerializer) const
    {
    }

    void load(Serializer& rSerializer)
    {
    }


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

    ///@}

}; // Class MPIColoringUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template< class TDataType >
inline std::istream& operator >> (std::istream& rIStream,
                                  MPIColoringUtilities& rThis)
{
    return rIStream;
}

/// output stream function
template< class TDataType >
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MPIColoringUtilities& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MPI_COLORING_UTILITIES_H_INCLUDED  defined 


