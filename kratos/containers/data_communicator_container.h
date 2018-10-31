//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Jordi Cotela
//

#ifndef KRATOS_DATA_COMMUNICATOR_CONTAINER_H_INCLUDED
#define KRATOS_DATA_COMMUNICATOR_CONTAINER_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <unordered_map>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/data_communicator.h"

namespace Kratos
{
///@addtogroup Kratos Core
///@{

///@name Kratos Classes
///@{

/// A container class for DataCommunicator pointers.
/** This class is intended to be used to hold the prototypes for the registration of DataCommunicators.
 *  @see DataCommunicator, MPIDataCommunicator.
 */
class DataCommunicatorContainer
{
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DataCommunicatorContainer
    KRATOS_CLASS_POINTER_DEFINITION(DataCommunicatorContainer);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DataCommunicatorContainer();

    /// Destructor.
    ~DataCommunicatorContainer();

    ///@}
    ///@name Operations
    ///@{

    void Register(const std::string Name, const DataCommunicator& rPrototype);

    ///@}
    ///@name Access
    ///@{

    DataCommunicator& Get(const std::string Name) const;

    ///@}
    ///@name Inquiry
    ///@{

    bool Has(const std::string Name) const;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const;

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const;

    /// Print object's data.
    void PrintData(std::ostream &rOStream) const;

    ///@}

  private:

    ///@name Member Variables
    ///@{

    std::unordered_map<std::string, DataCommunicator*> mpDataCommunicators;

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    DataCommunicatorContainer &operator=(DataCommunicatorContainer const &rOther) = delete;

    /// Copy constructor.
    DataCommunicatorContainer(DataCommunicatorContainer const &rOther) = delete;

    ///@}

}; // Class DataCommunicatorContainer

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream &operator>>(std::istream &rIStream,
                                DataCommunicatorContainer &rThis);

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
                                const DataCommunicatorContainer &rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_DATA_COMMUNICATOR_CONTAINER_H_INCLUDED  defined
