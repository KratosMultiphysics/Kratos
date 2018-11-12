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

#ifndef KRATOS_PARALLEL_ENVIRONMENT_H_INCLUDED
#define KRATOS_PARALLEL_ENVIRONMENT_H_INCLUDED

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

/// Holder for general data related to MPI (or suitable serial equivalents for non-MPI runs).
/** This class manages a registry of DataCommunicators, which can be used to perform MPI communication.
 *  @see DataCommunicator, MPIDataCommunicator.
 */
class ParallelEnvironment
{
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ParallelEnvironment
    KRATOS_CLASS_POINTER_DEFINITION(ParallelEnvironment);

    constexpr static bool MakeDefault = true;
    constexpr static bool DoNotMakeDefault = false;

    ///@}
    ///@name Access
    ///@{

    static ParallelEnvironment& GetInstance();

    ///@}
    ///@name Operations
    ///@{

    void RegisterDataCommunicator(
        const std::string Name,
        const DataCommunicator& rPrototype,
        bool Default = DoNotMakeDefault);

    DataCommunicator& GetDataCommunicator(const std::string Name) const;

    DataCommunicator& GetDefaultDataCommunicator() const;

    void SetDefaultDataCommunicator(const std::string Name);

    ///@}
    ///@name Inquiry
    ///@{

    bool HasDataCommunicator(const std::string Name) const;

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

    ///@name Private Life Cycle
    ///@{

    /// Internal constructor.
    ParallelEnvironment();

    ///@}
    ///@name Member Variables
    ///@{

    std::unordered_map<std::string, DataCommunicator::UniquePointer> mDataCommunicators;

    std::unordered_map<std::string, DataCommunicator::UniquePointer>::iterator mDefaultCommunicator;

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    ParallelEnvironment &operator=(ParallelEnvironment const &rOther) = delete;

    /// Copy constructor.
    ParallelEnvironment(ParallelEnvironment const &rOther) = delete;

    ///@}

}; // Class ParallelEnvironment

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream &operator>>(std::istream &rIStream,
                                ParallelEnvironment &rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
                                const ParallelEnvironment &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_PARALLEL_ENVIRONMENT_H_INCLUDED  defined
