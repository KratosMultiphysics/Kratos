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
class KRATOS_API(KRATOS_CORE) ParallelEnvironment
{
  public:
    ///@name Type Definitions
    ///@{

    constexpr static bool MakeDefault = true;
    constexpr static bool DoNotMakeDefault = false;

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Operations
    ///@{

    static void RegisterDataCommunicator(
        const std::string& rName,
        const DataCommunicator& rPrototype,
        const bool Default = DoNotMakeDefault);

    static DataCommunicator& GetDataCommunicator(const std::string& rName);

    static DataCommunicator& GetDefaultDataCommunicator();

    static void SetDefaultDataCommunicator(const std::string& rName);

    ///@}
    ///@name Inquiry
    ///@{

    static bool HasDataCommunicator(const std::string& rName);

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    static std::string Info();

    /// Print information about this object.
    static void PrintInfo(std::ostream &rOStream);

    /// Print object's data.
    static void PrintData(std::ostream &rOStream);

    ///@}

  private:

    ///@name Private Life Cycle
    ///@{

    /// Internal constructor.
    ParallelEnvironment();

    ///@}
    ///@name Private Operations
    ///@{

    void RegisterDataCommunicatorDetail(
        const std::string& Name,
        const DataCommunicator& rPrototype,
        const bool Default = DoNotMakeDefault);

    DataCommunicator& GetDataCommunicatorDetail(const std::string& rName) const;

    DataCommunicator& GetDefaultDataCommunicatorDetail() const;

    void SetDefaultDataCommunicatorDetail(const std::string& rName);

    ///@}
    ///@name Private Inquiry
    ///@{

    bool HasDataCommunicatorDetail(const std::string& rName) const;

    ///@}
    ///@name Private Access
    ///@{

    static ParallelEnvironment& GetInstance();

    ///@}
    ///@name Private Input and output
    ///@{

    /// Turn back information as a string.
    std::string InfoDetail() const;

    /// Print information about this object.
    void PrintInfoDetail(std::ostream &rOStream) const;

    /// Print object's data.
    void PrintDataDetail(std::ostream &rOStream) const;

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

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_PARALLEL_ENVIRONMENT_H_INCLUDED  defined
