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

#ifndef KRATOS_MPI_DATA_COMMUNICATOR_H_INCLUDED
#define KRATOS_MPI_DATA_COMMUNICATOR_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <mpi.h>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/data_communicator.h"

namespace Kratos
{
///@addtogroup Kratos MPI Core
///@{

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
  */
class MPIDataCommunicator: public DataCommunicator
{
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MPIDataCommunicator
    KRATOS_CLASS_POINTER_DEFINITION(MPIDataCommunicator);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor accepting an MPI_Comm object.
    MPIDataCommunicator(MPI_Comm TheMPIComm)
    : DataCommunicator()
    , mComm(TheMPIComm)
    {}

    /// Destructor.
    ~MPIDataCommunicator() override {}

    ///@}
    ///@name Operations
    ///@{

    void Barrier() const override
    {
        MPI_Barrier(mComm);
    }

    ///@}
    ///@name Access
    ///@{

    MPI_Comm GetMPICommunicator() const override
    {
        return mComm;
    }

    ///@}
    ///@name Inquiry
    ///@{

    int Rank() const override
    {
        int rank;
        MPI_Comm_rank(mComm, &rank);
        return rank;
    }

    int Size() const override
    {
        int size;
        MPI_Comm_size(mComm, &size);
        return size;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        PrintInfo(buffer);
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "MPIDataCommunicator";
    }

    /// Print object's data.
    void PrintData(std::ostream &rOStream) const override {}

    ///@}

  protected:
    ///@name Protected LifeCycle
    ///@{

    ///@}

  private:
    ///@name Member Variables
    ///@{

    MPI_Comm mComm;

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    MPIDataCommunicator &operator=(MPIDataCommunicator const &rOther) = delete;

    /// Copy constructor.
    MPIDataCommunicator(MPIDataCommunicator const &rOther) = delete;

    ///@}

}; // Class MPIDataCommunicator

///@}

///@name Input and output
///@{

/// input stream function
inline std::istream &operator>>(std::istream &rIStream,
                                MPIDataCommunicator &rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
                                const MPIDataCommunicator &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_MPI_DATA_COMMUNICATOR_H_INCLUDED  defined
