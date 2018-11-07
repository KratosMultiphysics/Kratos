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

#ifndef KRATOS_DATA_COMMUNICATOR_H_INCLUDED
#define KRATOS_DATA_COMMUNICATOR_H_INCLUDED

// System includes
#include <string>
#include <iostream>
//#ifdef KRATOS_USING_MPI
//#include <mpi.h>
//#endif

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_components.h"

namespace Kratos
{
///@addtogroup Kratos Core
///@{

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
  */
class DataCommunicator
{
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DataCommunicator
    KRATOS_CLASS_POINTER_DEFINITION(DataCommunicator);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DataCommunicator() {}

    /// Copy constructor.
    DataCommunicator(DataCommunicator const &rOther) {};

    /// Destructor.
    virtual ~DataCommunicator() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    virtual void Barrier() const {}

    virtual void SumAll(int& rValue) const {}

    virtual void SumAll(double& rValue) const {}

    virtual void SumAll(array_1d<double,3>& rValue) const {}

    virtual void MinAll(int& rValue) const {}

    virtual void MinAll(double& rValue) const {}

    virtual void MinAll(array_1d<double,3>& rValue) const {}

    virtual void MaxAll(int& rValue) const {}

    virtual void MaxAll(double& rValue) const {}

    virtual void MaxAll(array_1d<double,3>& rValue) const {}

    ///@}
    ///@name Inquiry
    ///@{

    virtual int Rank() const
    {
        return 0;
    }

    virtual int Size() const
    {
        return 1;
    }

    virtual bool IsDistributed() const
    {
        return false;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        PrintInfo(buffer);
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << "DataCommunicator";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const {}

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
    DataCommunicator &operator=(DataCommunicator const &rOther) = delete;

    ///@}

}; // Class DataCommunicator

template class KRATOS_API(KRATOS_CORE) KratosComponents<DataCommunicator >;

//void KRATOS_API(KRATOS_CORE) AddKratosComponent(std::string const& Name, DataCommunicator const& ThisComponent);

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream &operator>>(std::istream &rIStream,
                                DataCommunicator &rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
                                const DataCommunicator &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_DATA_COMMUNICATOR_H_INCLUDED  defined
