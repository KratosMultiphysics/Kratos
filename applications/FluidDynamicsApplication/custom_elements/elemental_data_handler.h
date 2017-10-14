//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#if !defined(KRATOS_ELEMENTAL_DATA_HANDLER_H_INCLUDED)
#define KRATOS_ELEMENTAL_DATA_HANDLER_H_INCLUDED

// System includes
#include <iostream>
#include <string>

// External includes
#include "boost/numeric/ublas/matrix_proxy.hpp"

// Project includes
#include "containers/variable.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/process_info.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
 */
template <class TDataType, class TStorageType>
class DataHandler
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DataHandler
    KRATOS_CLASS_POINTER_DEFINITION(DataHandler);

    typedef TDataType ValueType;

    typedef TStorageType StorageType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DataHandler(const Variable<TDataType>& rVariable):
    mrVariable(rVariable)
    {
    }

    /// Destructor.
    virtual ~DataHandler()
    {
    }

    ///@}
    ///@name Operations
    ///@{

    virtual void Initialize(const Element& rElement, const ProcessInfo& rProcessInfo)
    {
        KRATOS_ERROR << "Accessing base class DataHandler::Initialize() method." << std::endl;
    }

	virtual TDataType Interpolate(boost::numeric::ublas::matrix_row< Matrix >& rN, Element* pElement)
	{
		KRATOS_ERROR << "Accessing base class DataHandler::Interpolate() method." << std::endl;
	}
    
    ///@}
    ///@name Access
    ///@{

    virtual void Set(const Element& rElement, const ProcessInfo& rProcessInfo)
    {
        KRATOS_ERROR << "Accessing base class DataHandler::Set() method." << std::endl;
    }

    virtual TStorageType& Get()
    {
        KRATOS_ERROR << "Accessing base class DataHandler::Get() method." << std::endl;
    }

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
        buffer << "DataHandler for " << mrVariable.Name() << std::endl;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "DataHandler for " << mrVariable.Name() << std::endl;
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        PrintInfo(rOStream);
    }

    ///@}

protected:
    ///@name Protected member Variables
    ///@{

    const Variable<TDataType>& mrVariable;

    ///@}

private:
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    DataHandler& operator=(DataHandler const& rOther)
    {
    }

    /// Copy constructor.
    DataHandler(DataHandler const& rOther)
    {
    }

    ///@}

}; // Class DataHandler

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template< class TDataType, class TStorageType >
inline std::istream& operator>>(std::istream& rIStream, DataHandler<TDataType, TStorageType>& rThis)
{
    return rIStream;
}

/// output stream function
template< class TDataType, class TStorageType >
inline std::ostream& operator<<(std::ostream& rOStream, const DataHandler<TDataType, TStorageType>& rThis)
{
    rThis.PrintData(rOStream);
    return rOStream;
}
///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_FILENAME_H_INCLUDED  defined
