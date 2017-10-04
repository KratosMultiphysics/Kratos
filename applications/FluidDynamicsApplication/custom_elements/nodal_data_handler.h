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

#if !defined(KRATOS_NODAL_DATA_HANDLER_H_INCLUDED)
#define KRATOS_NODAL_DATA_HANDLER_H_INCLUDED

// System includes
#include <iostream>
#include <string>

// External includes

// Project includes
#include "includes/define.h"

// Application includes
#include "custom_elements/elemental_data_handler.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
 */
template <class TDataType, unsigned int TNumNodes, class TStorageType>
class NodalDataHandler : public DataHandler<TDataType, TStorageType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NodalDataHandler
    KRATOS_CLASS_POINTER_DEFINITION(NodalDataHandler);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    NodalDataHandler(const Variable<TDataType>& rVariable);

    /// Destructor.
    virtual ~NodalDataHandler();

    ///@}
    ///@name Operations
    ///@{

    virtual void Initialize(const Element& rElement, const ProcessInfo& rProcessInfo);

    ///@}
    ///@name Access
    ///@{

    virtual const TStorageType& Get() const; // NOTE: can this method be const? see embedded problems

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
        buffer << "NodalDataHandler for " << this->mrVariable.Name() << std::endl;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "NodalDataHandler for " << this->mrVariable.Name() << std::endl;
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

    TStorageType mValues;

    ///@}

private:
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    NodalDataHandler& operator=(NodalDataHandler const& rOther);

    /// Copy constructor.
    NodalDataHandler(NodalDataHandler const& rOther);

    ///@}

}; // Class NodalDataHandler

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template <class TDataType, unsigned int TNumNodes, class TStorageType>
inline std::istream& operator>>(std::istream& rIStream, NodalDataHandler<TDataType,TNumNodes,TStorageType>& rThis)
{
    return rIStream;
}

/// output stream function
template <class TDataType, unsigned int TNumNodes, class TStorageType>
inline std::ostream& operator<<(std::ostream& rOStream, const NodalDataHandler<TDataType,TNumNodes,TStorageType>& rThis)
{
    rThis.PrintData(rOStream);
    return rOStream;
}
///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_FILENAME_H_INCLUDED  defined
