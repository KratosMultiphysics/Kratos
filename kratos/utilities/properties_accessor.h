//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//

#if !defined(KRATOS_PROPERTIES_ACCESSOR_H)
#define  KRATOS_PROPERTIES_ACCESSOR_H

// System includes

// External includes

// Project includes
#include "includes/properties.h"


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
class PropertiesAccessor
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PropertiesAccessor
    KRATOS_CLASS_POINTER_DEFINITION(PropertiesAccessor);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PropertiesAccessor(){}

    /// Destructor.
    virtual ~PropertiesAccessor(){}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    template<class TVariableType>
    typename TVariableType::Type& operator[](const TVariableType& rV,
                                             const Geometry::Pointer pGeometry,
                                             const DataValueContainer& rDataContainer,
                                             const ProcessInfo& rCurrentProcessInfo,
                                             const int GaussPointIndex=0) const
    {
        // Perform some operations depending on the configuration of the Accessor

        return GetValue(rV);
    }

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
        buffer << "PropertiesAccessor" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "PropertiesAccessor";}

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
    PropertiesAccessor& operator=(PropertiesAccessor const& rOther){}

    /// Copy constructor.
    PropertiesAccessor(PropertiesAccessor const& rOther){}


    ///@}

}; // Class PropertiesAccessor

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                PropertiesAccessor& rThis){}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const PropertiesAccessor& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_PROPERTIES_ACCESSOR_H  defined


