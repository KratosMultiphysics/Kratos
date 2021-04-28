//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt

#if !defined(KRATOS_INTEGRATION_FLAGS_H_INCLUDED)
#define  KRATOS_INTEGRATION_FLAGS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "containers/flags.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/// Flags to control the integration
class KRATOS_API(KRATOS_CORE) IntegrationFlags
    : public Flags
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MapperFlags
    KRATOS_CLASS_POINTER_DEFINITION(IntegrationFlags);

    /// Local Flags
    KRATOS_DEFINE_LOCAL_FLAG(DO_NOT_CREATE_TESSELLATION_ON_SLAVE);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Destructor.
    virtual ~IntegrationFlags() = default;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "IntegrationFlags" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "IntegrationFlags";
    }

    ///@}

}; // Class MapperFlags

///@}

}  // namespace Kratos.

#endif // KRATOS_INTEGRATION_FLAGS_H_INCLUDED  defined