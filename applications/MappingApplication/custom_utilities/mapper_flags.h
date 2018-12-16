//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

#if !defined(KRATOS_MAPPER_FLAGS_H_INCLUDED)
#define  KRATOS_MAPPER_FLAGS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "containers/flags.h"


namespace Kratos
{
///@addtogroup MappingApplication
///@{

///@name Kratos Classes
///@{

/// Flags needed used in the MappingApplication
/** This class implements the flags needed in the MappingApplication.
* Some of them are exposed to Python
*/
class MapperFlags
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MapperFlags
    KRATOS_CLASS_POINTER_DEFINITION(MapperFlags);

    /// Local Flags
    KRATOS_DEFINE_LOCAL_FLAG( SWAP_SIGN );
    KRATOS_DEFINE_LOCAL_FLAG( ADD_VALUES );
    KRATOS_DEFINE_LOCAL_FLAG( REMESHED );
    KRATOS_DEFINE_LOCAL_FLAG( USE_TRANSPOSE );
    KRATOS_DEFINE_LOCAL_FLAG( ORIGIN_ONLY );
    KRATOS_DEFINE_LOCAL_FLAG( DESTINATION_ONLY );
    KRATOS_DEFINE_LOCAL_FLAG( TO_NON_HISTORICAL );
    KRATOS_DEFINE_LOCAL_FLAG( FROM_NON_HISTORICAL );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MapperFlags() = delete;

    /// Destructor.
    virtual ~MapperFlags() = default;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "MapperFlags" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MapperFlags";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

    ///@}

}; // Class MapperFlags

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MAPPER_FLAGS_H_INCLUDED  defined