//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Teschemacher
//


#if !defined(KRATOS_IGA_FLAGS_H_INCLUDED )
#define  KRATOS_IGA_FLAGS_H_INCLUDED

// System includes
#include "includes/define.h"
#include "containers/flags.h"
// External includes

// Project includes

namespace Kratos
{
///@name Kratos Classes
///@{

/// Flags needed used in the IGAApplication
/** This class implements the flags needed in the IGAApplication.
* Some of them are exposed to Python
*/
class IgaFlags
{
public:
    ///@name Type Definitions
    ///@{
    /// Local Flags
    KRATOS_DEFINE_LOCAL_FLAG(FIX_DISPLACEMENT_X);
    KRATOS_DEFINE_LOCAL_FLAG(FIX_DISPLACEMENT_Y);
    KRATOS_DEFINE_LOCAL_FLAG(FIX_DISPLACEMENT_Z);
    KRATOS_DEFINE_LOCAL_FLAG(FIX_ROTATION_X);
    KRATOS_DEFINE_LOCAL_FLAG(FIX_ROTATION_Y);
    KRATOS_DEFINE_LOCAL_FLAG(FIX_ROTATION_Z);
    ///@}

}; // Class IGAFlags

}  // namespace Kratos.

#endif // KRATOS_IGA_FLAGS_H_INCLUDED  defined