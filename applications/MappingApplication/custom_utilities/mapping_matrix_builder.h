
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//


#if !defined(KRATOS_MAPPING_MATRIX_BUILDER_H_INCLUDED )
#define  KRATOS_MAPPING_MATRIX_BUILDER_H_INCLUDED


// System includes


// External includes


// Project includes
#include "includes/define.h"


namespace Kratos
{
///@addtogroup MappingApplication
///@{

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
template<bool IsDistributed>
class KRATOS_API(MAPPING_APPLICATION) MappingMatrixBuilder
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MappingMatrixBuilder
    KRATOS_CLASS_POINTER_DEFINITION(MappingMatrixBuilder);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MappingMatrixBuilder();

    /// Destructor.
    virtual ~MappingMatrixBuilder() = default;

    /// Copy constructor.
    MappingMatrixBuilder(MappingMatrixBuilder const& rOther) = delete;

    /// Assignment operator.
    MappingMatrixBuilder& operator=(MappingMatrixBuilder const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{


    ///@}

private:
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}

}; // Class MappingMatrixBuilder

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MAPPING_MATRIX_BUILDER_H_INCLUDED  defined
