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

#if !defined(KRATOS_MAPPER_AXIS_ALIGNED_BOUNDING_BOH_H)
#define  KRATOS_MAPPER_AXIS_ALIGNED_BOUNDING_BOH_H

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "mapper_bounding_box.h"


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
class MapperAxisAlignedBoundingBox : public MapperBoundingBox
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MapperAxisAlignedBoundingBox
    KRATOS_CLASS_POINTER_DEFINITION(MapperAxisAlignedBoundingBox);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MapperAxisAlignedBoundingBox();

    /// Destructor.
    virtual ~MapperAxisAlignedBoundingBox() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


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
    std::string Info() const override
    {
        return "MapperAxisAlignedBoundingBox";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {}

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {}


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
    // MapperAxisAlignedBoundingBox& operator=(MapperAxisAlignedBoundingBox const& rOther) {}

    /// Copy constructor.
    MapperAxisAlignedBoundingBox(MapperAxisAlignedBoundingBox const& rOther) {}


    ///@}

}; // Class MapperAxisAlignedBoundingBox

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MAPPER_AXIS_ALIGNED_BOUNDING_BOH_H  defined
