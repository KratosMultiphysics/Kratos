//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

#if !defined(KRATOS_NODAL_DATA_H_INCLUDED )
#define  KRATOS_NODAL_DATA_H_INCLUDED

// System includes
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/// Stores all data and dofs which are stored in each elements
/** This class is the container for nodal data storing:
 *  Id : The Id of the node 
*/
class NodalData
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NodalData
    KRATOS_CLASS_POINTER_DEFINITION(NodalData);

    using IndexType = std::size_t;

    ///@}
    ///@name Life Cycle
    ///@{

    NodalData(IndexType TheId);

    /// Destructor.
    virtual ~NodalData();

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    NodalData& operator=(NodalData const& rOther);

    ///@}
    ///@name Operations
    ///@{


    ///@}
    ///@name Access
    ///@{


    IndexType Id() const
    {
        return mId;
    }

    IndexType GetId() const
    {
        return mId;
    }

    void SetId(IndexType NewId)
    {
        mId = NewId;
    }

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const;


    ///@}
    ///@name Friends
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

        IndexType mId;

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

    /// Copy constructor.
    NodalData(NodalData const& rOther);


    ///@}

}; // Class NodalData

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                NodalData& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const NodalData& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_NODAL_DATA_H_INCLUDED  defined
