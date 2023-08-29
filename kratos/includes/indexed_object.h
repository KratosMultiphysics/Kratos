//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//
//

#pragma once

// System includes
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"

namespace Kratos
{
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

/**
 * @class IndexedObject
 * @ingroup KratosCore
 * @brief This object defines an indexed object
 * @details It is used as base of any class which must be indexed
 * @author Riccardo Rossi
 */
class IndexedObject
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of IndexedObject
    KRATOS_CLASS_POINTER_DEFINITION(IndexedObject);

    /// The definition of the index type
    typedef std::size_t IndexType;

    /// The definition of the result_type
    typedef std::size_t result_type;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    explicit IndexedObject(IndexType NewId = 0) : mId(NewId) {}

    /// Destructor.
    virtual ~IndexedObject() {}

    /// Copy constructor.
    IndexedObject(IndexedObject const& rOther) : mId(rOther.mId) {}

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    IndexedObject& operator=(IndexedObject const& rOther)
    {
        mId = rOther.mId;

        return *this;
    }

    template<class TObjectType>
    IndexType operator()(TObjectType const& rThisObject) const
    {
        return rThisObject.Id();
    }

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

    virtual void SetId(IndexType NewId)
    {
        mId = NewId;
    }

    /// TODO: remove this function when removing data_file_io object.
    IndexType& DepricatedIdAccess()
    {
        return mId;
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
        buffer << "indexed object # "
               << mId;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }

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

    IndexType mId; /// The Id held by the object

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        rSerializer.save("Id",mId);
    }

    virtual void load(Serializer& rSerializer)
    {
        rSerializer.load("Id",mId);
    }

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class IndexedObject

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  IndexedObject& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const IndexedObject& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

